///////////////////////////////////////////////////////////////////////////////
///
///	\file    TimestepSchemeARKode.cpp
///	\author  David J. Gardner
///	\version February 4, 2016
///
///	<remarks>
///		Copyright 2000-2010 Paul Ullrich, Jorge Guerra, 
///             Daniel R. Reynolds, David J. Gardner
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

///////////////////////////////////////////////////////////////////////////////

// require SUNDIALS for compilation
#ifdef USE_SUNDIALS

// #define DEBUG_OUTPUT
// #define STATISTICS_OUTPUT

//#define DSS_INPUT
#define DSS_OUTPUT

#include "TimestepSchemeARKode.h"
#include "Model.h"
#include "Grid.h"
#include "HorizontalDynamics.h"
#include "VerticalDynamics.h"
#include "HorizontalDynamicsFEM.h"
#include "VerticalDynamicsFEM.h"
#include "Announce.h"

double tempdt;
int numsteps;
void * TimestepSchemeARKode::arkode_mem = NULL;

///////////////////////////////////////////////////////////////////////////////

TimestepSchemeARKode::TimestepSchemeARKode(
        Model & model,
	ARKodeCommandLineVariables & ARKodeVars
) :
        TimestepScheme(model),
	m_iNVectors(ARKodeVars.nvectors+2), // Add two extra temp vectors for applying hyperdiffusion
	m_iARKodeButcherTable(ARKodeVars.ARKodeButcherTable),
	m_strButcherTable(ARKodeVars.ButcherTable),
	m_strStepOut(ARKodeVars.StepOut),
	m_dRelTol(ARKodeVars.rtol),
	m_dAbsTol(ARKodeVars.atol),
	m_fFullyExplicit(ARKodeVars.FullyExplicit),
	m_fFullyImplicit(ARKodeVars.FullyImplicit),
	m_fDynamicStepSize(ARKodeVars.DynamicStepSize),
	m_iErrController(ARKodeVars.ErrController),
	m_dDynamicDeltaT(0.0),
	m_fAAFP(ARKodeVars.AAFP),
	m_iAAFPAccelVec(ARKodeVars.AAFPAccelVec),
	m_iNonlinIters(ARKodeVars.NonlinIters),
	m_iLinIters(ARKodeVars.LinIters),
	m_iPredictor(ARKodeVars.Predictor),
	m_dVAtol_vel(ARKodeVars.VAtol_vel),
	m_dVAtol_rho(ARKodeVars.VAtol_rho),
	m_dVAtol_theta(ARKodeVars.VAtol_theta),
	m_fWriteDiagnostics(ARKodeVars.WriteDiagnostics),
	m_fUsePreconditioning(ARKodeVars.UsePreconditioning),
	m_fColumnSolver(ARKodeVars.ColumnSolver),
	m_tOutT(ARKodeVars.OutputTime)
{
        // error flag
        int ierr = 0;

        // Allocate ARKode memory
        // arkode_mem = ARKodeCreate();
	// if (arkode_mem == NULL) _EXCEPTIONT("ERROR: ARKodeCreate returned NULL");

	// Set number of NVectors to use (default is 50)
	ierr = SetMaxTempestNVectorRegistryLength(m_iNVectors);
	if (ierr < 0) 
	  _EXCEPTIONT("ERROR: arkode_nvectors+2 > MAX_TEMPEST_NVECTORS");

#ifdef DSS_INPUT
	Announce("Apply DSS to RHS INPUT");
#endif
#ifdef DSS_OUTPUT
	Announce("Apply DSS to RHS OUTPUT");
#endif
}
					  
///////////////////////////////////////////////////////////////////////////////

void TimestepSchemeARKode::Initialize() {
  
  AnnounceStartBlock("Initializing ARKode");

  // error flag
  int ierr = 0;

  int pq = 0;

  numsteps = 0;

  // Get a copy of the grid
  Grid * pGrid = m_model.GetGrid();

  // Create a Tempest NVector for the model state (registry index 0)
  m_Y = N_VNew_Tempest(*pGrid, m_model);
  if (m_Y == NULL) _EXCEPTIONT("ERROR: N_VNew_Tempest returned NULL");

  int iY = NV_INDEX_TEMPEST(m_Y);
  if (iY != 0) _EXCEPTION1("ERROR: iY != 0, iY = %i",iY);

  // Reserve next two places (registry index 1 and 2) in nvector registry
  // to use as temporary states with hyperdiffusion
  int iRegistryIdx;

  iRegistryIdx = ReserveNextTempestNVectorRegistryIdx();
  if (iRegistryIdx < 0) _EXCEPTIONT("ERROR: NVector Registry is full");

  iRegistryIdx = ReserveNextTempestNVectorRegistryIdx();
  if (iRegistryIdx < 0) _EXCEPTIONT("ERROR: NVector Registry is full");

  // Get current simulation time in seconds
  Time timeCurrentT = m_model.GetCurrentTime();
  double dCurrentT  = timeCurrentT.GetSeconds();

  // Initialize ARKode
  if (m_fFullyExplicit) {
    Announce("Fully explicit not allowed");
  } else if (m_fFullyImplicit) {
    Announce("Fully implicit not allowed");
  } else {
    Announce("Running ARKode IMEX");
    arkode_mem = IMEXGARKStepCreate(ARKodeExplicitRHS, ARKodeImplicitRHS, dCurrentT, m_Y);
  }
  if (ierr < 0) _EXCEPTION1("ERROR: IMXGARKStepCreate, ierr = %i",ierr);

  // Set stop time 
  Time timeEndT = m_model.GetEndTime();
  double dEndT  = timeEndT.GetSeconds();

  dNextOut = m_tOutT.GetSeconds();

  ierr = IMEXGARKStepSetStopTime(arkode_mem, m_tNextOutT.GetSeconds());
  if (ierr < 0) _EXCEPTION1("ERROR: IMEXGARKStepSetStopTime, ierr = %i",ierr);  

  // Set function to post process arkode steps
  ierr = IMEXGARKStepSetPostprocessStepFn(arkode_mem, ARKodePostProcessStep);
  if (ierr < 0) _EXCEPTION1("ERROR: IMEXGARKStepSetPostprocessStepFn, ierr = %i",ierr);

  // Set adaptive or fixed time stepping  
  if (m_fDynamicStepSize) {

    // set dynamic timestepping flag to true
    m_model.SetDynamicTimestepping(m_fDynamicStepSize);

    ierr = IMEXGARKStepSetAdaptivityMethod(arkode_mem, m_iErrController, 1, pq, NULL);
    if (ierr < 0) _EXCEPTION1("ERROR: IMEXGARKStepSetAdaptivityMethod, ierr = %i", ierr);
    ierr = IMEXGARKStepSetInitStep(arkode_mem, 1.0);
    if (ierr < 0) _EXCEPTION1("ERROR: IMEXGARKStepSetInitStep, ierr = %i", ierr);

  } else {

    // Set fixed step size in seconds
    Time timeDeltaT = m_model.GetDeltaT();
    double dDeltaT  = timeDeltaT.GetSeconds();
   
    ierr = IMEXGARKStepSetFixedStep(arkode_mem, dDeltaT);   
    if (ierr < 0) _EXCEPTION1("ERROR: IMEXGARKStepSetFixedStep, ierr = %i",ierr);
  }
  
  // Set Butcher table
  if (m_strButcherTable != "") {  
    SetButcherTable();
  }  
  
  // Specify tolerances
  
  #if !defined(USE_COMPONENT_WISE_TOLERANCES)
    ierr = IMEXGARKStepSStolerances(arkode_mem, m_dRelTol, m_dAbsTol);
    if (ierr < 0) _EXCEPTION1("ERROR: IMEXGARKStepSStolerances, ierr = %i",ierr);

  #else
    m_T = N_VNew_Tempest(*pGrid, m_model);

    AssignComponentWiseTolerances();

    ierr = IMEXGARKStepSVtolerances(arkode_mem, m_dRelTol, m_T);
    if (ierr < 0) _EXCEPTION1("ERROR: IMEXGARKStepSVtolerances, ierr = %i",ierr);

    Announce("Component-wise tolerances will be assigned");
  
  #endif

  // Nonlinear Solver Settings
  if (!m_fFullyExplicit) {
    
    // Set predictor method
    ierr = IMEXGARKStepSetPredictorMethod(arkode_mem, m_iPredictor);
    if (ierr < 0) _EXCEPTION1("ERROR: IMEXGARKStepSetPredictorMethod, ierr = %i",ierr);

    // Newton iteration
    if (m_fColumnSolver) {

      /* Set custom solver functions into ARKode memory structure, 
         and set relevant parameters so that it is called appropriately */

      SUNMatrix A;
      SUNLinearSolver LS;

      A = SUNMatrix_Tempest();
      LS = SUNLinSol_Tempest(arkode_mem);
      ierr = IMEXGARKStepSetLinearSolver(arkode_mem, LS, A);
      if (ierr < 0) _EXCEPTION1("ERROR: IMEXGARKStepSetLinearSolver, ierr = %i",ierr);

      // set lin sys fn for custom linear solver
      ierr = IMEXGARKStepSetLinSysFn(arkode_mem, ARKodeLinSysFn);
      if (ierr < 0) _EXCEPTION1("ERROR: IMEXGARKStepSetLinSysFn, ierr = %i",ierr);

    } else {
      
      // We are using SPGMR
      int precflag = PREC_NONE;
      if (m_fUsePreconditioning)  precflag = PREC_RIGHT;

      // Linear Solver Settings
      SUNLinearSolver LS = SUNLinSol_SPGMR(m_Y, precflag, m_iLinIters);
      ierr = IMEXGARKStepSetLinearSolver(arkode_mem, LS, NULL);
      if (ierr < 0) _EXCEPTION1("ERROR: IMEXGARKStepSetLinearSolver, ierr = %i",ierr);

      if (m_fUsePreconditioning) {
        ARKLsPrecSolveFn psolve = ARKodePreconditionerSolve;
        ierr = IMEXGARKStepSetPreconditioner(arkode_mem, NULL, psolve);
        if (ierr < 0) _EXCEPTION1("ERROR: IMEXGARKStepSetPreconditioner, ierr = %i",ierr);
      }
    }  

    // if negative nonlinear iterations are specified, switch to linear-implicit mode
    if (m_iNonlinIters < 0) {
      
      ierr = IMEXGARKStepSetLinear(arkode_mem, 1);
      if (ierr < 0) _EXCEPTION1("ERROR: IMEXGARKStepSetLinear, ierr = %i",ierr);
      
      // otherwise, set the Max nonlinear solver iterations
    } else { 
      
      ierr = IMEXGARKStepSetMaxNonlinIters(arkode_mem, m_iNonlinIters);
      if (ierr < 0) _EXCEPTION1("ERROR: IMEXGARKStepSetMaxNonlinIters, ierr = %i",ierr);
    }
  }

  // Set diagnostics output file
  if (m_fWriteDiagnostics) {

    // Get process rank
    int iRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &iRank);

    // Root process writes diagnostic file
    if (iRank == 0) {
      FILE * pFile;
      pFile = fopen("ARKode_Diagnostics.txt","w");
      
      ierr = IMEXGARKStepSetDiagnostics(arkode_mem, pFile);
      if (ierr < 0) _EXCEPTION1("ERROR: IMEXGARKStepSetDiagnostics, ierr = %i",ierr);
    }
  }

  AnnounceEndBlock("Done");
  int iRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &iRank);
  if (iRank == 0) {
    int nnn = m_strStepOut.length();
    char char_array[nnn+1];
    strcpy(char_array, m_strStepOut.c_str());
    m_fStep_Profile = fopen(char_array,"w");
    if (m_fStep_Profile == NULL) {
      Announce("Error creating the output file for the timestep profile.\n");
    }
  }


#ifdef DEBUG_OUTPUT
  // Get process rank
//  int iRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &iRank);

  if (iRank == 0) {
    std::cout << std::endl
	      << " Proc "
	      << iRank
	      << " N_Vector Registry Length After Init: "
	      << GetTempestNVectorRegistryLength()
	      << std::endl;
  }
#endif

#define NO_NVECTOR_TESTING
#ifdef NVECTOR_TESTING  
  // Call NVector "test" routine on m_Y and then halt simulation
  N_VTest_Tempest(m_Y);
  _EXCEPTIONT("Halt: NVector Testing complete");
#endif

}

///////////////////////////////////////////////////////////////////////////////

void TimestepSchemeARKode::Step(
	bool fFirstStep,
	bool fLastStep,
	const Time & time,
	double dDeltaT
) {
#ifdef DEBUG_OUTPUT
  AnnounceStartBlock("ARKode Step");
#endif

  // error flag
  int ierr = 0;

  int pq = 0;

  // Adjust last time step size if using fixed step sizes
  if (!m_fDynamicStepSize && fLastStep) {
    ierr = IMEXGARKStepSetFixedStep(arkode_mem, dDeltaT);
    if (ierr < 0) _EXCEPTION1("ERROR: IMEXGARKStepSetFixedStep, ierr = %i",ierr);
  }
  // Get a copy of the grid
  Grid * pGrid = m_model.GetGrid();

  // Get a copy of the HorizontalDynamics
  HorizontalDynamicsFEM * pHorizontalDynamicsFEM 
    = dynamic_cast<HorizontalDynamicsFEM*>(m_model.GetHorizontalDynamics());

  // Get a copy of the VerticalDynamics
  VerticalDynamicsFEM * pVerticalDynamicsFEM 
    = dynamic_cast<VerticalDynamicsFEM*>(m_model.GetVerticalDynamics());

  // Current time in seconds
  double dCurrentT = time.GetSeconds();

  // Next time in seconds  
  // Time timeDeltaT = m_model.GetDeltaT();
  // Time timeNextT  = time;

  // timeDeltaT is possibly too large on the last step because the value that 
  // GetDeltaT returns is not adjusted if necessary for the last step,
  // stop time in arkode will adjust time step to match end time
  // timeNextT += timeDeltaT; 

  int stepmode = ARK_ONE_STEP;

  // max next time is where you need output from
  double dNextT = dNextOut;

  // if at the very end, cut it off
  Time tEndT = m_model.GetEndTime();
  double dEndT = tEndT.GetSeconds();

  // if we've gone past the final time
  if (dNextT > dEndT) {
    dNextT = dEndT;
  }
  ierr = IMEXGARKStepSetStopTime(arkode_mem, dNextT);
  if (ierr < 0) _EXCEPTION1("ERROR: IMEXGARKStepSetStopTime, ierr = %i",ierr);
  
  // ARKode timestep
  Announce("curr t = %f, next t = %f\n", dCurrentT, dNextT);
  ierr = IMEXGARKStepEvolve(arkode_mem, dNextT, m_Y, &dCurrentT, stepmode);
  if (ierr < 0) _EXCEPTION1("ERROR: IMEXGARKStepEvolve, ierr = %i",ierr);
  numsteps++;
  
  // Output for step
  Announce("\n numsteps is %d\n", numsteps);
  int iJRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &iJRank);
  if (dCurrentT == dNextT) {
    dNextOut += m_tOutT.GetSeconds();
  }

#ifdef STATISTICS_OUTPUT
//  if (fLastStep || ierr < 0) {
 if (true) {   
    int iRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &iRank);
    
    if (iRank == 0 && !m_fFullyExplicit) {
      long int nsteps, expsteps, accsteps, step_attempts, nfe_evals, nfi_evals, nlinsetups, 
        netfails, nniters, nncfails, npsolves, nliters, nlcfails, nfevalsLS;
//      realtype hinused, hlast, hcur, tcur;
//      ierr = ARKodeGetIntegratorStats(arkode_mem, &nsteps, &expsteps, &accsteps, &step_attempts, 
//                                      &nfe_evals, &nfi_evals, &nlinsetups, &netfails, &hinused, 
//                                      &hlast, &hcur, &tcur);
      ierr = IMEXGARKStepGetTimestepperStats(arkode_mem, &nsteps, &expsteps, &accsteps, &step_attempts,
					&nfe_evals, &nfi_evals, &nlinsetupts, &netfails);
      if (ierr < 0) _EXCEPTION1("ERROR: IMEXGARKStepGetIntegratorStats, ierr = %i",ierr);

      ierr = IMEXGARKStepGetNonlinSolvStats(arkode_mem, &nniters, &nncfails);
      if (ierr < 0) _EXCEPTION1("ERROR: IMEXGARKStepGetNonlinSolvStats, ierr = %i",ierr);
      
      if (m_fColumnSolver) {
        
        std::cout << std::endl
                  << "TimestepSchemeARKode::Step cumulative stats:\n"
                  << "  steps: " << nsteps << " (" << step_attempts << " attempted)\n"
                  << "  step size: " << hlast << " previous, " << hcur << " next\n"
                  << "  fevals: " << nfe_evals << " exp, " << nfi_evals << " imp\n" 
                  << "  nonlinear: " << nniters << " iters, " << nncfails << " failures\n"
                  << std::endl;
        
      } else {
        
        ierr = ARKSpilsGetNumPrecSolves(arkode_mem, &npsolves);
        if (ierr < 0) _EXCEPTION1("ERROR: ARKodeGetNumPrecSolves, ierr = %i",ierr);
        
        ierr = ARKSpilsGetNumLinIters(arkode_mem, &nliters);
        if (ierr < 0) _EXCEPTION1("ERROR: ARKSpilsGetNumLinIters, ierr = %i",ierr);
        
        ierr = ARKSpilsGetNumConvFails(arkode_mem, &nlcfails);
        if (ierr < 0) _EXCEPTION1("ERROR: ARKSpilsGetNumConvFails, ierr = %i",ierr);
        
        ierr = ARKSpilsGetNumRhsEvals(arkode_mem, &nfevalsLS);
        if (ierr < 0) _EXCEPTION1("ERROR: ARKSpilsGetNumRhsEvals, ierr = %i",ierr);
        
        std::cout << std::endl
                  << "TimestepSchemeARKode::Step cumulative stats:\n"
                  << "  steps: " << nsteps << " (" << step_attempts << " attempted)\n"
                  << "  step size: " << hlast << " previous, " << hcur << " next\n"
                  << "  fevals: " << nfe_evals << " exp, " << nfi_evals << " imp, " << nfevalsLS << " lin solve\n" 
                  << "  nonlinear: " << nniters << " iters, " << nncfails << " failures\n"
                  << "  linear: " << nliters << " iters, " << nlcfails << " fails, " << npsolves << " prec solves\n"
                  << std::endl;
        
      } // column solver
    } // root proc
  } // last step
#endif

  // check ARKode return flag
  if (ierr < 0) _EXCEPTION1("ERROR: IMEXGARKStep, ierr = %i",ierr);

  ierr = IMEXGARKStepGetLastStep(arkode_mem, &m_dDynamicDeltaT);
  tempdt = m_dDynamicDeltaT;
  if (ierr < 0) _EXCEPTION1("ERROR: IMEXGARKStepGetLastStep, ierr = %i", ierr);
  int iRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &iRank);
  if (iRank == 0) {
    fprintf(m_fStep_Profile, "%f %f\n", dCurrentT - m_dDynamicDeltaT, m_dDynamicDeltaT);
  }
  if (iRank == 0 && fLastStep) {
    fclose(m_fStep_Profile);
  }
  
#ifdef DEBUG_OUTPUT
  AnnounceEndBlock("Done");

  // Get process rank
  if (fFirstStep) {
    int iRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &iRank);
    
    if (iRank == 0) {
      std::cout << std::endl
		<< " Proc "
		<< iRank
		<< " N_Vector Registry Length After First Step: "
		<< GetTempestNVectorRegistryLength()
		<< std::endl;
    }
  }
#endif
}

///////////////////////////////////////////////////////////////////////////////

static int ARKodePostProcessStep(
	realtype time, 
	N_Vector Y, 
	void * user_data
) {
#ifdef DEBUG_OUTPUT
  AnnounceStartBlock("Post Process Step");
#endif
  
  // error flag
  int ierr;

  // Get idex of Tempest NVector in the grid
  int iY = NV_INDEX_TEMPEST(Y);

  // Get a copy of the grid
  Grid * pGrid = NV_GRID_TEMPEST(Y);

  // Get a copy of the model
  Model * pModel = NV_MODEL_TEMPEST(Y);

  // Get a copy of the HorizontalDynamics
  HorizontalDynamicsFEM * pHorizontalDynamicsFEM 
    = dynamic_cast<HorizontalDynamicsFEM*>(pModel->GetHorizontalDynamics());

  // Get a copy of the VerticalDynamics
  VerticalDynamicsFEM * pVerticalDynamicsFEM 
    = dynamic_cast<VerticalDynamicsFEM*>(pModel->GetVerticalDynamics());

  // Filter negative tracers
  pHorizontalDynamicsFEM->FilterNegativeTracers(iY);
  pVerticalDynamicsFEM->FilterNegativeTracers(iY);
  
  // Perform DSS (average values at shared nodes)
  pGrid->PostProcessSubstage(iY, DataType_State);
  pGrid->PostProcessSubstage(iY, DataType_Tracers);

  // Get last time step size <<< NEED TO ADJUST STEP SIZE WITH ADAPTIVE STEPPING
  Time timeT     = pModel->GetCurrentTime();  // model still has old time
  double dOldT   = timeT.GetSeconds();
//  double dDeltaT = time - dOldT;
  double dDeltaT = tempdt;  
 
  // Apply hyperdiffusion (initial, update, temp)
  pGrid->CopyData(iY, 2, DataType_State);
  pGrid->CopyData(iY, 2, DataType_Tracers);
  
  pHorizontalDynamicsFEM->StepAfterSubCycle(2, 1, iY, timeT, dDeltaT);
  
  pGrid->CopyData(1, iY, DataType_State);
  pGrid->CopyData(1, iY, DataType_Tracers);

#ifdef DEBUG_OUTPUT
  AnnounceEndBlock("Done");
#endif

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

static int ARKodeExplicitRHS(
	realtype time, 
	N_Vector Y, 
	N_Vector Ydot, 
	void * user_data
) {
//Announce("jab explicit rhs begin\n");
#ifdef DEBUG_OUTPUT
  AnnounceStartBlock("Explicit RHS");
#endif

  // model time (not used in RHS)
  Time timeT; //= Time(0,0,0,0,0); // only works for integers

  // index of input data in registry
  int iY = NV_INDEX_TEMPEST(Y);

  // index of output data in registry
  int iYdot = NV_INDEX_TEMPEST(Ydot);

  // Get a copy of the grid
  Grid * pGrid = NV_GRID_TEMPEST(Y);

  // Get a copy of the model
  Model * pModel = NV_MODEL_TEMPEST(Y);

  // Get a copy of the HorizontalDynamics
  HorizontalDynamicsFEM * pHorizontalDynamicsFEM 
    = dynamic_cast<HorizontalDynamicsFEM*>(pModel->GetHorizontalDynamics());

  // Get a copy of the VerticalDynamics
  VerticalDynamicsFEM * pVerticalDynamicsFEM 
    = dynamic_cast<VerticalDynamicsFEM*>(pModel->GetVerticalDynamics());

  // Filter negative tracers
  pHorizontalDynamicsFEM->FilterNegativeTracers(iY);
  pVerticalDynamicsFEM->FilterNegativeTracers(iY);

#ifdef DSS_INPUT
  // Perform DSS (average values at shared nodes)
  pGrid->PostProcessSubstage(iY, DataType_State);
  pGrid->PostProcessSubstage(iY, DataType_Tracers);
#endif

  // zero out iYdot
  pGrid->ZeroData(iYdot, DataType_State);
  pGrid->ZeroData(iYdot, DataType_Tracers);
  // Compute explicit RHS
  pHorizontalDynamicsFEM->StepExplicit(iY, iYdot, timeT, 1.0);
  pVerticalDynamicsFEM->StepExplicit(iY, iYdot, timeT, 1.0);
#ifdef DSS_OUTPUT
  // Perform DSS on RHS (average values at shared nodes)
  pGrid->PostProcessSubstage(iYdot, DataType_State);
  pGrid->PostProcessSubstage(iYdot, DataType_Tracers);
#endif
#ifdef DEBUG_OUTPUT
  // output ||fe||_max for sanity check
  //Announce("||fe|| = %g", N_VMaxNorm(Ydot));
  AnnounceEndBlock("Done");
#endif

//  Announce("jab explicit rhs finish\n"); 
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

static int ARKodeImplicitRHS(
	realtype time, 
	N_Vector Y, 
	N_Vector Ydot, 
	void * user_data
) {
 
#ifdef DEBUG_OUTPUT
  AnnounceStartBlock("Implicit RHS");
#endif

  // model time (not used in RHS)
  Time timeT; //= Time(0,0,0,0,0); // only works for integers

  // index of input data in registry
  int iY = NV_INDEX_TEMPEST(Y);

  // index of output data in registry
  int iYdot = NV_INDEX_TEMPEST(Ydot);

  // Get a copy of the grid
  Grid * pGrid = NV_GRID_TEMPEST(Y);

  // Get a copy of the model
  Model * pModel = NV_MODEL_TEMPEST(Y);

  // Get a copy of the HorizontalDynamics
  HorizontalDynamicsFEM * pHorizontalDynamicsFEM 
    = dynamic_cast<HorizontalDynamicsFEM*>(pModel->GetHorizontalDynamics());

  // Get a copy of the VerticalDynamics
  VerticalDynamicsFEM * pVerticalDynamicsFEM 
    = dynamic_cast<VerticalDynamicsFEM*>(pModel->GetVerticalDynamics());

  // Filter negative tracers
  pHorizontalDynamicsFEM->FilterNegativeTracers(iY);
  pVerticalDynamicsFEM->FilterNegativeTracers(iY);

#ifdef DSS_INPUT
  // Perform DSS (average values at shared nodes)
  pGrid->PostProcessSubstage(iY, DataType_State);
  pGrid->PostProcessSubstage(iY, DataType_Tracers);
#endif

  // zero out iYdot
  pGrid->ZeroData(iYdot, DataType_State);
  pGrid->ZeroData(iYdot, DataType_Tracers);

  // Compute implicit RHS
  pVerticalDynamicsFEM->StepImplicitTermsExplicitly(iY, iYdot, timeT, 1.0);

#ifdef DSS_OUTPUT
  // Perform DSS on RHS (average values at shared nodes)
  pGrid->PostProcessSubstage(iYdot, DataType_State);
  pGrid->PostProcessSubstage(iYdot, DataType_Tracers);
#endif


#ifdef DEBUG_OUTPUT
  // output ||fi||_max for sanity check
  //Announce("||fi|| = %g", N_VMaxNorm(Ydot));
  AnnounceEndBlock("Done");
#endif

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

static int ARKodeFullRHS(
	realtype time, 
	N_Vector Y, 
	N_Vector Ydot, 
	void * user_data
) {

#ifdef DEBUG_OUTPUT
  AnnounceStartBlock("Full RHS");

  // Get process rank
  int iRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &iRank);

  if (iRank == 0) {
    std::cout << std::endl
	      << " Proc "
	      << iRank
	      << " Y Vector Registry Index: "
	      << NV_INDEX_TEMPEST(Y)
	      << " Ydot Vector Registry Index: "
	      << NV_INDEX_TEMPEST(Ydot) 
	      << std::endl;
  }
#endif

  // model time (not used in RHS)
  Time timeT; //= Time(0,0,0,0,0); // only works for integers

  // index of input data in registry
  int iY = NV_INDEX_TEMPEST(Y);

  // index of output data in registry
  int iYdot = NV_INDEX_TEMPEST(Ydot);

  // Get a copy of the grid
  Grid * pGrid = NV_GRID_TEMPEST(Y);

  // Get a copy of the model
  Model * pModel = NV_MODEL_TEMPEST(Y);

  // Get a copy of the HorizontalDynamics
  HorizontalDynamicsFEM * pHorizontalDynamicsFEM 
    = dynamic_cast<HorizontalDynamicsFEM*>(pModel->GetHorizontalDynamics());

  // Get a copy of the VerticalDynamics
  VerticalDynamicsFEM * pVerticalDynamicsFEM 
    = dynamic_cast<VerticalDynamicsFEM*>(pModel->GetVerticalDynamics());

  // Filter negative tracers
  pHorizontalDynamicsFEM->FilterNegativeTracers(iY);
  pVerticalDynamicsFEM->FilterNegativeTracers(iY);

#ifdef DSS_INPUT
  // Perform DSS (average values at shared nodes)
  pGrid->PostProcessSubstage(iY, DataType_State);
  pGrid->PostProcessSubstage(iY, DataType_Tracers);
#endif
    
  // zero out iYdot
  pGrid->ZeroData(iYdot, DataType_State);
  pGrid->ZeroData(iYdot, DataType_Tracers);

  // Compute full RHS
  pHorizontalDynamicsFEM->StepExplicit(iY, iYdot, timeT, 1.0);
  pVerticalDynamicsFEM->StepExplicit(iY, iYdot, timeT, 1.0);

#ifdef DSS_OUTPUT
  // Perform DSS on RHS (average values at shared nodes)
  pGrid->PostProcessSubstage(iYdot, DataType_State);
  pGrid->PostProcessSubstage(iYdot, DataType_Tracers);
#endif

#ifdef DEBUG_OUTPUT
  AnnounceEndBlock("Done");
#endif

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

static int ARKodePreconditionerSolve(
	realtype time, 
	N_Vector Y, 
	N_Vector F, 
	N_Vector R,
	N_Vector Z,
	realtype gamma,
	realtype delta,
	int lr,
	void *user_data
//	N_Vector TMP
) {

#ifdef DEBUG_OUTPUT
  AnnounceStartBlock("Preconditioner Solve");
#endif

  // model time (not used in preconditioning)
  Time timeT;

  // index of various N_Vector arguments in registry
  int iY = NV_INDEX_TEMPEST(Y);
  int iF = NV_INDEX_TEMPEST(F);
  int iR = NV_INDEX_TEMPEST(R);
  int iZ = NV_INDEX_TEMPEST(Z);

  // Get a copy of the grid
  Grid * pGrid = NV_GRID_TEMPEST(Y);

  // Get a copy of the model
  Model * pModel = NV_MODEL_TEMPEST(Y);

  // Get a copy of the VerticalDynamics
  VerticalDynamicsFEM * pVerticalDynamicsFEM 
    = dynamic_cast<VerticalDynamicsFEM*>(pModel->GetVerticalDynamics());

  // Copy right-hand side into solution N_Vector
  N_VScale_Tempest(1.0, R, Z);

  // Call column-wise linear solver (iZ holds RHS on input, solution on output)
  pVerticalDynamicsFEM->SolveImplicit(iY, iZ, timeT, gamma);


  /*
  // test for preconditioner accuracy
  //   fill R with ones
  N_VConst_Tempest(1.0, R);
  //   perform solve using scaled time step of 1
  pVerticalDynamicsFEM->SolveImplicit(iY, iR, timeT, 1.0);
  //   perform finite-difference matrix-vector product
  N_Vector V = N_VClone_Tempest(R);
  double sigma = 1.e-8;
  N_VLinearSum_Tempest(1.0, Y, sigma, R, V);
  N_Vector FV = N_VClone_Tempest(R);
  ARKodeImplicitRHS(time, V, FV, user_data);
  N_VLinearSum_Tempest(1.0/sigma, FV, -1.0/sigma, F, Z);   // Z = (F(Y+sigma*R)-F(Y))/sigma
  N_VLinearSum_Tempest(1.0, R, -1.0, Z, V);                // V = R-Z = (I-Ji)*R = (I-Ji)*(I-Ji)^{-1}*R
  N_VConst_Tempest(1.0, R);
  N_VLinearSum_Tempest(1.0, R, -1.0, V, Z);                // Z = R-V = R - (I-Ji)*(I-Ji)^{-1}*R
  printf("Estimated Preconditioner error = %g\n",sqrt(N_VDotProd_Tempest(Z, Z)));
  N_VDestroy_Tempest(V);
  N_VDestroy_Tempest(FV);
  */

#ifdef DEBUG_OUTPUT
  AnnounceEndBlock("Done");
#endif

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

// constructor for matrix for custom solve
SUNMatrix SUNMatrix_Tempest()
{
  SUNMatrix M = SUNMatNewEmpty();
  if (M == NULL) return(NULL);
  M->ops->getid = SUNMatGetID_Tempest;
  M->ops->clone = SUNMatClone_Tempest;
  return M;
}

SUNMatrix SUNMatClone_Tempest(SUNMatrix M)
{
  SUNMatrix B = SUNMatrix_Tempest();
  return B;
}

SUNMatrix_ID SUNMatGetID_Tempest(SUNMatrix M)
{
  return SUNMATRIX_CUSTOM;
}

///////////////////////////////////////////////////////////////////////////////

static int ARKodeLinSysFn(
	realtype t, N_Vector y,
	N_Vector fy, SUNMatrix A,
	SUNMatrix M, booleantype jok,
	booleantype *jcur, realtype gamma,
	void *user_data, N_Vector tmp1,
	N_Vector tmp2, N_Vector tmp3
) {
  return(ARKLS_SUCCESS);
}

///////////////////////////////////////////////////////////////////////////////

SUNLinearSolver SUNLinSol_Tempest(void* ark_mem)
{
  SUNLinearSolver S = SUNLinSolNewEmpty();

  S->ops->gettype = ARKodeColumnLType;
  S->ops->solve = ARKodeColumnLSolve;
  S->ops->free = ARKodeColumnLFree;

  S->content = ark_mem;

  return S;
}

///////////////////////////////////////////////////////////////////////////////

int ARKodeColumnLSolve(
        SUNLinearSolver S,
	SUNMatrix A,
	N_Vector x,
	N_Vector b,
	realtype tol
) {
#ifdef DEBUG_OUTPUT
  AnnounceStartBlock("ARKodeColumnLSolve Start");
#endif

  int ierr;

  // model time (not used in preconditioning)
  Time timeT;

  double t_gamma;
  N_Vector ark_ycurr;

  void *ark_mem;
  ark_mem = S->content;
  ierr = IMEXGARKStepGetCurrentGamma(ark_mem, &t_gamma);
  if (ierr < 0) _EXCEPTION1("ERROR: IMEXGARKStepGetCurrentGamma, ierr = %i",ierr);
  ierr = IMEXGARKStepGetCurrentState(ark_mem, &ark_ycurr);
  if (ierr < 0) _EXCEPTION1("ERROR: IMEXGARKStepGetCurrentState, ierr = %i",ierr);

  // index of relevant N_Vector arguments in registry
  int iB = NV_INDEX_TEMPEST(b);
  int iY = NV_INDEX_TEMPEST(ark_ycurr);

  // Get a copy of the grid
  Grid * pGrid = NV_GRID_TEMPEST(b);

  // Get a copy of the model
  Model * pModel = NV_MODEL_TEMPEST(b);

  // Get a copy of the VerticalDynamics
  VerticalDynamicsFEM * pVerticalDynamicsFEM 
    = dynamic_cast<VerticalDynamicsFEM*>(pModel->GetVerticalDynamics());

  // Call column-wise linear solver (iB holds RHS on input, solution on output)
  pVerticalDynamicsFEM->SolveImplicit(iY, iB, timeT, t_gamma);

  N_VScale(1.0, b, x);

#ifdef DEBUG_OUTPUT
  AnnounceEndBlock("Done");
#endif
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

SUNLinearSolver_Type ARKodeColumnLType(SUNLinearSolver S)
{
  return SUNLINEARSOLVER_DIRECT;
}

///////////////////////////////////////////////////////////////////////////////

// this function does nothing since no persistent linear solver memory is used
int ARKodeColumnLFree(SUNLinearSolver S)
{

  free(S->content);
  S->content = NULL;

  free(S->ops);
  S->ops = NULL;

  SUNLinSolFreeEmpty(S);

#ifdef DEBUG_OUTPUT
  AnnounceStartBlock("ARKodeColumnLFree Start");
#endif

#ifdef DEBUG_OUTPUT
  AnnounceEndBlock("Done");
#endif

  return(0);
}

///////////////////////////////////////////////////////////////////////////////

void TimestepSchemeARKode::SetButcherTable()
{
  int ierr = 0; // error flag  

  int iStages;  // number of RK stages
  int iQorder;  // global order of accuracy for the method
  int iPorder;  // global order of accuracy for the embedding

  double * pci  = NULL; // implicit stage times
  double * pce  = NULL; // explicit stage times
  double * pAi  = NULL; // A matrix for the implicit method
  double * pAe  = NULL; // A matrix for the explicit method
  double * pbi  = NULL; // implicit b coefficient array
  double * pbe  = NULL; // explicit b coefficient array
  double * pb2i = NULL; // implicit b embedding array
  double * pb2e = NULL; // explicit b embedding array

 if (m_strButcherTable == "asg_222_sdirk") {

      Announce("Timestepping with ASG 222 SDIRK GARK method");

      iStages = 2;
      iQorder = 2;
      iPorder = -1;
      realtype g_g = 1.0 - sqrt(2.0)/2.0;

      double * Aee = (realtype *) calloc( iStages*iStages, sizeof(realtype*) );
      double * be  = (realtype *) calloc( iStages, sizeof(realtype) );
      double * ce  = (realtype *) calloc( iStages, sizeof(realtype) );

      Aee[0] = 0.0;
      Aee[1] = 0.0;
      Aee[2] = 1.0 / (2.0 * g_g);
      Aee[3] = 0.0;

      be[0] = 1.0 - g_g;
      be[1] = g_g;

      ce[0] = 0.0;
      ce[1] = 1.0 / g_g;

      double * Aei = (realtype *) calloc( iStages*iStages, sizeof(realtype*) );
      Aei[0] = 0.0;
      Aei[1] = 0.0;
      Aei[2] = 1.0 / (2.0 * g_g);
      Aei[3] = 0.0;

      double * Aii = (realtype *) calloc( iStages*iStages, sizeof(realtype*) );
      double * bi  = (realtype *) calloc( iStages, sizeof(realtype) );
      double * ci  = (realtype *) calloc( iStages, sizeof(realtype) );

      Aii[0] = g_g;
      Aii[1] = 0.0;
      Aii[2] = 1.0 - g_g;
      Aii[3] = g_g;

      bi[0] = 1.0 - g_g;
      bi[1] = g_g;

      ci[0] = g_g;
      ci[1] = 1.0;

      double * Aie = (realtype *) calloc( iStages*iStages, sizeof(realtype*) );
      Aie[0] = g_g;
      Aie[1] = 0.0;
      Aie[2] = 1.0 - g_g;
      Aie[3] = g_g;

      ierr = IMEXGARKStepSetTables(arkode_mem, iStages, 
				   iQorder, iPorder,
				   ce, ci,
				   Aee, Aei,
				   Aie, Aii,
				   be, bi,
				   NULL, NULL);
      
  }
  if (m_strButcherTable == "asg_332_sdirk") {

      Announce("Timestepping with ASG 332 SDIRK Gark method");

        iStages = 3;
        iQorder = 3;
        iPorder = -1;

	realtype ce2 = -0.1684172003560462642161441431010876;
	realtype ce3 =  0.8981606511689358909417299678385599;

        double * Aee = (realtype *) calloc( iStages*iStages, sizeof(realtype*) );
        double * be  = (realtype *) calloc( iStages, sizeof(realtype) );
        double * ce  = (realtype *) calloc( iStages, sizeof(realtype) );

        Aee[0] = 0.0;
        Aee[1] = 0.0;
        Aee[2] = 0.0;
        Aee[3] = ce2;
        Aee[4] = 0.0;
        Aee[5] = 0.0;
        Aee[6] = ce3*(3.0*ce2*ce2-3.0*ce2+ce3)/(ce2*(3.0*ce2-2.0));
        Aee[7] = ce3*(ce2-ce3)/(ce2*(3.0*ce2-2.0));
	Aee[8] = 0.0;

        be[0] = (6.0*ce2*ce3-3.0*ce2-3.0*ce3+2.0)/(6.0*ce2*ce3);
        be[1] = (2.0-3.0*ce3)/(6.0*ce2*ce2-6.0*ce2*ce3);
	be[2] = -(2.0-3.0*ce2)/(6.0*ce2-ce3-6.0*ce3*ce3);

	ce[0] = 0.0;
	ce[1] = ce2;
	ce[2] = ce3;

        double * Aei = (realtype *) calloc( iStages*iStages, sizeof(realtype*) );
        Aei[0] = 0.0;
        Aei[1] = 0.0;
	Aei[2] = 0.0;
        Aei[3] = ce2;
	Aei[4] = 0.0;
	Aei[5] = 0.0;
	Aei[6] = ce3*((3.0-sqrt(3.0))*ce2+(3.0+sqrt(3.0))*ce3-4.0)/(6.0*ce2-4.0);
	Aei[7] = (3.0+sqrt(3.0))*(ce2-ce3)*ce3/(6.0*ce2-4.0);
	Aei[8] = 0.0;

        double * Aii = (realtype *) calloc( iStages*iStages, sizeof(realtype*) );
        double * bi  = (realtype *) calloc( iStages, sizeof(realtype) );
        double * ci  = (realtype *) calloc( iStages, sizeof(realtype) );

        Aii[0] = (3.0+sqrt(3.0))/6.0;
        Aii[1] = 0.0;
        Aii[2] = 0.0;
        Aii[3] = -sqrt(3.0)/3.0;
	Aii[4] = (3.0+sqrt(3.0))/6.0;
	Aii[5] = 0.0;
	Aii[6] = 0.0;
	Aii[7] = 0.0;
	Aii[8] = 0.0;

        bi[0] = 0.5;
        bi[1] = 0.5;
	bi[2] = 0.0;

        ci[0] = (3.0+sqrt(3.0))/6.0;
        ci[1] = (3.0-sqrt(3.0))/6.0;
	ci[2] = 0.0;

        double * Aie = (realtype *) calloc( iStages*iStages, sizeof(realtype*) );
        Aie[0] = (3.0+sqrt(3.0))/6.0;
        Aie[1] = 0.0;
        Aie[2] = 0.0;
        Aie[3] = ((3.0-sqrt(3.0))*ce2-2.0)/(6.0*ce2);
	Aie[4] = 1.0/(3.0*ce2);
	Aie[5] = 0.0;
	Aie[6] = 0.0;
	Aie[7] = 0.0;
	Aie[8] = 0.0;

        ierr = IMEXGARKStepSetTables(arkode_mem, iStages, 
  				   iQorder, iPorder,
  				   ce, ci,
				   Aee, Aei,
				   Aie, Aii,
				   be, bi,
				   NULL, NULL);
      
  }
  else if (m_strButcherTable == "asg_333_sdirk") {
      
        Announce("Timestepping with ASG 333 SDIRK GARK method");
  
        iStages = 3;
        iQorder = 3;
        iPorder = -1;

	realtype g_g = 0.4358665215084589994160194511935568425292940929384256424857;
	realtype ce2 = -0.1684172003560462642161441431010876;
	realtype ce3 =  0.8981606511689358909417299678385599;
	realtype ci2 =  (6.0*g_g*g_g-9.0*g_g+2.0)/(6.0*g_g*g_g-12.0*g_g+3.0);
	realtype bi1 = (1.0-4.0*g_g)/(-12.0*g_g*g_g*g_g + 36.0*g_g*g_g - 24.0*g_g + 4.0);
	realtype bi2 = -3.0*(2.0*g_g*g_g-4.0*g_g+1.0)*(2.0*g_g*g_g-4.0*g_g+1.0)/(4.0*(3.0*g_g*g_g*g_g-9.0*g_g*g_g+6.0*g_g-1.0));
	realtype g_alpha = 2.0*(9.0*g_g*g_g*g_g*g_g-30.0*g_g*g_g*g_g+27.0*g_g*g_g-9.0*g_g+1.0)/(9.0*(2.0*g_g*g_g-4.0*g_g+1.0)*(2.0*g_g*g_g-4.0*g_g+1.0)*ce2);
	realtype g_beta = ce3*(4.0*(-3.0*g_g*g_g*g_g+9.0*g_g*g_g-6.0*g_g+1.0)-3.0*(4.0*g_g*g_g-5.0*g_g+1.0)*ce2+3.0*(6.0*g_g*g_g*g_g-14.0*g_g*g_g+7.0*g_g-1.0)*ce3)/(2.0*(3.0*g_g*g_g*g_g-9.0*g_g*g_g+6.0*g_g-1.0)*(3.0*ce2-2.0));

        double * Aee = (realtype *) calloc( iStages*iStages, sizeof(realtype*) );
        double * be  = (realtype *) calloc( iStages, sizeof(realtype) );
        double * ce  = (realtype *) calloc( iStages, sizeof(realtype) );

        Aee[0] = 0.0;
        Aee[1] = 0.0;
        Aee[2] = 0.0;
        Aee[3] = ce2;
        Aee[4] = 0.0;
        Aee[5] = 0.0;
        Aee[6] = ce3*(3.0*ce2*ce2-3.0*ce2+ce3)/(ce2*(3.0*ce2-2.0));
        Aee[7] = ce3*(ce2-ce3)/(ce2*(3.0*ce2-2.0));
        Aee[8] = 0.0;
  
        be[0] = (6.0*ce2*ce3-3.0*ce2-3.0*ce3+2.0)/(6.0*ce2*ce3);
        be[1] = (2.0-3.0*ce3)/(6.0*ce2*ce2-6.0*ce2*ce3);
	be[2] = -(2.0-3.0*ce2)/(6.0*ce2-ce3-6.0*ce3*ce3);

	ce[0] = 0.0;
	ce[1] = ce2;
	ce[2] = ce3;

        double * Aei = (realtype *) calloc( iStages*iStages, sizeof(realtype*) );
        Aei[0] = 0.0;
        Aei[1] = 0.0;
        Aei[2] = 0.0;
        Aei[3] = ce2;
        Aei[4] = 0.0;
        Aei[5] = 0.0;
        Aei[6] = g_beta;
        Aei[7] = ce3-g_beta;
        Aei[8] = 0.0;

        double * Aii = (realtype *) calloc( iStages*iStages, sizeof(realtype*) );
        double * bi  = (realtype *) calloc( iStages, sizeof(realtype) );
        double * ci  = (realtype *) calloc( iStages, sizeof(realtype) );

        Aii[0] = g_g;
        Aii[1] = 0.0;
        Aii[2] = 0.0;
        Aii[3] = ci2-g_g;
        Aii[4] = g_g;
        Aii[5] = 0.0;
        Aii[6] = bi1;
        Aii[7] = bi2;
        Aii[8] = g_g;

	bi[0] = bi1;
	bi[1] = bi2;
	bi[2] = g_g;

	ci[0] = g_g;
	ci[1] = ci2; 
	ci[2] = 1.0;

	double * Aie = (realtype *) calloc( iStages*iStages, sizeof(realtype*) );
	
        Aie[0] = g_g;
        Aie[1] = 0.0;
        Aie[2] = 0.0;
        Aie[3] = ci2 - g_alpha;
        Aie[4] = g_alpha;
        Aie[5] = 0.0;
	Aie[6] = be[0];
	Aie[7] = be[1];
	Aie[8] = be[2];

        ierr = IMEXGARKStepSetTables(arkode_mem, iStages, 
		  		     iQorder, iPorder,
				     ce, ci,
				     Aee, Aei,
				     Aie, Aii,
				     be, bi,
				     NULL, NULL);

  }
  else if (m_strButcherTable == "asg_333_esdirk") {
      
        Announce("Timestepping with ASG 333 ESDIRK GARK method");
  
        iStages = 3;
        iQorder = 3;
        iPorder = -1;

	realtype ce2 = -0.1684172003560462642161441431010876;
	realtype ce3 =  0.8981606511689358909417299678385599;
	realtype alpha1 = 0.0;
	realtype alpha2 = 0.0;

        double * Aee = (realtype *) calloc( iStages*iStages, sizeof(realtype*) );
        double * be  = (realtype *) calloc( iStages, sizeof(realtype) );
        double * ce  = (realtype *) calloc( iStages, sizeof(realtype) );

        Aee[0] = 0.0;
        Aee[1] = 0.0;
        Aee[2] = 0.0;
        Aee[3] = ce2;
        Aee[4] = 0.0;
        Aee[5] = 0.0;
        Aee[6] = ce3*(3.0*ce2*ce2-3.0*ce2+ce3)/(ce2*(3.0*ce2-2.0));
        Aee[7] = ce3*(ce2-ce3)/(ce2*(3.0*ce2-2.0));
        Aee[8] = 0.0;
  
        be[0] = (6.0*ce2*ce3-3.0*ce2-3.0*ce3+2.0)/(6.0*ce2*ce3);
        be[1] = (2.0-3.0*ce3)/(6.0*ce2*ce2-6.0*ce2*ce3);
	be[2] = -(2.0-3.0*ce2)/(6.0*ce2*ce3-6.0*ce3*ce3);

	ce[0] = 0.0;
	ce[1] = ce2;
	ce[2] = ce3;

        double * Aei = (realtype *) calloc( iStages*iStages, sizeof(realtype*) );
        Aei[0] = 0.0;
        Aei[1] = 0.0;
        Aei[2] = 0.0;
        Aei[3] = ce2;
        Aei[4] = 0.0;
        Aei[5] = 0.0;
        Aei[6] = ce3*(3.0*(sqrt(3.0)-2.0)*ce2-3.0*ce3-2.0*sqrt(3.0)+6.0)/((sqrt(3.0)-3.0)*(3.0*ce2-2.0));
        Aei[7] = 3.0*ce3*(ce3-ce2)/((sqrt(3.0)-3.0)*(3.0*ce2-2.0));
        Aei[8] = 0.0;

        double * Aii = (realtype *) calloc( iStages*iStages, sizeof(realtype*) );
        double * bi  = (realtype *) calloc( iStages, sizeof(realtype) );
        double * ci  = (realtype *) calloc( iStages, sizeof(realtype) );

        Aii[0] = 0.0;
        Aii[1] = 0.0;
        Aii[2] = 0.0;
        Aii[3] = (3.0-sqrt(3.0))/6.0;
        Aii[4] = (3.0-sqrt(3.0))/6.0;
        Aii[5] = 0.0;
        Aii[6] = (3.0-sqrt(3.0))/12.0;
        Aii[7] = (1.0+sqrt(3.0))/4.0;
        Aii[8] = (3.0-sqrt(3.0))/6.0;

	bi[0] = (3.0-sqrt(3.0))/12.0;
	bi[1] = (1.0+sqrt(3.0))/4.0;
	bi[2] = (3.0-sqrt(3.0))/6.0;

	ci[0] = 0.0;
	ci[1] = (sqrt(3.0)-1.0)/3.0; 
	ci[2] = 1.0;

	double * Aie = (realtype *) calloc( iStages*iStages, sizeof(realtype*) );
	
        Aie[0] = 0.0;
        Aie[1] = 0.0;
        Aie[2] = 0.0;
        Aie[3] = (sqrt(3.0)-1.0)/3.0 - alpha1;
        Aie[4] = alpha1;
        Aie[5] = 0.0;
	Aie[6] = 1.0-alpha2+(-3.0*(1.0+sqrt(3.0))*alpha1*ce2+2.0*(sqrt(3.0)-3.0)*alpha2*ce2+2.0)/(2.0*(sqrt(3.0)-3.0)*ce3);
	Aie[7] = alpha2;
	Aie[8] = (3.0*(1.0+sqrt(3.0))*alpha1*ce2-2.0*(sqrt(3.0)-3.0)*alpha2*ce2-2.0)/(2.0*(sqrt(3.0)-3.0)*ce3);
  }
  else if (m_strButcherTable == "asg_455_sdirk") {

        Announce("Timestepping with ASG 455 SDIRK GARK method");
  
        iStages = 5;
        iQorder = 4;
        iPorder = -1;

        double * Aee = (realtype *) calloc( iStages*iStages, sizeof(realtype*) );
        double * be  = (realtype *) calloc( iStages, sizeof(realtype) );
        double * ce  = (realtype *) calloc( iStages, sizeof(realtype) );
	double * de  = (realtype *) calloc( iStages, sizeof(realtype) );

        Aee[0] = 0.0;
        Aee[1] = 0.0;
        Aee[2] = 0.0;
	Aee[3] = 0.0;
	Aee[4] = 0.0;
	Aee[5] = 2.0;
	Aee[6] = 0.0;
	Aee[7] = 0.0;
	Aee[8] = 0.0;
	Aee[9] = 0.0;
	Aee[10] = 43.0/64.0;
	Aee[11] = 5.0/64.0;
	Aee[12] = 0.0;
	Aee[13] = 0.0;
	Aee[14] = 0.0;
	Aee[15] = 3.0/2.0;
	Aee[16] = 3.0/10.0;
	Aee[17] = -4.0/5.0;
	Aee[18] = 0.0;
	Aee[19] = 0.0;
	Aee[20] = 1.0/4.0;
	Aee[21] = 1.0/60.0;
	Aee[22] = 16.0/15.0;
	Aee[23] = -1.0/3.0;
	Aee[24] = 0.0;
        
	be[0] = 1.0/4.0;
        be[1] = 1.0/60.0;
	be[2] = 16.0/15.0;
	be[3] = -1.0/3.0; 
	be[4] = 0.0;

	ce[0] = 0.0;
	ce[1] = 2.0;
	ce[2] = 3.0/4.0;
	ce[3] = 1.0;
	ce[4] = 1.0;

	de[0] = 1.0/3.0;
	de[1] = -1.0/30.0;
	de[2] = 8.0/15.0;
	de[3] = 0.0;
	de[4] = 1.0/6.0;

        double * Aei = (realtype *) calloc( iStages*iStages, sizeof(realtype*) );

        Aei[0] = 0.0;
        Aei[1] = 0.0;
        Aei[2] = 0.0;
	Aei[3] = 0.0;
	Aei[4] = 0.0;
	Aei[5] = 2.0;
	Aei[6] = 0.0;
	Aei[7] = 0.0;
	Aei[8] = 0.0;
	Aei[9] = 0.0;
	Aei[10] = 75.0/92.0;
	Aei[11] = -3.0/46.0;
	Aei[12] = 0.0;
	Aei[13] = 0.0;
	Aei[14] = 0.0;
	Aei[15] = -463.0/414.0;
	Aei[16] = -264.0/1127.0;
	Aei[17] = 2075.0/882.0;
	Aei[18] = 0.0;
	Aei[19] = 0.0;
	Aei[20] = 0.0;
	Aei[21] = 0.0;
	Aei[22] = 10.0/3.0;
	Aei[23] = -7.0/3.0;
	Aei[24] = 0.0;

        double * Aii = (realtype *) calloc( iStages*iStages, sizeof(realtype*) );
        double * bi  = (realtype *) calloc( iStages, sizeof(realtype) );
        double * ci  = (realtype *) calloc( iStages, sizeof(realtype) );
	double * di  = (realtype *) calloc( iStages, sizeof(realtype) );

        Aii[0] = 5.0/8.0;
        Aii[1] = 0.0;
        Aii[2] = 0.0;
	Aii[3] = 0.0;
	Aii[4] = 0.0;
	Aii[5] = 115.0/32.0;
	Aii[6] = 5.0/8.0;
	Aii[7] = 0.0;
	Aii[8] = 0.0;
	Aii[9] = 0.0;
	Aii[10] = 448929.0/1909000.0;
	Aii[11] = -4851.0/477250.0;
	Aii[12] = 5.0/8.0;
	Aii[13] = 0.0;
	Aii[14] = 0.0;
	Aii[15] = 1437831.0/22554260.0;
	Aii[16] = -77127636.0/3039186535.0;
	Aii[17] = 71170425.0/211421672.0;
	Aii[18] = 5.0/8.0;
	Aii[19] = 0.0;
	Aii[20] = 31688.0/9315.0;
	Aii[21] = -1536.0/580405.0;
	Aii[22] = -20750.0/3969.0;
	Aii[23] = 49031.0/22248.0;
	Aii[24] = 5.0/8.0;
        
	bi[0] = 31688.0/9315.0;
        bi[1] = -1536.0/580405.0;
	bi[2] = -20750.0/3969.0;	
	bi[3] = 49031.0/22248.0; 
	bi[4] = 5.0/8.0;

	ci[0] = 5.0/8.0;
	ci[1] = 135.0/32.0;
	ci[2] = 17.0/20.0;
	ci[3] = 1.0;
	ci[4] = 1.0;

	de[0] = 11117363.0/64981440.0;
	de[1] = 62612719.0/2783622380.0;
	de[2] = 2063909125.0/609130368.0;
	de[3] = -49178093.0/19400256.0;
	de[4] = -655.0/13952.0;

	double * Aie = (realtype *) calloc( iStages*iStages, sizeof(realtype*) );
	
        Aie[0] = 5.0/8.0;
        Aie[1] = 0.0;
        Aie[2] = 0.0;
	Aie[3] = 0.0;
	Aie[4] = 0.0;
	Aie[5] = -4585.0/2048.0;
	Aie[6] = 13225.0/2048.0;
	Aie[7] = 0.0;
	Aie[8] = 0.0;
	Aie[9] = 0.0;
	Aie[10] = 1.0;
	Aie[11] = 13257.0/83000.0;
	Aie[12] = -25707.0/83000.0;
	Aie[13] = 0.0;
	Aie[14] = 0.0;
	Aie[15] = 13909039.0/9610076.0;
	Aie[16] = 17988381.0/48050380.0;
	Aie[17] = -11006374.0/12012595.0;
	Aie[18] = 4635.0/49031.0;
	Aie[19] = 0.0;
	Aie[20] = 1.0/4.0;
	Aie[21] = 1.0/60.0;
	Aie[22] = 16.0/15.0;
	Aie[23] = -1.0/3.0;
	Aie[24] = 0.0;

        ierr = IMEXGARKStepSetTables(arkode_mem, iStages, 
		  		     iQorder, iPorder,
				     ce, ci,
				     Aee, Aei,
				     Aie, Aii,
				     be, bi,
				     de, di);
  }
  else {
      _EXCEPTIONT("ERROR: Invalid IMEX Butcher table name");
  }   

  if (ierr < 0) _EXCEPTION1("ERROR: SetButcherTable, ierr = %i",ierr);
 
}

///////////////////////////////////////////////////////////////////////////////

void TimestepSchemeARKode::AssignComponentWiseTolerances() {

	// index of various N_Vector arguments in registry
	int iT = NV_INDEX_TEMPEST(m_T);

	// Get a copy of the grid
	Grid * pGrid = NV_GRID_TEMPEST(m_T);
	Announce("tolerances are %f, %f, %f\n", m_dVAtol_vel, m_dVAtol_rho, m_dVAtol_theta);

	// Tolerances are assigned in Tempest (objects Grid and GridPatch)
	pGrid->AssignComponentWiseTolerances(iT, m_dVAtol_vel, m_dVAtol_rho, m_dVAtol_theta);

} 

///////////////////////////////////////////////////////////////////////////////
#endif
