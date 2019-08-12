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

 if (m_strButcherTable == "jab_gark") {

      Announce("Timestepping with JAB Gark method");

      iStages = 2;
      iQorder = 2;
      iPorder = -1;
      realtype gark_gamma = 1.0 - sqrt(2.0)/2.0;

      double * Aee = (realtype *) calloc( iStages*iStages, sizeof(realtype*) );
      double * be  = (realtype *) calloc( iStages, sizeof(realtype) );
      double * ce  = (realtype *) calloc( iStages, sizeof(realtype) );

      Aee[0] = 0.0;
      Aee[1] = 0.0;
      Aee[2] = 1.0 / (2.0 * gark_gamma);
      Aee[3] = 0.0;

      be[0] = 1.0 - gark_gamma;
      be[1] = gark_gamma;

      ce[0] = 0.0;
      ce[1] = 1.0 / gark_gamma;

      double * Aei = (realtype *) calloc( iStages*iStages, sizeof(realtype*) );
      Aei[0] = 0.0;
      Aei[1] = 0.0;
      Aei[2] = 1.0 / (2.0 * gark_gamma);
      Aei[3] = 0.0;

      double * Aii = (realtype *) calloc( iStages*iStages, sizeof(realtype*) );
      double * bi  = (realtype *) calloc( iStages, sizeof(realtype) );
      double * ci  = (realtype *) calloc( iStages, sizeof(realtype) );

      Aii[0] = gark_gamma;
      Aii[1] = 0.0;
      Aii[2] = 1.0 - gark_gamma;
      Aii[3] = gark_gamma;

      bi[0] = 1.0 - gark_gamma;
      bi[1] = gark_gamma;

      ci[0] = gark_gamma;
      ci[1] = 1.0;

      double * Aie = (realtype *) calloc( iStages*iStages, sizeof(realtype*) );
      Aie[0] = gark_gamma;
      Aie[1] = 0.0;
      Aie[2] = 1.0 - gark_gamma;
      Aie[3] = gark_gamma;

      ierr = IMEXGARKStepSetTables(arkode_mem, iStages, 
				   iQorder, iPorder,
				   ce, ci,
				   Aee, Aei,
				   Aie, Aii,
				   be, bi,
				   NULL, NULL);
      
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
