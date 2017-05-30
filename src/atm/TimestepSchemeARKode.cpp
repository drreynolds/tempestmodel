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

//#define DEBUG_OUTPUT
#define STATISTICS_OUTPUT

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
	m_dRelTol(ARKodeVars.rtol),
	m_dAbsTol(ARKodeVars.atol),
	m_fFullyExplicit(ARKodeVars.FullyExplicit),
	m_fFullyImplicit(ARKodeVars.FullyImplicit),
	m_fDynamicStepSize(ARKodeVars.DynamicStepSize),
	m_dDynamicDeltaT(0.0),
	m_fAAFP(ARKodeVars.AAFP),
	m_iAAFPAccelVec(ARKodeVars.AAFPAccelVec),
	m_iNonlinIters(ARKodeVars.NonlinIters),
	m_iLinIters(ARKodeVars.LinIters),
	m_iPredictor(ARKodeVars.Predictor),
	m_fWriteDiagnostics(ARKodeVars.WriteDiagnostics),
	m_fUsePreconditioning(ARKodeVars.UsePreconditioning),
	m_fColumnSolver(ARKodeVars.ColumnSolver)
{
        // error flag
        int ierr = 0;

        // Allocate ARKode memory
        arkode_mem = ARKodeCreate();
	if (arkode_mem == NULL) _EXCEPTIONT("ERROR: ARKodeCreate returned NULL");

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
    Announce("Running ARKode fully explicit");
    ierr = ARKodeInit(arkode_mem, ARKodeFullRHS, NULL, dCurrentT, m_Y);
  } else if (m_fFullyImplicit) {
    Announce("Running ARKode fully implicit");
    ierr = ARKodeInit(arkode_mem, NULL, ARKodeFullRHS, dCurrentT, m_Y);
  } else {
    Announce("Running ARKode IMEX");
    ierr = ARKodeInit(arkode_mem, ARKodeExplicitRHS, ARKodeImplicitRHS, dCurrentT, m_Y);
  }

  if (ierr < 0) _EXCEPTION1("ERROR: ARKodeInit, ierr = %i",ierr);

  // Set stop time 
  Time timeEndT = m_model.GetEndTime();
  double dEndT  = timeEndT.GetSeconds();

  ierr = ARKodeSetStopTime(arkode_mem, dEndT);
  if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetStopTime, ierr = %i",ierr);  

  // Set function to post process arkode steps
  ierr = ARKodeSetPostprocessStepFn(arkode_mem, ARKodePostProcessStep);
  if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetPostprocessStepFn, ierr = %i",ierr);

  // Set adaptive or fixed time stepping  
  if (m_fDynamicStepSize) {

    // set dynamic timestepping flag to true
    m_model.SetDynamicTimestepping(m_fDynamicStepSize);

  } else {

    // Set fixed step size in seconds
    Time timeDeltaT = m_model.GetDeltaT();
    double dDeltaT  = timeDeltaT.GetSeconds();
   
    ierr = ARKodeSetFixedStep(arkode_mem, dDeltaT);   
    if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetFixedStep, ierr = %i",ierr);
  }

  // Set Butcher table
  if (m_strButcherTable != "") {  
    SetButcherTable();
  }  
  
  // Specify tolerances
  
#if !defined(USE_COMPONENT_WISE_TOLERANCES)
  ierr = ARKodeSStolerances(arkode_mem, m_dRelTol, m_dAbsTol);
  if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSStolerances, ierr = %i",ierr);

#else
  m_T = N_VNew_Tempest(*pGrid,m_model);
  
  AssignComponentWiseTolerances();

  ierr = ARKodeSVtolerances(arkode_mem, m_dRelTol, m_T);
  if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSVtolerances, ierr = %i",ierr);

  Announce("Component-wise tolerances will be assigned");
#endif

  // Nonlinear Solver Settings
  if (!m_fFullyExplicit) {
    
    // Set predictor method
    ierr = ARKodeSetPredictorMethod(arkode_mem, m_iPredictor);
    if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetPredictorMethod, ierr = %i",ierr);

    if (m_fAAFP) { // Anderson accelerated fixed point solver
      
      ierr = ARKodeSetFixedPoint(arkode_mem, m_iAAFPAccelVec);
      if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetFixedPoint, ierr = %i",ierr);

    } else { // Newton iteration

      if (m_fColumnSolver) {

        /* Set custom solver functions into ARKode memory structure, 
           and set relevant parameters so that it is called appropriately */
        ARKodeMem ark_mem;
        //ark_mem = (ARKodeMem) arkode_mem;
        ark_mem = static_cast<ARKodeMem>(arkode_mem);
        ark_mem->ark_linit  = ARKodeColumnLInit;
        ark_mem->ark_lsetup = ARKodeColumnLSetup;
        ark_mem->ark_lsolve = ARKodeColumnLSolve;
        ark_mem->ark_lfree  = ARKodeColumnLFree;
        ark_mem->ark_lsolve_type = 4;
        ark_mem->ark_setupNonNull = FALSE;
        ark_mem->ark_lmem = NULL;

      } else {

        int precflag = PREC_NONE;
        if (m_fUsePreconditioning)  precflag = PREC_RIGHT;

        // Linear Solver Settings
        ierr = ARKSpgmr(arkode_mem, precflag, m_iLinIters);      
        if (ierr < 0) _EXCEPTION1("ERROR: ARKSpgmr, ierr = %i",ierr);

        // Use preconditioning if requested
        if (m_fUsePreconditioning) {
          ierr = ARKSpilsSetPreconditioner(arkode_mem, NULL, ARKodePreconditionerSolve);
          if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetPreconditioner, ierr = %i",ierr);
        }
      }  
    }

    // if negative nonlinear iterations are specified, switch to linear-implicit mode
    if (m_iNonlinIters < 0) {
      
      ierr = ARKodeSetLinear(arkode_mem, 1);
      if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetLinear, ierr = %i",ierr);
      
      // otherwise, set the Max nonlinear solver iterations
    } else { 
      
      ierr = ARKodeSetMaxNonlinIters(arkode_mem, m_iNonlinIters);
      if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetMaxNonlinIters, ierr = %i",ierr);
      
      //ierr = ARKodeSetNonlinConvCoef(arkode_mem, 0.001);
      //if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetNonlinConvCoef, ierr = %i",ierr);
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
      
      ierr = ARKodeSetDiagnostics(arkode_mem, pFile);
      if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetDiagnostics, ierr = %i",ierr);
    }
  }

  AnnounceEndBlock("Done");

#ifdef DEBUG_OUTPUT
  // Get process rank
  int iRank;
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

  // Adjust last time step size if using fixed step sizes
  if (!m_fDynamicStepSize && fLastStep) {
    ierr = ARKodeSetFixedStep(arkode_mem, dDeltaT);
    if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetFixedStep, ierr = %i",ierr);
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
  Time timeDeltaT = m_model.GetDeltaT();
  Time timeNextT  = time;

  // timeDeltaT is possibly too large on the last step because the value that 
  // GetDeltaT returns is not adjusted if necessary for the last step,
  // stop time in arkode will adjust time step to match end time
  timeNextT += timeDeltaT; 

  double dNextT = timeNextT.GetSeconds();

  // set up call to ARKode to evolve to dNextT, using either a single step or adaptivity
  int stepmode = ARK_ONE_STEP;
  if (m_fDynamicStepSize) {
    stepmode = ARK_NORMAL;

    Time timeEndT = m_model.GetEndTime();
    double dEndT  = timeEndT.GetSeconds();

    if (dNextT > dEndT) dNextT = dEndT;

    ierr = ARKodeSetStopTime(arkode_mem, dNextT);
    if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetStopTime, ierr = %i",ierr);
  } else {
    stepmode = ARK_ONE_STEP;
  }

  // ARKode timestep
  ierr = ARKode(arkode_mem, dNextT, m_Y, &dCurrentT, stepmode);
  if (ierr < 0) _EXCEPTION1("ERROR: ARKode, ierr = %i",ierr);

  // // With dynamic stepping, get the last step size to update model time
  // if (m_fDynamicStepSize) {
  //   ierr = ARKodeGetLastStep(arkode_mem, &m_dDynamicDeltaT);
  //   if (ierr < 0) _EXCEPTION1("ERROR: ARKodeGetLastStep, ierr = %i",ierr);
  // }

#ifdef STATISTICS_OUTPUT
  if (fLastStep) {
    
    int iRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &iRank);
    
    if (iRank == 0 && !m_fFullyExplicit) {
      long int nsteps, expsteps, accsteps, step_attempts, nfe_evals, nfi_evals, nlinsetups, 
        netfails, nniters, nncfails, npsolves, nliters, nlcfails, nfevalsLS;
      realtype hinused, hlast, hcur, tcur;
      ierr = ARKodeGetIntegratorStats(arkode_mem, &nsteps, &expsteps, &accsteps, &step_attempts, 
                                      &nfe_evals, &nfi_evals, &nlinsetups, &netfails, &hinused, 
                                      &hlast, &hcur, &tcur);
      if (ierr < 0) _EXCEPTION1("ERROR: ARKodeGetIntegratorStats, ierr = %i",ierr);

      ierr = ARKodeGetNonlinSolvStats(arkode_mem, &nniters, &nncfails);
      if (ierr < 0) _EXCEPTION1("ERROR: ARKodeGetNonlinSolvStats, ierr = %i",ierr);
      
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
  double dDeltaT = time - dOldT;
   
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
	void *user_data,
	N_Vector TMP
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

// this function does nothing since all relevant initialization done earlier
int ARKodeColumnLInit(ARKodeMem ark_mem) {

#ifdef DEBUG_OUTPUT
  AnnounceStartBlock("ARKodeColumnLInit Start");
#endif

#ifdef DEBUG_OUTPUT
  AnnounceEndBlock("Done");
#endif

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

// this function is a placeholder only (should not be called)
int ARKodeColumnLSetup(
        ARKodeMem ark_mem,
        int convfail,
        N_Vector ypred,
        N_Vector fpred,
        booleantype *jcurPtr,
        N_Vector vtemp1,
        N_Vector vtemp2,
        N_Vector vtemp3
) {

#ifdef DEBUG_OUTPUT
  AnnounceStartBlock("ARKodeColumnLSetup Start");
#endif

#ifdef DEBUG_OUTPUT
  AnnounceEndBlock("Done");
#endif

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int ARKodeColumnLSolve(
        ARKodeMem ark_mem,
        N_Vector b,
        N_Vector weight,
        N_Vector ycur,
        N_Vector fcur
) {

#ifdef DEBUG_OUTPUT
  AnnounceStartBlock("ARKodeColumnLSolve Start");
#endif

  // model time (not used in preconditioning)
  Time timeT;

  // index of relevant N_Vector arguments in registry
  int iB = NV_INDEX_TEMPEST(b);
  int iY = NV_INDEX_TEMPEST(ycur);

  // Get a copy of the grid
  Grid * pGrid = NV_GRID_TEMPEST(b);

  // Get a copy of the model
  Model * pModel = NV_MODEL_TEMPEST(b);

  // Get a copy of the VerticalDynamics
  VerticalDynamicsFEM * pVerticalDynamicsFEM 
    = dynamic_cast<VerticalDynamicsFEM*>(pModel->GetVerticalDynamics());

  // Call column-wise linear solver (iB holds RHS on input, solution on output)
  pVerticalDynamicsFEM->SolveImplicit(iY, iB, timeT, ark_mem->ark_gamma);

#ifdef DEBUG_OUTPUT
  AnnounceEndBlock("Done");
#endif

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

// this function does nothing since no persistent linear solver memory is used
int ARKodeColumnLFree(ARKodeMem ark_mem)
{

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

  if (m_fFullyExplicit) {

    // ============================================================================
    // Fully Explicit Methods
    // ============================================================================
    
    if (m_strButcherTable == "heun_euler_2_1_2") {
      ierr = ARKodeSetERKTableNum(arkode_mem, HEUN_EULER_2_1_2);

    } else if (m_strButcherTable == "bogacki_shampine_4_2_3") {
      ierr = ARKodeSetERKTableNum(arkode_mem, BOGACKI_SHAMPINE_4_2_3);

    } else if (m_strButcherTable == "ark324l2sa_erk_4_2_3") {
      ierr = ARKodeSetERKTableNum(arkode_mem, ARK324L2SA_ERK_4_2_3);

    } else if (m_strButcherTable == "zonneveld_5_3_4") {
      ierr = ARKodeSetERKTableNum(arkode_mem, ZONNEVELD_5_3_4);

    } else if (m_strButcherTable == "ark436l2sa_erk_6_3_4") {
      ierr = ARKodeSetERKTableNum(arkode_mem, ARK436L2SA_ERK_6_3_4);

    } else if (m_strButcherTable == "sayfy_aburub_6_3_4") {
      ierr = ARKodeSetERKTableNum(arkode_mem, SAYFY_ABURUB_6_3_4);

    } else if (m_strButcherTable == "cash_karp_6_4_5") {
      ierr = ARKodeSetERKTableNum(arkode_mem, CASH_KARP_6_4_5);

    } else if (m_strButcherTable == "fehlberg_6_4_5") {
      ierr = ARKodeSetERKTableNum(arkode_mem, FEHLBERG_6_4_5);

    } else if (m_strButcherTable == "dormand_prince_7_4_5") {
      ierr = ARKodeSetERKTableNum(arkode_mem, DORMAND_PRINCE_7_4_5);

    } else if (m_strButcherTable == "ark548l2sa_erk_8_4_5") {
      ierr = ARKodeSetERKTableNum(arkode_mem, ARK548L2SA_ERK_8_4_5);

    } else if (m_strButcherTable == "verner_8_5_6") {
      ierr = ARKodeSetERKTableNum(arkode_mem, VERNER_8_5_6);

    } else if (m_strButcherTable == "fehlberg_13_7_8") {
      ierr = ARKodeSetERKTableNum(arkode_mem, FEHLBERG_13_7_8);

    } else if (m_strButcherTable == "forward_euler") {

      // ------------------------------------------------------------------------
      // Forward Euler - 2 stages, 1st order
      // ------------------------------------------------------------------------
      Announce("Timestepping with Forward Euler");

      iStages = 2;
      iQorder = 1;
      iPorder = 0;
      
      pAe  = new double [iStages * iStages];
      pce  = new double [iStages];
      pbe  = new double [iStages];
      
      pAe[0]  = 0.0;  pAe[1]  = 0.0; 
      pAe[2]  = 1.0;  pAe[3]  = 0.0; 
      
      pce[0] = 0.0;
      pce[1] = 1.0;
      
      pbe[0] = 1.0;
      pbe[1] = 0.0;
            
      ierr = ARKodeSetERKTable(arkode_mem, iStages, iQorder, iPorder, 
			       pce, pAe, pbe, NULL);
      
      delete[] pce;
      delete[] pAe;
      delete[] pbe;
      
    } else if (m_strButcherTable == "kgu63") {

      // ------------------------------------------------------------------------
      // Kinnmark Gray Ullrich - 6 stages 3rd order
      // ------------------------------------------------------------------------
      Announce("Timestepping with KGU(6,3)");
      
      iStages = 6;
      iQorder = 3;
      iPorder = 0;
      
      pAe  = new double [iStages * iStages];
      pce  = new double [iStages];
      pbe  = new double [iStages];
      
      pAe[0]  = 0.0;  pAe[1]  = 0.0; pAe[2]  = 0.0;       pAe[3]  = 0.0;       pAe[4]  = 0.0;  pAe[5]  = 0.0;
      pAe[6]  = 0.2;  pAe[7]  = 0.0; pAe[8]  = 0.0;       pAe[9]  = 0.0;       pAe[10] = 0.0;  pAe[11] = 0.0;
      pAe[12] = 0.0;  pAe[13] = 0.2; pAe[14] = 0.0;       pAe[15] = 0.0;       pAe[16] = 0.0;  pAe[17] = 0.0;
      pAe[18] = 0.0;  pAe[19] = 0.0; pAe[20] = 1.0 / 3.0; pAe[21] = 0.0;       pAe[22] = 0.0;  pAe[23] = 0.0;
      pAe[24] = 0.0;  pAe[25] = 0.0; pAe[26] = 0.0;       pAe[27] = 2.0 / 3.0; pAe[28] = 0.0;  pAe[29] = 0.0;
      pAe[30] = 0.25; pAe[31] = 0.0; pAe[32] = 0.0;       pAe[33] = 0.0;       pAe[34] = 0.75; pAe[35] = 0.0;
      
      pce[0] = 0.0;
      pce[1] = 0.2;
      pce[2] = 0.2;
      pce[3] = 1.0 / 3.0;
      pce[4] = 2.0 / 3.0;
      pce[5] = 1.0;
      
      pbe[0] = 0.25;
      pbe[1] = 0.0;
      pbe[2] = 0.0;
      pbe[3] = 0.0;
      pbe[4] = 0.75;
      pbe[5] = 0.0;
             
      ierr = ARKodeSetERKTable(arkode_mem, iStages, iQorder, iPorder, 
			       pce, pAe, pbe, NULL);
      
      delete[] pce;
      delete[] pAe;
      delete[] pbe;

    } else if (m_strButcherTable == "ssprk54") {

      // ------------------------------------------------------------------------
      // Strong Stability Preserving ERK - 5 stages, 4th order
      // ------------------------------------------------------------------------
      Announce("Timestepping with SSPRK(5,4)");
      
      iStages = 5;
      iQorder = 4;
      iPorder = 0;
      
      pAe  = new double [iStages * iStages];
      pce  = new double [iStages];
      pbe  = new double [iStages];
      
      double alpha1 = 0.555629506348765;
      double alpha2 = 0.379898148511597;
      double alpha3 = 0.821920045606868;
      
      double beta1 = 0.517231671970585;
      double beta2 = 0.096059710526147;
      double beta3 = 0.386708617503259;
      
      pAe[0]  = 0.0;               pAe[1]  = 0.0;               pAe[2]  = 0.0;               pAe[3]  = 0.0;               pAe[4]  = 0.0;
      pAe[5]  = 0.391752226571890; pAe[6]  = 0.0;               pAe[7]  = 0.0;               pAe[8]  = 0.0;               pAe[9]  = 0.0;
      pAe[10] = alpha1 * pAe[5];   pAe[11] = 0.368410593050371; pAe[12] = 0.0;               pAe[13] = 0.0;               pAe[14] = 0.0;
      pAe[15] = alpha2 * pAe[10];  pAe[16] = alpha2 * pAe[11];  pAe[17] = 0.251891774271694; pAe[18] = 0.0;               pAe[19] = 0.0;
      pAe[20] = alpha3 * pAe[15];  pAe[21] = alpha3 * pAe[16];  pAe[22] = alpha3 * pAe[17];  pAe[23] = 0.544974750228521; pAe[24] = 0.0; 
      
      pce[0] = 0.0;
      pce[1] = pAe[5];
      pce[2] = pAe[10] + pAe[11];
      pce[3] = pAe[15] + pAe[16] + pAe[17];
      pce[4] = pAe[20] + pAe[21] + pAe[22] + pAe[23];
      
      pbe[0] = beta1 * pAe[10] + beta2 * pAe[15] + beta3 * pAe[20];
      pbe[1] = beta1 * pAe[11] + beta2 * pAe[16] + beta3 * pAe[21];
      pbe[2] = beta2 * pAe[17] + beta3 * pAe[22];
      pbe[3] = beta3 * pAe[23] + 0.063692468666290;
      pbe[4] = 0.226007483236906;
            
      ierr = ARKodeSetERKTable(arkode_mem, iStages, iQorder, iPorder, 
			       pce, pAe, pbe, NULL);      
      
      delete[] pce;
      delete[] pAe;
      delete[] pbe;

    } else {     
      _EXCEPTIONT("ERROR: Invalid explicit Butcher table name");
    }

    
  } else if (m_fFullyImplicit) {

    // ==========================================================================
    // Fully Implicit Methods
    // ==========================================================================
    
    if (m_strButcherTable == "sdirk_2_1_2") {
      ierr = ARKodeSetIRKTableNum(arkode_mem, SDIRK_2_1_2);
      
    } else if (m_strButcherTable == "billington_3_3_2") {
      ierr = ARKodeSetIRKTableNum(arkode_mem, BILLINGTON_3_3_2);
      
    } else if (m_strButcherTable == "trbdf2_3_3_2") {
      ierr = ARKodeSetIRKTableNum(arkode_mem, TRBDF2_3_3_2);

    } else if (m_strButcherTable == "kvaerno_4_2_3") {
      ierr = ARKodeSetIRKTableNum(arkode_mem, KVAERNO_4_2_3);

    } else if (m_strButcherTable == "ark324l2sa_dirk_4_2_3") {
      ierr = ARKodeSetIRKTableNum(arkode_mem, ARK324L2SA_DIRK_4_2_3);

    } else if (m_strButcherTable == "cash_5_2_4") {
      ierr = ARKodeSetIRKTableNum(arkode_mem, CASH_5_2_4);

    } else if (m_strButcherTable == "cash_5_3_4") {
      ierr = ARKodeSetIRKTableNum(arkode_mem, CASH_5_3_4);

    } else if (m_strButcherTable == "sdirk_5_3_4") {
      ierr = ARKodeSetIRKTableNum(arkode_mem, SDIRK_5_3_4);

    } else if (m_strButcherTable == "kvaerno_5_3_4") {
      ierr = ARKodeSetIRKTableNum(arkode_mem, KVAERNO_5_3_4);

    } else if (m_strButcherTable == "ark436l2sa_dirk_6_3_4") {
      ierr = ARKodeSetIRKTableNum(arkode_mem, ARK436L2SA_DIRK_6_3_4);

    } else if (m_strButcherTable == "kvaerno_7_4_5") {
      ierr = ARKodeSetIRKTableNum(arkode_mem, KVAERNO_7_4_5);

    } else if (m_strButcherTable == "ark548l2sa_dirk_8_4_5") {
      ierr = ARKodeSetIRKTableNum(arkode_mem, ARK548L2SA_DIRK_8_4_5);

    } else {     
      _EXCEPTIONT("ERROR: Invalid implicit Butcher table name");
    }

    // ierr = ARKodeSetIRKTable(arkode_mem, iStages, iQorder, iPorder, pci, pAi, pbi, pb2i)
    // if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetIRKTable, ierr = %i",ierr);

  } else {

    // ==========================================================================
    // IMEX Methods
    // ==========================================================================

    if (m_strButcherTable == "ark324l2sa_erk_4_2_3" ||
	m_strButcherTable == "ark324l2sa_dirk_4_2_3" ) {     
      ierr = ARKodeSetARKTableNum(arkode_mem, 
				  ARK324L2SA_DIRK_4_2_3, ARK324L2SA_ERK_4_2_3);

    } else if (m_strButcherTable == "ark436l2sa_erk_6_3_4" ||
	       m_strButcherTable == "ark436l2sa_dirk_6_3_4" ) {
      ierr = ARKodeSetARKTableNum(arkode_mem, 
				  ARK436L2SA_DIRK_6_3_4, ARK436L2SA_ERK_6_3_4);

    } else if (m_strButcherTable == "ark548l2sa_erk_8_4_5" ||
	       m_strButcherTable == "ark548l2sa_dirk_8_4_5" ) {
      ierr = ARKodeSetARKTableNum(arkode_mem, 
				  ARK548L2SA_DIRK_8_4_5, ARK548L2SA_ERK_8_4_5);

    } else if (m_strButcherTable == "ars233") {      

      // ------------------------------------------------------------------------
      // ARS233 - 2 implicit stages, 3 explicit stages, 3rd order
      //
      // Ascher, Ruuth, and Spiteri, Implicit-explicit Runge-Kutta methods for 
      // time-dependent partial differential equations, 1997. (section 2.4)
      //
      // Note: Weller, Lock, and Wood, Runge-Kutta IMEX schemes for the 
      // Horizontally Explicit/vertically Implicit (HEVI) soltuion of wave 
      // equations, 2013 uses different table for ARS233 than the Ascher, Ruuth,
      // and Spiteri paper. 
      // ------------------------------------------------------------------------
      Announce("Timestepping with ARS233");
      
      iStages = 3;
      iQorder = 3;
      iPorder = 0;
      
      pci  = new double [iStages];
      pce  = new double [iStages];
      pAi  = new double [iStages * iStages];
      pAe  = new double [iStages * iStages];
      pbi  = new double [iStages];
      pbe  = new double [iStages];
      
      double gamma = (3.0 + std::sqrt(3.0)) / 6.0;
      
      // Implicit table
      pci[0] = 0.0;
      pci[1] = gamma;
      pci[2] = 1.0 - gamma;
      
      pAi[0] = 0.0; pAi[1] = 0.0;               pAi[2] = 0.0;
      pAi[3] = 0.0; pAi[4] = gamma;             pAi[5] = 0.0;
      pAi[6] = 0.0; pAi[7] = 1.0 - 2.0 * gamma; pAi[8] = gamma;

      pbi[0] = 0.0;
      pbi[1] = 0.5;
      pbi[2] = 0.5;     

      // Explicit table
      pce[0] = 0.0;
      pce[1] = gamma;
      pce[2] = 1.0 - gamma;

      pAe[0] = 0.0;         pAe[1] = 0.0;               pAe[2] = 0.0;
      pAe[3] = gamma;       pAe[4] = 0.0;               pAe[5] = 0.0;
      pAe[6] = gamma - 1.0; pAe[7] = 2.0*(1.0 - gamma); pAe[8] = 0.0;
      
      pbe[0] = 0.0;
      pbe[1] = 0.5;
      pbe[2] = 0.5;
            
      // arkode memory, stages, order, emdedding order, ci, ce, Ai, Ae, bi, be, b2i, b2e
      ierr = ARKodeSetARKTables(arkode_mem, iStages, iQorder, iPorder, 
				pci, pce, pAi, pAe, pbi, pbe, NULL, NULL);

      delete[] pci;     
      delete[] pce;
      delete[] pAi;
      delete[] pAe;
      delete[] pbi;
      delete[] pbe;

    } else if (m_strButcherTable == "ars232") {      

      // ------------------------------------------------------------------------
      // ARS232 - 2 implicit stages, 3 explicit stages, 2nd order
      //          L-stable
      //
      // Ascher, Ruuth, and Spiteri, Implicit-explicit Runge-Kutta methods for 
      // time-dependent partial differential equations, 1997. (section 2.5)
      // ------------------------------------------------------------------------
      Announce("Timestepping with ARS232");
      
      iStages = 3;
      iQorder = 2;
      iPorder = 0;
      
      pci  = new double [iStages];
      pce  = new double [iStages];
      pAi  = new double [iStages * iStages];
      pAe  = new double [iStages * iStages];
      pbi  = new double [iStages];
      pbe  = new double [iStages];
      
      double gamma = 1.0 - 1.0/std::sqrt(2.0); // (2-sqrt(2))/2
      double delta = -2.0 * std::sqrt(2.0) / 3.0;
      
      // Implicit table
      pci[0] = 0.0;
      pci[1] = gamma;
      pci[2] = 1.0;
      
      pAi[0] = 0.0; pAi[1] = 0.0;         pAi[2] = 0.0;
      pAi[3] = 0.0; pAi[4] = gamma;       pAi[5] = 0.0;
      pAi[6] = 0.0; pAi[7] = 1.0 - gamma; pAi[8] = gamma;
      
      pbi[0] = 0.0;
      pbi[1] = 1.0 - gamma;
      pbi[2] = gamma;

      // Explicit table
      pce[0] = 0.0;
      pce[1] = gamma;
      pce[2] = 1.0;

      pAe[0] = 0.0;   pAe[1] = 0.0;         pAe[2] = 0.0;
      pAe[3] = gamma; pAe[4] = 0.0;         pAe[5] = 0.0;
      pAe[6] = delta; pAe[7] = 1.0 - delta; pAe[8] = 0.0;
      
      pbe[0] = 0.0;
      pbe[1] = 1.0 - gamma;
      pbe[2] = gamma;
            
      // arkode memory, stages, order, emdedding order, ci, ce, Ai, Ae, bi, be, b2i, b2e
      ierr = ARKodeSetARKTables(arkode_mem, iStages, iQorder, iPorder, 
				pci, pce, pAi, pAe, pbi, pbe, NULL, NULL);
      
      delete[] pci;
      delete[] pce;
      delete[] pAi;
      delete[] pAe;
      delete[] pbi;
      delete[] pbe;

    } else if (m_strButcherTable == "ars222") {      

      // ------------------------------------------------------------------------
      // ARS222 - 2 implicit stages, 2 explicit stages, 2nd order
      //          L-stable
      //
      // Ascher, Ruuth, and Spiteri, Implicit-explicit Runge-Kutta methods for 
      // time-dependent partial differential equations, 1997. (section 2.6)
      // ------------------------------------------------------------------------
      Announce("Timestepping with ARS222");
      
      iStages = 3;
      iQorder = 2;
      iPorder = 0;
      
      pci  = new double [iStages];
      pce  = new double [iStages];
      pAi  = new double [iStages * iStages];
      pAe  = new double [iStages * iStages];
      pbi  = new double [iStages];
      pbe  = new double [iStages];
      
      double gamma = 1.0 - 1.0/std::sqrt(2.0); // (2-sqrt(2))/2
      double delta = 1.0 - 1.0/(2.0*gamma);
      
      // Implicit table
      pci[0] = 0.0;
      pci[1] = gamma;
      pci[2] = 1.0;
      
      pAi[0] = 0.0; pAi[1] = 0.0;         pAi[2] = 0.0;
      pAi[3] = 0.0; pAi[4] = gamma;       pAi[5] = 0.0;
      pAi[6] = 0.0; pAi[7] = 1.0 - gamma; pAi[8] = gamma;
      
      pbi[0] = 0.0;
      pbi[1] = 1.0 - gamma;
      pbi[2] = gamma;

      // Explicit table
      pce[0] = 0.0;
      pce[1] = gamma;
      pce[2] = 1.0;

      pAe[0] = 0.0;   pAe[1] = 0.0;         pAe[2] = 0.0;
      pAe[3] = gamma; pAe[4] = 0.0;         pAe[5] = 0.0;
      pAe[6] = delta; pAe[7] = 1.0 - delta; pAe[8] = 0.0;
      
      pbe[0] = delta;
      pbe[1] = 1.0 - delta;
      pbe[2] = 0.0;
            
      // arkode memory, stages, order, emdedding order, ci, ce, Ai, Ae, bi, be, b2i, b2e
      ierr = ARKodeSetARKTables(arkode_mem, iStages, iQorder, iPorder, 
				pci, pce, pAi, pAe, pbi, pbe, NULL, NULL);
      
      delete[] pci;
      delete[] pce;
      delete[] pAi;
      delete[] pAe;
      delete[] pbi;
      delete[] pbe;

    } else if (m_strButcherTable == "ars343") {      

      // ------------------------------------------------------------------------
      // ARS343 - 3 implicit stages, 4 explicit stages, 3rd order
      //          L-stable
      //
      // Ascher, Ruuth, and Spiteri, Implicit-explicit Runge-Kutta methods for 
      // time-dependent partial differential equations, 1997. (section 2.7)
      // ------------------------------------------------------------------------
      Announce("Timestepping with ARS343");
      
      iStages = 4;
      iQorder = 3;
      iPorder = 0;
      
      pci  = new double [iStages];
      pce  = new double [iStages];
      pAi  = new double [iStages * iStages];
      pAe  = new double [iStages * iStages];
      pbi  = new double [iStages];
      pbe  = new double [iStages];
     
      double gamma  = 0.4358665215084590;
      double gamma2 = gamma * gamma;

      double b1 = -1.5 * gamma2 + 4.0 * gamma - 0.25;
      double b2 =  1.5 * gamma2 - 5.0 * gamma + 1.25;

      double a42 = 0.5529291480359398; // double check value 
      double a43 = 0.5529291480359398;

      double a31 = (1.0 - 4.5 * gamma + 1.5 * gamma2) * a42 
	+ (2.75 - 10.5 * gamma + 3.75 * gamma2) * a43 
	- 3.5 + 13 * gamma - 4.5 * gamma2;

      double a32 = (-1.0 + 4.5 * gamma - 1.5 * gamma2) * a42
	+ (-2.75 + 10.5 * gamma - 3.75 * gamma2) * a43
	+ 4.0 - 12.5 * gamma + 4.5 * gamma2;

      double a41 = 1.0 - a42 - a43;

      // double gamma = 0.4358665215;
      // double b1    = 1.208496649;
      // double b2    = -0.644363171;
      // double a31   = 0.3212788860;
      // double a32   = 0.3966543747;
      // double a41   = -0.105858296;
      // double a42   = 0.5529291479;
      // double a43   = 0.5529291479;
      
      // Implicit table
      pci[0] = 0.0;
      pci[1] = gamma;
      pci[2] = (1.0+gamma)/2.0;
      pci[3] = 1.0;
      
      pAi[0]  = 0.0;  pAi[1]  = 0.0;            pAi[2]  = 0.0;    pAi[3]  = 0.0; 
      pAi[4]  = 0.0;  pAi[5]  = gamma;          pAi[6]  = 0.0;    pAi[7]  = 0.0;
      pAi[8]  = 0.0;  pAi[9]  = (1-gamma)/2.0;  pAi[10] = gamma;  pAi[11] = 0.0; 
      pAi[12] = 0.0;  pAi[13] = b1;             pAi[14] = b2;     pAi[15] = gamma; 

      pbi[0] = 0.0;
      pbi[1] = b1;
      pbi[2] = b2;
      pbi[3] = gamma;

      // Explicit table
      pce[0] = 0.0;
      pce[1] = gamma;
      pce[2] = (1.0+gamma)/2.0;
      pce[3] = 1.0;
      
      pAe[0]  = 0.0;    pAe[1]  = 0.0;  pAe[2]  = 0.0;  pAe[3]  = 0.0; 
      pAe[4]  = gamma;  pAe[5]  = 0.0;  pAe[6]  = 0.0;  pAe[7]  = 0.0;
      pAe[8]  = a31;    pAe[9]  = a32;  pAe[10] = 0.0;  pAe[11] = 0.0; 
      pAe[12] = a41;    pAe[13] = a42;  pAe[14] = a43;  pAe[15] = 0.0; 

      pbe[0] = 0.0;     
      pbe[1] = b1;
      pbe[2] = b2;
      pbe[3] = gamma;
            
      // arkode memory, stages, order, emdedding order, ci, ce, Ai, Ae, bi, be, b2i, b2e
      ierr = ARKodeSetARKTables(arkode_mem, iStages, iQorder, iPorder, 
				pci, pce, pAi, pAe, pbi, pbe, NULL, NULL);

      delete[] pci;     
      delete[] pce;
      delete[] pAi;
      delete[] pAe;
      delete[] pbi;
      delete[] pbe;

    } else if (m_strButcherTable == "ars443") {      

      // ------------------------------------------------------------------------
      // ARS443 - 4 implicit stages, 4 explicit stages, 3rd order
      //          L-stable
      //
      // Ascher, Ruuth, and Spiteri, Implicit-explicit Runge-Kutta methods for 
      // time-dependent partial differential equations, 1997. (section 2.8)
      // ------------------------------------------------------------------------
      Announce("Timestepping with ARS443");
      
      iStages = 5;
      iQorder = 3;
      iPorder = 0;
      
      pci  = new double [iStages];
      pce  = new double [iStages];
      pAi  = new double [iStages * iStages];
      pAe  = new double [iStages * iStages];
      pbi  = new double [iStages];
      pbe  = new double [iStages];
            
      // Implicit table
      pci[0] = 0.0;
      pci[1] = 0.5;
      pci[2] = 2.0/3.0;
      pci[3] = 0.5;
      pci[4] = 1.0;
      
      pAi[0]  = 0.0;  pAi[1]  = 0.0;      pAi[2]  = 0.0;   pAi[3]  = 0.0;  pAi[4]  = 0.0;
      pAi[5]  = 0.0;  pAi[6]  = 0.5;      pAi[7]  = 0.0;   pAi[8]  = 0.0;  pAi[9]  = 0.0;
      pAi[10] = 0.0;  pAi[11] = 1.0/6.0;  pAi[12] = 0.5;   pAi[13] = 0.0;  pAi[14] = 0.0;
      pAi[15] = 0.0;  pAi[16] = -0.5;     pAi[17] = 0.5;   pAi[18] = 0.5;  pAi[19] = 0.0;
      pAi[20] = 0.0;  pAi[21] = 1.5;      pAi[22] = -1.5;  pAi[23] = 0.5;  pAi[24] = 0.5;
      
      pbi[0] = 0.0;
      pbi[1] = 1.5;
      pbi[2] = -1.5;
      pbi[3] = 0.5;
      pbi[4] = 0.5;

      // Explicit table
      pce[0] = 0.0;
      pce[1] = 0.5;
      pce[2] = 2.0/3.0;
      pce[3] = 0.5;
      pce[4] = 1.0;

      pAe[0]  = 0.0;        pAe[1]  = 0.0;       pAe[2]  = 0.0;   pAe[3]  = 0.0;       pAe[4]  = 0.0;
      pAe[5]  = 0.5;        pAe[6]  = 0.0;       pAe[7]  = 0.0;   pAe[8]  = 0.0;       pAe[9]  = 0.0;
      pAe[10] = 11.0/18.0;  pAe[11] = 1.0/18.0;  pAe[12] = 0.0;   pAe[13] = 0.0;       pAe[14] = 0.0;
      pAe[15] = 5.0/6.0;    pAe[16] = -5.0/6.0;  pAe[17] = 0.5;   pAe[18] = 0.0;       pAe[19] = 0.0;
      pAe[20] = 0.25;       pAe[21] = 7.0/4.0;   pAe[22] = 0.75;  pAe[23] = -7.0/4.0;  pAe[24] = 0.0;
      
      pbe[0] = 0.25;
      pbe[1] = 7.0/4.0;
      pbe[2] = 0.75;
      pbe[3] = -7.0/4.0;
      pbe[4] = 0.0;
            
      // arkode memory, stages, order, emdedding order, ci, ce, Ai, Ae, bi, be, b2i, b2e
      ierr = ARKodeSetARKTables(arkode_mem, iStages, iQorder, iPorder, 
				pci, pce, pAi, pAe, pbi, pbe, NULL, NULL);

      delete[] pci;     
      delete[] pce;
      delete[] pAi;
      delete[] pAe;
      delete[] pbi;
      delete[] pbe;

    } else if (m_strButcherTable == "ark232") {      

      // ------------------------------------------------------------------------
      // ARK232 - 2 implicit stages, 3 explicit stages, 2nd order
      //          L-stable
      //
      // Giraldo, Kelly, and Constantinescu, Implicit-Explicit formulations of a 
      // three dimensional nonhydrostatic unificed model of the atmosphere 
      // (NUMA), 2013. (equation 3.10)
      // ------------------------------------------------------------------------
      Announce("Timestepping with ARK232");
      
      iStages = 3;
      iQorder = 2;
      iPorder = 0;

      pci  = new double [iStages];     
      pce  = new double [iStages];
      pAi  = new double [iStages * iStages];
      pAe  = new double [iStages * iStages];
      pbi  = new double [iStages];
      pbe  = new double [iStages];
      
      double gamma = 1.0 - 1.0/std::sqrt(2.0);
      double alpha = 1.0/6.0 * (3.0 + 2.0 * std::sqrt(2.0)); // 0.5 + sqrt(2)/3
      double delta = 1.0 / (2.0 * std::sqrt(2.0)); // sqrt(2)/4

      double twogamma = 2.0 * gamma; // 2 - sqrt(2)
      
      // Implicit table
      pci[0] = 0.0;
      pci[1] = twogamma;
      pci[2] = 1.0;
      
      pAi[0] = 0.0;   pAi[1] = 0.0;   pAi[2] = 0.0;
      pAi[3] = gamma; pAi[4] = gamma; pAi[5] = 0.0;
      pAi[6] = delta; pAi[7] = delta; pAi[8] = gamma;

      pbi[0] = delta;
      pbi[1] = delta;
      pbi[2] = gamma;

      // Explicit table
      pce[0] = 0.0;
      pce[1] = twogamma;
      pce[2] = 1.0;
      
      pAe[0] = 0.0;         pAe[1] = 0.0;   pAe[2] = 0.0;
      pAe[3] = twogamma;    pAe[4] = 0.0;   pAe[5] = 0.0;
      pAe[6] = 1.0 - alpha; pAe[7] = alpha; pAe[8] = 0.0;
      
      pbe[0] = delta;
      pbe[1] = delta;
      pbe[2] = gamma;
            
      // arkode memory, stages, order, emdedding order, ci, ce, Ai, Ae, bi, be, b2i, b2e
      ierr = ARKodeSetARKTables(arkode_mem, iStages, iQorder, iPorder, 
				pci, pce, pAi, pAe, pbi, pbe, NULL, NULL);     

      delete[] pci;     
      delete[] pce;
      delete[] pAi;
      delete[] pAe;
      delete[] pbi;
      delete[] pbe;

    } else if (m_strButcherTable == "ssp2_222") {

      // ------------------------------------------------------------------------
      // ssp2(2,2,2) - 2 implicit stages, 2 explicit stages, 2nd order
      // 
      // Pareschi and Russo, Implicit-explicit Runge-Kutta schemes and 
      // application to hyperbolic systems with relaxation, 2005. (Table 2)     
      // ------------------------------------------------------------------------
      Announce("Timestepping with SSP2(2,2,2)");
      
      iStages = 2;
      iQorder = 2;
      iPorder = 0;
      
      pci  = new double [iStages];
      pce  = new double [iStages];
      pAi  = new double [iStages * iStages];
      pAe  = new double [iStages * iStages];
      pbi  = new double [iStages];
      pbe  = new double [iStages];
      
      double gamma = 1.0 - 1.0 / std::sqrt(2.0);
      
      // Implicit Table
      pci[0] = gamma;
      pci[1] = 1.0 - gamma;
      
      pAi[0] = gamma;           pAi[1] = 0.0;
      pAi[2] = 1 - 2.0 * gamma; pAi[3] = gamma;

      pbi[0] = 0.5;
      pbi[1] = 0.5;

      // Explicit Table
      pce[0] = 0.0;
      pce[1] = 1.0;
      
      pAe[0] = 0.0; pAe[1] = 0.0;
      pAe[2] = 1.0; pAe[3] = 0.0;
      
      pbe[0] = 0.5;
      pbe[1] = 0.5;
         
      // arkode memory, stages, order, emdedding order, ci, ce, Ai, Ae, bi, be, b2i, b2e
      ierr = ARKodeSetARKTables(arkode_mem, iStages, iQorder, iPorder, 
				pci, pce, pAi, pAe, pbi, pbe, NULL, NULL);

      delete[] pci;      
      delete[] pce;
      delete[] pAi;
      delete[] pAe;
      delete[] pbi;
      delete[] pbe;

    } else if (m_strButcherTable == "ssp2_332a") {

      // ------------------------------------------------------------------------
      // ssp2(3,3,2)a - 3 implicit stages, 3 explicit stages, 2nd order
      //
      // Pareschi and Russo, Implicit-explicit Runge-Kutta schemes and 
      // application to hyperbolic systems with relaxation, 2005. (Table 4)
      // ------------------------------------------------------------------------
      Announce("Timestepping with SSP2(3,3,2)a");
      
      iStages = 3;
      iQorder = 2;
      iPorder = 0;
      
      pci  = new double [iStages];
      pce  = new double [iStages];
      pAi  = new double [iStages * iStages];
      pAe  = new double [iStages * iStages];
      pbi  = new double [iStages];
      pbe  = new double [iStages];
            
      // Implicit Table
      pci[0] = 0.25;
      pci[1] = 0.25;
      pci[2] = 1.0;
      
      pAi[0] = 0.25;       pAi[1] = 0.0;        pAi[2] = 0.0;
      pAi[3] = 0.0;        pAi[4] = 0.25;       pAi[5] = 0.0;
      pAi[6] = 1.0 / 3.0;  pAi[7] = 1.0 / 3.0;  pAi[8] = 1.0 / 3.0;

      pbi[0] = 1.0 / 3.0;
      pbi[1] = 1.0 / 3.0;
      pbi[2] = 1.0 / 3.0;
      
      // Explicit Table
      pce[0] = 0.0;
      pce[1] = 0.5;
      pce[2] = 1.0;

      pAe[0] = 0.0;  pAe[1] = 0.0;  pAe[2] = 0.0;
      pAe[3] = 0.5;  pAe[4] = 0.0;  pAe[5] = 0.0;
      pAe[6] = 0.5;  pAe[7] = 0.5;  pAe[8] = 0.0;
      
      pbe[0] = 1.0 / 3.0;
      pbe[1] = 1.0 / 3.0;
      pbe[2] = 1.0 / 3.0;
            
      // arkode memory, stages, order, emdedding order, ci, ce, Ai, Ae, bi, be, b2i, b2e
      ierr = ARKodeSetARKTables(arkode_mem, iStages, iQorder, iPorder, 
				pci, pce, pAi, pAe, pbi, pbe, NULL, NULL);

      delete[] pci;      
      delete[] pce;
      delete[] pAi;
      delete[] pAe;
      delete[] pbi;
      delete[] pbe;

    } else if (m_strButcherTable == "ssp3_332") {

      // ------------------------------------------------------------------------
      // ssp3(3,3,2) - 3 implicit stages, 3 explicit stages, 2nd order, L-stable
      //
      // Pareschi and Russo, Implicit-explicit Runge-Kutta schemes and 
      // application to hyperbolic systems with relaxation, 2005. (Table 5)
      // ------------------------------------------------------------------------
      Announce("Timestepping with SSP3(3,3,2)");
      
      iStages = 3;
      iQorder = 2;
      iPorder = 0;
      
      pci  = new double [iStages];
      pce  = new double [iStages];
      pAi  = new double [iStages * iStages];
      pAe  = new double [iStages * iStages];
      pbi  = new double [iStages];
      pbe  = new double [iStages];

      double gamma = 1.0 - 1.0 / std::sqrt(2.0); 
            
      // Implicit table
      pci[0] = gamma;
      pci[1] = 1.0 - gamma;
      pci[2] = 0.5;
      
      pAi[0] = gamma;              pAi[1] = 0.0;    pAi[2] = 0.0;
      pAi[3] = 1.0 - 2.0 * gamma;  pAi[4] = gamma;  pAi[5] = 0.0;
      pAi[6] = 0.5 - gamma;        pAi[7] = 0.0;    pAi[8] = gamma;

      pbi[0] = 1.0 / 6.0;
      pbi[1] = 1.0 / 6.0;
      pbi[2] = 2.0 / 3.0;

      // Explicit table
      pce[0] = 0.0;
      pce[1] = 1.0;
      pce[2] = 0.5;
      
      pAe[0] = 0.0;   pAe[1] = 0.0;   pAe[2] = 0.0;
      pAe[3] = 1.0;   pAe[4] = 0.0;   pAe[5] = 0.0;
      pAe[6] = 0.25;  pAe[7] = 0.25;  pAe[8] = 0.0;    

      pbe[0] = 1.0 / 6.0;
      pbe[1] = 1.0 / 6.0;
      pbe[2] = 2.0 / 3.0;
            
      // arkode memory, stages, order, emdedding order, ci, ce, Ai, Ae, bi, be, b2i, b2e
      ierr = ARKodeSetARKTables(arkode_mem, iStages, iQorder, iPorder, 
				pci, pce, pAi, pAe, pbi, pbe, NULL, NULL);

      delete[] pci;      
      delete[] pce;
      delete[] pAi;
      delete[] pAe;
      delete[] pbi;
      delete[] pbe;

    } else if (m_strButcherTable == "ssp3_433") {

      // ------------------------------------------------------------------------
      // ssp3(4,3,3) - 4 implicit stages, 3 explicit stages, 3rd order, L-stable
      //
      // Pareschi and Russo, Implicit-explicit Runge-Kutta schemes and 
      // application to hyperbolic systems with relaxation, 2005. (Table 6) 
      // ------------------------------------------------------------------------
      Announce("Timestepping with SSP3(4,3,3)");
      
      iStages = 4;
      iQorder = 3;
      iPorder = 0;
      
      pci  = new double [iStages];
      pce  = new double [iStages];
      pAi  = new double [iStages * iStages];
      pAe  = new double [iStages * iStages];
      pbi  = new double [iStages];
      pbe  = new double [iStages];

      double alpha = 0.24169426078821; 
      double beta  = 0.06042356519705;
      double eta   = 0.12915286960590;
      double delta = 0.5 - beta - eta - alpha;    

      // Implicit table
      pci[0] = alpha;
      pci[1] = 0.0;
      pci[2] = 1.0;
      pci[3] = 0.5;
      
      pAi[0]  = alpha;   pAi[1]  = 0.0;        pAi[2]  = 0.0;    pAi[3]  = 0.0; 
      pAi[4]  = -alpha;  pAi[5]  = alpha;      pAi[6]  = 0.0;    pAi[7]  = 0.0;
      pAi[8]  = 0.0;     pAi[9]  = 1 - alpha;  pAi[10] = alpha;  pAi[11] = 0.0; 
      pAi[12] = beta;    pAi[13] = eta;        pAi[14] = delta;  pAi[15] = alpha; 

      pbi[0] = 0.0;
      pbi[1] = 1.0 / 6.0;
      pbi[2] = 1.0 / 6.0;
      pbi[3] = 2.0 / 3.0;

      // Explicit table
      pce[0] = 0.0;
      pce[1] = 0.0;
      pce[2] = 1.0;
      pce[3] = 0.5;
      
      pAe[0]  = 0.0;  pAe[1]  = 0.0;   pAe[2]  = 0.0;   pAe[3]  = 0.0; 
      pAe[4]  = 0.0;  pAe[5]  = 0.0;   pAe[6]  = 0.0;   pAe[7]  = 0.0;
      pAe[8]  = 0.0;  pAe[9]  = 1.0;   pAe[10] = 0.0;   pAe[11] = 0.0; 
      pAe[12] = 0.0;  pAe[13] = 0.25;  pAe[14] = 0.25;  pAe[15] = 0.0; 

      pbe[0] = 0.0;     
      pbe[1] = 1.0 / 6.0;
      pbe[2] = 1.0 / 6.0;
      pbe[3] = 2.0 / 3.0;
            
      // arkode memory, stages, order, emdedding order, ci, ce, Ai, Ae, bi, be, b2i, b2e
      ierr = ARKodeSetARKTables(arkode_mem, iStages, iQorder, iPorder, 
				pci, pce, pAi, pAe, pbi, pbe, NULL, NULL);

      delete[] pci;      
      delete[] pce;
      delete[] pAi;
      delete[] pAe;
      delete[] pbi;
      delete[] pbe;

    } else if (m_strButcherTable == "ssp2_332b") {

      // ------------------------------------------------------------------------
      // ssp2(3,3,2)b - 3 implicit stages, 3 explicit stages, 2nd order
      //                modification of ssp2(3,3,2)a
      //
      // Higueras, Strong stability for additive Runge-Kutta methods, 2006.
      // ------------------------------------------------------------------------
      Announce("Timestepping with SSP2(3,3,2)b");
      
      iStages = 3;
      iQorder = 2;
      iPorder = 0;
      
      pci  = new double [iStages];
      pce  = new double [iStages];
      pAi  = new double [iStages * iStages];
      pAe  = new double [iStages * iStages];
      pbi  = new double [iStages];
      pbe  = new double [iStages];

      // Implicit table
      pci[0] = 0.2;
      pci[1] = 0.3;
      pci[2] = 1.0;
      
      pAi[0] = 0.2;      pAi[1] = 0.0;      pAi[2] = 0.0;   
      pAi[3] = 0.1;      pAi[4] = 0.2;      pAi[5] = 0.0;  
      pAi[6] = 1.0/3.0;  pAi[7] = 1.0/3.0;  pAi[8] = 1.0/3.0; 

      pbi[0] = 1.0 / 3.0;
      pbi[1] = 1.0 / 3.0;
      pbi[2] = 1.0 / 3.0;

      // Explicit table
      pce[0] = 0.0;
      pce[1] = 0.5;
      pce[2] = 1.0;
      
      pAe[0] = 0.0;  pAe[1] = 0.0;  pAe[2] = 0.0;   
      pAe[3] = 0.5;  pAe[4] = 0.0;  pAe[5] = 0.0;  
      pAe[6] = 0.5;  pAe[7] = 0.5;  pAe[8] = 0.0; 

      pbe[0] = 1.0 / 3.0;
      pbe[1] = 1.0 / 3.0;
      pbe[2] = 1.0 / 3.0;
            
      // arkode memory, stages, order, emdedding order, ci, ce, Ai, Ae, bi, be, b2i, b2e
      ierr = ARKodeSetARKTables(arkode_mem, iStages, iQorder, iPorder, 
				pci, pce, pAi, pAe, pbi, pbe, NULL, NULL);

      delete[] pci;      
      delete[] pce;
      delete[] pAi;
      delete[] pAe;
      delete[] pbi;
      delete[] pbe;

    } else if (m_strButcherTable == "ssp3_333") {

      // ------------------------------------------------------------------------
      // ssp3(3,3,3) - 3 implicit stages, 3 explicit stages, 3rd order
      //
      // Higueras, Characterizing strong stability preserving additive 
      // Runge-Kutta methos, 2009.
      // ------------------------------------------------------------------------
      Announce("Timestepping with SSP3(3,3,3)");
      
      iStages = 3;
      iQorder = 3;
      iPorder = 0;
      
      pci  = new double [iStages];
      pce  = new double [iStages];
      pAi  = new double [iStages * iStages];
      pAe  = new double [iStages * iStages];
      pbi  = new double [iStages];
      pbe  = new double [iStages];

      // Implicit table
      pci[0] = 0.0;
      pci[1] = 1.0;
      pci[2] = 0.5;
      
      pAi[0] = 0.0;        pAi[1] = 0.0;       pAi[2] = 0.0;   
      pAi[3] = 14.0/15.0;  pAi[4] = 1.0/15.0;  pAi[5] = 0.0;  
      pAi[6] = 7.0/30.0;   pAi[7] = 0.2;       pAi[8] = 1.0/15.0; 

      pbi[0] = 1.0 / 6.0;
      pbi[1] = 1.0 / 6.0;
      pbi[2] = 2.0 / 3.0;

      // Explicit table
      pce[0] = 0.0;
      pce[1] = 1.0;
      pce[2] = 0.5;
      
      pAe[0] = 0.0;   pAe[1] = 0.0;   pAe[2] = 0.0;   
      pAe[3] = 1.0;   pAe[4] = 0.0;   pAe[5] = 0.0;  
      pAe[6] = 0.25;  pAe[7] = 0.25;  pAe[8] = 0.0; 

      pbe[0] = 1.0 / 6.0;
      pbe[1] = 1.0 / 6.0;
      pbe[2] = 2.0 / 3.0;
            
      // arkode memory, stages, order, emdedding order, ci, ce, Ai, Ae, bi, be, b2i, b2e
      ierr = ARKodeSetARKTables(arkode_mem, iStages, iQorder, iPorder, 
				pci, pce, pAi, pAe, pbi, pbe, NULL, NULL);

      delete[] pci;      
      delete[] pce;
      delete[] pAi;
      delete[] pAe;
      delete[] pbi;
      delete[] pbe;

    } else if (m_strButcherTable == "ssp2_332_lspum") {

      // ------------------------------------------------------------------------
      // ssp2(3,3,2)-LSPUM - 3 implicit stages, 3 explicit stages, 2nd order
      //                     L-stable
      //
      // Higueras, Happenhofer, Koch, and Kupka, Optimized strong stability 
      // preserving IMEX Runge-Kutta methods, 2014. (equation 17)
      // ------------------------------------------------------------------------
      Announce("Timestepping with SSP2(3,3,2)-LSPUM");
      
      iStages = 3;
      iQorder = 2;
      iPorder = 0;
      
      pci  = new double [iStages];
      pce  = new double [iStages];
      pAi  = new double [iStages * iStages];
      pAe  = new double [iStages * iStages];
      pbi  = new double [iStages];
      pbe  = new double [iStages];

      // Implicit table
      pci[0] = 2.0/11.0;
      pci[1] = 289.0/462.0;
      pci[2] = 751.0/924.0;
      
      pAi[0] = 2.0/11.0;       pAi[1] = 0.0;         pAi[2] = 0.0;   
      pAi[3] = 205.0/462.0;    pAi[4] = 2.0/11.0;    pAi[5] = 0.0;  
      pAi[6] = 2033.0/4620.0;  pAi[7] = 21.0/110.0;  pAi[8] = 2.0/11.0; 

      pbi[0] = 24.0/55.0;
      pbi[1] = 0.2;
      pbi[2] = 4.0/11.0;

      // Explicit table
      pce[0] = 0.0;
      pce[1] = 5.0/6.0;
      pce[2] = 11.0/12.0;
      
      pAe[0] = 0.0;        pAe[1] = 0.0;        pAe[2] = 0.0;   
      pAe[3] = 5.0/6.0;    pAe[4] = 0.0;        pAe[5] = 0.0;  
      pAe[6] = 11.0/24.0;  pAe[7] = 11.0/24.0;  pAe[8] = 0.0; 

      pbe[0] = 24.0/55.0;
      pbe[1] = 0.2;
      pbe[2] = 4.0/11.0;
            
      // arkode memory, stages, order, emdedding order, ci, ce, Ai, Ae, bi, be, b2i, b2e
      ierr = ARKodeSetARKTables(arkode_mem, iStages, iQorder, iPorder, 
				pci, pce, pAi, pAe, pbi, pbe, NULL, NULL);

      delete[] pci;      
      delete[] pce;
      delete[] pAi;
      delete[] pAe;
      delete[] pbi;
      delete[] pbe;

    } else if (m_strButcherTable == "ssp2_332_lpum") {

      // ------------------------------------------------------------------------
      // ssp2(3,3,2)-LPUM - 3 implicit stages, 3 explicit stages, 2nd order
      //                    L-stable
      //
      // Higueras, Happenhofer, Koch, and Kupka, Optimized strong stability 
      // preserving IMEX Runge-Kutta methods, 2014. (equation 20)
      // ------------------------------------------------------------------------
      Announce("Timestepping with SSP2(3,3,2)-LPUM");
      
      iStages = 3;
      iQorder = 2;
      iPorder = 0;
      
      pci  = new double [iStages];
      pce  = new double [iStages];
      pAi  = new double [iStages * iStages];
      pAe  = new double [iStages * iStages];
      pbi  = new double [iStages];
      pbe  = new double [iStages];

      // Implicit table
      pci[0] = 2.0/11.0;
      pci[1] = 69.0/154.0;
      pci[2] = 67.0/77.0;
      
      pAi[0] = 2.0/11.0;     pAi[1] = 0.0;         pAi[2] = 0.0;   
      pAi[3] = 41.0/154.0;   pAi[4] = 2.0/11.0;    pAi[5] = 0.0;  
      pAi[6] = 289.0/847.0;  pAi[7] = 42.0/121.0;  pAi[8] = 2.0/11.0; 

      pbi[0] = 1.0/3.0;
      pbi[1] = 1.0/3.0;
      pbi[2] = 1.0/3.0;

      // Explicit table
      pce[0] = 0.0;
      pce[1] = 0.5;
      pce[2] = 1.0;
      
      pAe[0] = 0.0;  pAe[1] = 0.0;  pAe[2] = 0.0;   
      pAe[3] = 0.5;  pAe[4] = 0.0;  pAe[5] = 0.0;  
      pAe[6] = 0.5;  pAe[7] = 0.5;  pAe[8] = 0.0; 

      pbe[0] = 1.0/3.0;
      pbe[1] = 1.0/3.0;
      pbe[2] = 1.0/3.0;
            
      // arkode memory, stages, order, emdedding order, ci, ce, Ai, Ae, bi, be, b2i, b2e
      ierr = ARKodeSetARKTables(arkode_mem, iStages, iQorder, iPorder, 
				pci, pce, pAi, pAe, pbi, pbe, NULL, NULL);

      delete[] pci;      
      delete[] pce;
      delete[] pAi;
      delete[] pAe;
      delete[] pbi;
      delete[] pbe;

    } else if (m_strButcherTable == "ssp2_332_lpm1") {

      // ------------------------------------------------------------------------
      // ssp2(3,3,2)-LPM1 - 3 implicit stages, 3 explicit stages, 2nd order
      //                    L-stable
      //
      // Higueras, Happenhofer, Koch, and Kupka, Optimized strong stability 
      // preserving IMEX Runge-Kutta methods, 2014. (equation 22)
      // ------------------------------------------------------------------------
      Announce("Timestepping with SSP2(3,3,2)-LPM1");
      
      iStages = 3;
      iQorder = 2;
      iPorder = 0;
      
      pci  = new double [iStages];
      pce  = new double [iStages];
      pAi  = new double [iStages * iStages];
      pAe  = new double [iStages * iStages];
      pbi  = new double [iStages];
      pbe  = new double [iStages];

      // Implicit table
      pci[0] = 2.0/11.0;
      pci[1] = 4523.0/9317.0;
      pci[2] = 15517.0/18634.0;
      
      pAi[0] = 2.0/11.0;           pAi[1] = 0.0;       pAi[2] = 0.0;   
      pAi[3] = 2829.0/9317.0;      pAi[4] = 2.0/11.0;  pAi[5] = 0.0;  
      pAi[6] = 148529.0/428582.0;  pAi[7] = 7.0/23.0;  pAi[8] = 2.0/11.0; 

      pbi[0] = 1.0/3.0;
      pbi[1] = 1.0/3.0;
      pbi[2] = 1.0/3.0;

      // Explicit table
      pce[0] = 0.0;
      pce[1] = 0.5;
      pce[2] = 1.0;
      
      pAe[0] = 0.0;  pAe[1] = 0.0;  pAe[2] = 0.0;   
      pAe[3] = 0.5;  pAe[4] = 0.0;  pAe[5] = 0.0;  
      pAe[6] = 0.5;  pAe[7] = 0.5;  pAe[8] = 0.0; 

      pbe[0] = 1.0/3.0;
      pbe[1] = 1.0/3.0;
      pbe[2] = 1.0/3.0;
            
      // arkode memory, stages, order, emdedding order, ci, ce, Ai, Ae, bi, be, b2i, b2e
      ierr = ARKodeSetARKTables(arkode_mem, iStages, iQorder, iPorder, 
				pci, pce, pAi, pAe, pbi, pbe, NULL, NULL);

      delete[] pci;      
      delete[] pce;
      delete[] pAi;
      delete[] pAe;
      delete[] pbi;
      delete[] pbe;

    } else if (m_strButcherTable == "ssp2_332_lpm2") {

      // ------------------------------------------------------------------------
      // ssp2(3,3,2)-LPM2 - 3 implicit stages, 3 explicit stages, 2nd order
      //                    L-stable
      //
      // Higueras, Happenhofer, Koch, and Kupka, Optimized strong stability 
      // preserving IMEX Runge-Kutta methods, 2014. (equation 23)
      // ------------------------------------------------------------------------
      Announce("Timestepping with SSP2(3,3,2)-LPM2");
      
      iStages = 3;
      iQorder = 2;
      iPorder = 0;
      
      pci  = new double [iStages];
      pce  = new double [iStages];
      pAi  = new double [iStages * iStages];
      pAe  = new double [iStages * iStages];
      pbi  = new double [iStages];
      pbe  = new double [iStages];

      // Implicit table
      pci[0] = 2.0/11.0;
      pci[1] = 5003.0/13310.0;
      pci[2] = 6271.0/6655.0;
      
      pAi[0] = 2.0/11.0;          pAi[1] = 0.0;        pAi[2] = 0.0;   
      pAi[3] = 2583.0/13310.0;    pAi[4] = 2.0/11.0;   pAi[5] = 0.0;  
      pAi[6] = 39731.0/139755.0;  pAi[7] = 10.0/21.0;  pAi[8] = 2.0/11.0; 

      pbi[0] = 1.0/3.0;
      pbi[1] = 1.0/3.0;
      pbi[2] = 1.0/3.0;

      // Explicit table
      pce[0] = 0.0;
      pce[1] = 0.5;
      pce[2] = 1.0;
      
      pAe[0] = 0.0;  pAe[1] = 0.0;  pAe[2] = 0.0;   
      pAe[3] = 0.5;  pAe[4] = 0.0;  pAe[5] = 0.0;  
      pAe[6] = 0.5;  pAe[7] = 0.5;  pAe[8] = 0.0; 

      pbe[0] = 1.0/3.0;
      pbe[1] = 1.0/3.0;
      pbe[2] = 1.0/3.0;
            
      // arkode memory, stages, order, emdedding order, ci, ce, Ai, Ae, bi, be, b2i, b2e
      ierr = ARKodeSetARKTables(arkode_mem, iStages, iQorder, iPorder, 
				pci, pce, pAi, pAe, pbi, pbe, NULL, NULL);

      delete[] pci;      
      delete[] pce;
      delete[] pAi;
      delete[] pAe;
      delete[] pbi;
      delete[] pbe;

    } else {
      _EXCEPTIONT("ERROR: Invalid IMEX Butcher table name");
    }   
  }

  if (ierr < 0) _EXCEPTION1("ERROR: SetButcherTable, ierr = %i",ierr);
 
}

///////////////////////////////////////////////////////////////////////////////

void TimestepSchemeARKode::AssignComponentWiseTolerances() {
  
      // index of various N_Vector arguments in registry
      int iT = NV_INDEX_TEMPEST(m_T);

      // Get a copy of the grid
      Grid * pGrid = NV_GRID_TEMPEST(m_T);

      // Tolerances are assigned in Tempest (objects Grid and GridPatch)
      pGrid->AssignComponentWiseTolerances(iT);

}

///////////////////////////////////////////////////////////////////////////////

#endif
