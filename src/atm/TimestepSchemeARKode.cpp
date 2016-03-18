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

#include "TimestepSchemeARKode.h"
#include "Model.h"
#include "Grid.h"
#include "HorizontalDynamics.h"
#include "VerticalDynamics.h"
#include "HorizontalDynamicsFEM.h"
#include "VerticalDynamicsFEM.h"
#include "Announce.h"

void * TimestepSchemeARKode::ARKodeMem = NULL;

///////////////////////////////////////////////////////////////////////////////

TimestepSchemeARKode::TimestepSchemeARKode(
        Model & model,
	ARKodeCommandLineVariables & ARKodeVars
) :
        TimestepScheme(model),
	m_iNVectors(ARKodeVars.nvectors),
	m_iARKodeButcherTable(ARKodeVars.ARKodeButcherTable),
	m_iSetButcherTable(ARKodeVars.SetButcherTable),
	m_dRelTol(ARKodeVars.rtol),
	m_dAbsTol(ARKodeVars.atol),
	m_fFullyExplicit(ARKodeVars.FullyExplicit),
	m_fFullyImplicit(false),
	m_fDynamicStepSize(ARKodeVars.DynamicStepSize),
	m_dDynamicDeltaT(0.0),
	m_fAAFP(ARKodeVars.AAFP),
	m_iAAFPAccelVec(ARKodeVars.AAFPAccelVec),
	m_iNonlinIters(ARKodeVars.NonlinIters),
	m_iLinIters(ARKodeVars.LinIters),
	m_fWriteDiagnostics(ARKodeVars.WriteDiagnostics),
	m_fUsePreconditioning(ARKodeVars.UsePreconditioning)
{
        // Allocate ARKode memory
        ARKodeMem = ARKodeCreate();

	if (ARKodeMem == NULL) _EXCEPTIONT("ERROR: ARKodeCreate returned NULL");

	// Check input paramters
	if (m_iARKodeButcherTable >= 0 && m_iSetButcherTable >= 0)
	  _EXCEPTIONT("ERROR: ARKodeButcherTable and SetButcherTable are both >= 0.");
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
    ierr = ARKodeInit(ARKodeMem, ARKodeFullRHS, NULL, dCurrentT, m_Y);
  } else if (m_fFullyImplicit) {
    ierr = ARKodeInit(ARKodeMem, NULL, ARKodeFullRHS, dCurrentT, m_Y);
  } else {
    ierr = ARKodeInit(ARKodeMem, ARKodeExplicitRHS, ARKodeImplicitRHS, dCurrentT, m_Y);
  }

  if (ierr < 0) _EXCEPTION1("ERROR: ARKodeInit, ierr = %i",ierr);

  // Set stop time 
  Time timeEndT = m_model.GetEndTime();
  double dEndT  = timeEndT.GetSeconds();

  ierr = ARKodeSetStopTime(ARKodeMem, dEndT);
  if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetStopTime, ierr = %i",ierr);  

  // Set function to post process arkode steps
  ierr = ARKodeSetPostprocessStepFn(ARKodeMem, ARKodePostProcessStep);
  if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetPostprocessStepFn, ierr = %i",ierr);

  // Set adaptive or fixed time stepping  
  if (m_fDynamicStepSize) {

    // set dynamic timestepping flag to true
    m_model.SetDynamicTimestepping(m_fDynamicStepSize);

  } else {

    // Set fixed step size in seconds
    Time timeDeltaT = m_model.GetDeltaT();
    double dDeltaT  = timeDeltaT.GetSeconds();
   
    ierr = ARKodeSetFixedStep(ARKodeMem, dDeltaT);   
    if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetFixedStep, ierr = %i",ierr);
  }

  // Select ARKode Butcher table
  if (m_iARKodeButcherTable >= 0) {

    if (m_fFullyExplicit) {

      ierr = ARKodeSetERKTableNum(ARKodeMem, m_iARKodeButcherTable);
      if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetERKTableNum, ierr = %i",ierr);

    } else if (m_fFullyImplicit) {

      ierr = ARKodeSetIRKTableNum(ARKodeMem, m_iARKodeButcherTable);
      if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetIRKTableNum, ierr = %i",ierr);

    } else {

      if (m_iARKodeButcherTable == 2 || m_iARKodeButcherTable == 16) {

	ierr = ARKodeSetARKTableNum(ARKodeMem, 16, 2);

      } else if (m_iARKodeButcherTable == 4 || m_iARKodeButcherTable == 21) {

	ierr = ARKodeSetARKTableNum(ARKodeMem, 21, 4);

      } else if (m_iARKodeButcherTable == 9 || m_iARKodeButcherTable == 23) {

	ierr = ARKodeSetARKTableNum(ARKodeMem, 23, 9);

      } else {
	ierr = ARK_ILL_INPUT;
      }

      if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetARKTableNum, ierr = %i",ierr);
    }
  } 

  // Set a user supplied Butcher table
  if (m_iSetButcherTable >= 0) {  
    SetButcherTable();
  }  
  
  // Specify tolerances
  ierr = ARKodeSStolerances(ARKodeMem, m_dRelTol, m_dAbsTol);
  if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSStolerances, ierr = %i",ierr);

  // Nonlinear Solver Settings
  if (!m_fFullyExplicit) {
    
    if (m_fAAFP) { // Anderson accelerated fixed point solver
      
      ierr = ARKodeSetFixedPoint(ARKodeMem, m_iAAFPAccelVec);
      if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetFixedPoint, ierr = %i",ierr);

    } else { // Newton iteration

      int precflag = PREC_NONE;
      if (m_fUsePreconditioning)  precflag = PREC_RIGHT;

      // Linear Solver Settings
      ierr = ARKSpgmr(ARKodeMem, precflag, m_iLinIters);      
      if (ierr < 0) _EXCEPTION1("ERROR: ARKSpgmr, ierr = %i",ierr);

      // Use preconditioning if requested
      if (m_fUsePreconditioning) {
	ierr = ARKSpilsSetPreconditioner(ARKodeMem, NULL, ARKodePreconditionerSolve);
	if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetPreconditioner, ierr = %i",ierr);
      }
  
    }

    // if negative nonlinear iterations are specified, switch to linear-implicit mode
    if (m_iNonlinIters < 0) {
      
      ierr = ARKodeSetLinear(ARKodeMem, 1);
      if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetLinear, ierr = %i",ierr);
      
      // otherwise, set the Max nonlinear solver iterations
    } else { 
      
      ierr = ARKodeSetMaxNonlinIters(ARKodeMem, m_iNonlinIters);
      if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetMaxNonlinIters, ierr = %i",ierr);
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
      
      ierr = ARKodeSetDiagnostics(ARKodeMem, pFile);
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
	      << GetLengthTempestNVectorRegistry()
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
    ierr = ARKodeSetFixedStep(ARKodeMem, dDeltaT);
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

  // timeDeltaT is possibly too large on last step because value that 
  // GetDeltaT returns is not adjusted if necessary on the last step,
  // stop time in arkode will adjust time step to match end time
  timeNextT += timeDeltaT; 

  double dNextT = timeNextT.GetSeconds();

  // set up call to ARKode to evolve to dNextT, using either a single step or adaptivity
  int stepmode = ARK_ONE_STEP;
  if (m_fDynamicStepSize) {
    stepmode = ARK_NORMAL;
    ierr = ARKodeSetStopTime(ARKodeMem, dNextT);
    if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetStopTime, ierr = %i",ierr);
  } else {
    stepmode = ARK_ONE_STEP;
  }

  // ARKode timestep
  ierr = ARKode(ARKodeMem, dNextT, m_Y, &dCurrentT, stepmode);
  if (ierr < 0) _EXCEPTION1("ERROR: ARKode, ierr = %i",ierr);

  // // With dynamic stepping, get the last step size to update model time
  // if (m_fDynamicStepSize) {
  //   ierr = ARKodeGetLastStep(ARKodeMem, &m_dDynamicDeltaT);
  //   if (ierr < 0) _EXCEPTION1("ERROR: ARKodeGetLastStep, ierr = %i",ierr);
  // }

#ifdef STATISTICS_OUTPUT
  int iRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &iRank);
    
  if (iRank == 0) {
    long int nsteps, expsteps, accsteps, step_attempts, nfe_evals, nfi_evals, nlinsetups, 
      netfails, nniters, nncfails, npsolves, nliters, nlcfails, nfevalsLS;
    realtype hinused, hlast, hcur, tcur;
    ierr = ARKodeGetIntegratorStats(ARKodeMem, &nsteps, &expsteps, &accsteps, &step_attempts, 
				    &nfe_evals, &nfi_evals, &nlinsetups, &netfails, &hinused, 
				    &hlast, &hcur, &tcur);
    if (ierr < 0) _EXCEPTION1("ERROR: ARKodeGetIntegratorStats, ierr = %i",ierr);

    ierr = ARKodeGetNonlinSolvStats(ARKodeMem, &nniters, &nncfails);
    if (ierr < 0) _EXCEPTION1("ERROR: ARKodeGetNonlinSolvStats, ierr = %i",ierr);

    ierr = ARKSpilsGetNumPrecSolves(ARKodeMem, &npsolves);
    if (ierr < 0) _EXCEPTION1("ERROR: ARKodeGetNumPrecSolves, ierr = %i",ierr);

    ierr = ARKSpilsGetNumLinIters(ARKodeMem, &nliters);
    if (ierr < 0) _EXCEPTION1("ERROR: ARKSpilsGetNumLinIters, ierr = %i",ierr);

    ierr = ARKSpilsGetNumConvFails(ARKodeMem, &nlcfails);
    if (ierr < 0) _EXCEPTION1("ERROR: ARKSpilsGetNumConvFails, ierr = %i",ierr);

    ierr = ARKSpilsGetNumRhsEvals(ARKodeMem, &nfevalsLS);
    if (ierr < 0) _EXCEPTION1("ERROR: ARKSpilsGetNumRhsEvals, ierr = %i",ierr);

    std::cout << std::endl
              << "TimestepSchemeARKode::Step cumulative stats:\n"
              << "  steps: " << nsteps << " (" << step_attempts << " attempted)\n"
              << "  step size: " << hlast << " previous, " << hcur << " next\n"
              << "  fevals: " << nfe_evals << " exp, " << nfi_evals << " imp, " << nfevalsLS << " lin solve\n" 
              << "  nonlinear: " << nniters << " iters, " << nncfails << " failures\n"
              << "  linear: " << nliters << " iters, " << nlcfails << " fails, " << npsolves << " prec solves\n"
              << std::endl;
  }
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
		<< "N_Vector Registry Length After First Step: "
		<< GetLengthTempestNVectorRegistry()
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
  
  // Exchange
  pGrid->PostProcessSubstage(iY, DataType_State);
  pGrid->PostProcessSubstage(iY, DataType_Tracers);

  // Get last time step size
  Time timeT     = pModel->GetCurrentTime();  // model still has old time
  double dOldT   = timeT.GetSeconds();
  double dDeltaT = time - dOldT;
   
  // Apply hyperdiffusion (initial, update, temp)
  pGrid->CopyData(iY, 1, DataType_State);
  pGrid->CopyData(iY, 1, DataType_Tracers);
  
  pHorizontalDynamicsFEM->StepAfterSubCycle(iY, 1, 2, timeT, dDeltaT);
  
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

  // Exchange
  pGrid->PostProcessSubstage(iY, DataType_State);
  pGrid->PostProcessSubstage(iY, DataType_Tracers);

  // zero out iYdot
  pGrid->ZeroData(iYdot, DataType_State);
  pGrid->ZeroData(iYdot, DataType_Tracers);

  // Compute explicit RHS
  pHorizontalDynamicsFEM->StepExplicit(iY, iYdot, timeT, 1.0);

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

  // Exchange
  pGrid->PostProcessSubstage(iY, DataType_State);
  pGrid->PostProcessSubstage(iY, DataType_Tracers);

  // zero out iYdot
  pGrid->ZeroData(iYdot, DataType_State);
  pGrid->ZeroData(iYdot, DataType_Tracers);

  // Compute implicit RHS
  pVerticalDynamicsFEM->StepImplicitTermsExplicitly(iY, iYdot, timeT, 1.0);

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
  
  // Exchange
  pGrid->PostProcessSubstage(iY, DataType_State);
  pGrid->PostProcessSubstage(iY, DataType_Tracers);
    
  // zero out iYdot
  pGrid->ZeroData(iYdot, DataType_State);
  pGrid->ZeroData(iYdot, DataType_Tracers);

  // Compute full RHS
  pHorizontalDynamicsFEM->StepExplicit(iY, iYdot, timeT, 1.0);
  pVerticalDynamicsFEM->StepExplicit(iY, iYdot, timeT, 1.0);

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

  // model time (not used in RHS)
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

  // Exchange state and RHS vector data
  pGrid->PostProcessSubstage(iY, DataType_State);
  pGrid->PostProcessSubstage(iY, DataType_Tracers);
  pGrid->PostProcessSubstage(iR, DataType_State);
  pGrid->PostProcessSubstage(iR, DataType_Tracers);
    
  // Call column-wise linear solver (iR holds RHS on input, solution on output)
  pVerticalDynamicsFEM->SolveImplicit(iY, iR, timeT, gamma);

  // Copy solution into output N_Vector
  N_VScale_Tempest(1.0, R, Z);

#ifdef DEBUG_OUTPUT
  AnnounceEndBlock("Done");
#endif

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

void TimestepSchemeARKode::SetButcherTable()
{
  // error flag  
  int ierr = 0;

  // number of RK stages
  int iStages;

  // global order of accuracy for the method
  int iQorder;

  // global order of accuracy for the embedding
  int iPorder;

  // stage times
  double * pc = NULL;

  // A matrix for the implicit method
  double * pAi = NULL;

  // A matrix for the explicit method
  double * pAe = NULL;

  // b coefficient array
  double * pb = NULL;

  // b embedding array
  double * pbembed = NULL;

  if (m_fFullyExplicit) {

    if (m_iSetButcherTable == 0) {

      // Forward Euler 2 stages 1st order
      Announce("Timestepping with Forward Euler");

      iStages = 2;
      iQorder = 1;
      iPorder = 0;
      
      pAe = new double [iStages * iStages];
      pc  = new double [iStages];
      pb  = new double [iStages];
      pbembed = new double [iStages];
      
      pAe[0]  = 0.0;  pAe[1]  = 0.0; 
      pAe[2]  = 1.0;  pAe[3]  = 0.0; 
      
      pc[0] = 0.0;
      pc[1] = 1.0;
      
      pb[0] = 1.0;
      pb[1] = 0.0;
      
      pbembed[0] = 0.0;
      pbembed[1] = 0.0;
      
      ierr = ARKodeSetERKTable(ARKodeMem, iStages, iQorder, iPorder, pc, pAe, pb, pbembed);
      
      if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetERKTable, ierr = %i",ierr);
      
      delete[] pc;
      delete[] pAe;
      delete[] pb;
      delete[] pbembed; 

      
    } else if (m_iSetButcherTable == 1) {
      
      // Kinnmark Gray Ullrich ERK 6 stages 3rd order
      Announce("Timestepping with KGU(6,3)");
      
      iStages = 6;
      iQorder = 3;
      iPorder = 0;
      
      pAe = new double [iStages * iStages];
      pc  = new double [iStages];
      pb  = new double [iStages];
      pbembed = new double [iStages];
      
      pAe[0]  = 0.0;  pAe[1]  = 0.0; pAe[2]  = 0.0;       pAe[3]  = 0.0;       pAe[4]  = 0.0;  pAe[5]  = 0.0;
      pAe[6]  = 0.2;  pAe[7]  = 0.0; pAe[8]  = 0.0;       pAe[9]  = 0.0;       pAe[10] = 0.0;  pAe[11] = 0.0;
      pAe[12] = 0.0;  pAe[13] = 0.2; pAe[14] = 0.0;       pAe[15] = 0.0;       pAe[16] = 0.0;  pAe[17] = 0.0;
      pAe[18] = 0.0;  pAe[19] = 0.0; pAe[20] = 1.0 / 3.0; pAe[21] = 0.0;       pAe[22] = 0.0;  pAe[23] = 0.0;
      pAe[24] = 0.0;  pAe[25] = 0.0; pAe[26] = 0.0;       pAe[27] = 2.0 / 3.0; pAe[28] = 0.0;  pAe[29] = 0.0;
      pAe[30] = 0.25; pAe[31] = 0.0; pAe[32] = 0.0;       pAe[33] = 0.0;       pAe[34] = 0.75; pAe[35] = 0.0;
      
      pc[0] = 0.0;
      pc[1] = 0.2;
      pc[2] = 0.2;
      pc[3] = 1.0 / 3.0;
      pc[4] = 2.0 / 3.0;
      pc[5] = 1.0;
      
      pb[0] = 0.25;
      pb[1] = 0.0;
      pb[2] = 0.0;
      pb[3] = 0.0;
      pb[4] = 0.75;
      pb[5] = 0.0;
      
      pbembed[0] = 0.0;
      pbembed[1] = 0.0;
      pbembed[2] = 0.0;
      pbembed[3] = 0.0;
      pbembed[4] = 0.0;
      pbembed[5] = 0.0;
       
      ierr = ARKodeSetERKTable(ARKodeMem, iStages, iQorder, iPorder, pc, pAe, pb, pbembed);
      
      if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetERKTable, ierr = %i",ierr);
      
      delete[] pc;
      delete[] pAe;
      delete[] pb;
      delete[] pbembed; 

    } else if (m_iSetButcherTable == 2) {

      // Strong Stability Preserving ERK 5 stages 4th order
      Announce("Timestepping with SSPRK(5,4)");
      
      iStages = 5;
      iQorder = 4;
      iPorder = 0;
      
      pAe = new double [iStages * iStages];
      pc  = new double [iStages];
      pb  = new double [iStages];
      pbembed = new double [iStages];
      
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
      
      pc[0] = 0.0;
      pc[1] = pAe[5];
      pc[2] = pAe[10] + pAe[11];
      pc[3] = pAe[15] + pAe[16] + pAe[17];
      pc[4] = pAe[20] + pAe[21] + pAe[22] + pAe[23];
      
      pb[0] = beta1 * pAe[10] + beta2 * pAe[15] + beta3 * pAe[20];
      pb[1] = beta1 * pAe[11] + beta2 * pAe[16] + beta3 * pAe[21];
      pb[2] = beta2 * pAe[17] + beta3 * pAe[22];
      pb[3] = beta3 * pAe[23] + 0.063692468666290;
      pb[4] = 0.226007483236906;
      
      pbembed[0] = 0.0;
      pbembed[1] = 0.0;
      pbembed[2] = 0.0;
      pbembed[3] = 0.0;
      pbembed[4] = 0.0;
      
      ierr = ARKodeSetERKTable(ARKodeMem, iStages, iQorder, iPorder, pc, pAe, pb, pbembed);
      
      if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetERKTable, ierr = %i",ierr);
      
      delete[] pc;
      delete[] pAe;
      delete[] pb;
      delete[] pbembed;

    } else {     
      _EXCEPTIONT("ERROR: Invalid explicit Butcher table vaule");
    }
    
  } else if (m_fFullyImplicit) {
    _EXCEPTIONT("ERROR: SetButcherTable() not implemented for fully implicit");
    // ierr = ARKodeSetIRKTable(ARKodeMem, iStages, iQorder, iPorder, pc, pAi, pb, pbembed)
    // if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetIRKTableNum, ierr = %i",ierr);
  } else {

    if (m_iSetButcherTable == 0) {      

      // ARS232
      Announce("Timestepping with ARS232");
      
      iStages = 3;
      iQorder = 2;
      iPorder = 0;
      
      pc  = new double [iStages];
      pAi = new double [iStages * iStages];
      pAe = new double [iStages * iStages];
      pb  = new double [iStages];
      pbembed = new double [iStages];
      
      double gamma = 1.0 - 1.0/std::sqrt(2.0);
      double delta = -2.0 * std::sqrt(2.0) / 3.0;
      
      pc[0] = 0.0;
      pc[1] = gamma;
      pc[2] = 1.0;
      
      pAi[0] = 0.0; pAi[1] = 0.0;         pAi[2] = 0.0;
      pAi[3] = 0.0; pAi[4] = gamma;       pAi[5] = 0.0;
      pAi[6] = 0.0; pAi[7] = 1.0 - gamma; pAi[8] = gamma;
      
      pAe[0] = 0.0;   pAe[1] = 0.0;         pAe[2] = 0.0;
      pAe[3] = gamma; pAe[4] = 0.0;         pAe[5] = 0.0;
      pAe[6] = delta; pAe[7] = 1.0 - delta; pAe[8] = 0.0;
      
      pb[0] = 0.0;
      pb[1] = 1.0 - gamma;
      pb[2] = gamma;
      
      pbembed[0] = 0.0;
      pbembed[1] = 0.0;
      pbembed[2] = 0.0;
      
      ierr = ARKodeSetARKTables(ARKodeMem, iStages, iQorder, iPorder, pc, pAi, pAe, pb, pbembed);  
      
      if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetARKTableNum, ierr = %i",ierr);
      
      delete[] pc;
      delete[] pAi;
      delete[] pAe;
      delete[] pb;
      delete[] pbembed;

    } else {      
      _EXCEPTIONT("ERROR: Invalid IMEX Butcher table vaule");
    }
    
  } 
}

///////////////////////////////////////////////////////////////////////////////

#endif
