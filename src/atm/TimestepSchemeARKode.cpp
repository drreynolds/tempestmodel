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

//#define DEBUG_PRINT_ON
//#define STEP_LIKE_ARKODE
//#define TEST_STEP

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
	m_fFixedStepSize(true),
	m_fAAFP(ARKodeVars.AAFP),
	m_iAAFPAccelVec(ARKodeVars.AAFPAccelVec),
	m_iNonlinIters(ARKodeVars.NonlinIters),
	m_iLinIters(ARKodeVars.LinIters),
	m_dDynamicDeltaT(0.0),
	m_fWriteDiagnostics(ARKodeVars.WriteDiagnostics)
{
        // Allocate ARKode memory
        ARKodeMem = ARKodeCreate();

	if (ARKodeMem == NULL) _EXCEPTIONT("ERROR: ARKodeCreate returned NULL");

	// Check input paramters
	if (m_iARKodeButcherTable >= 0 && m_iSetButcherTable >= 0)
	  _EXCEPTIONT("ERROR: ARKodeButcherTable and SetButcherTable are both >= 0.");

	// Check if using adaptive time stepping
	Time timeDeltaT = m_model.GetDeltaT();
	
	if (timeDeltaT.IsZero()) {
	  m_fFixedStepSize = false;
	  Announce("Warning: Adaptive time stepping not yet tested.");
	}
}
					  
///////////////////////////////////////////////////////////////////////////////

void TimestepSchemeARKode::Initialize() {
  
  AnnounceStartBlock("Initializing ARKode");

#ifdef STEP_LIKE_ARKODE
  Announce("DEBUGGING: STEP_LIKE_ARKODE");
#endif

  // error flag
  int ierr = 0;

  // Get a copy of the grid
  Grid * pGrid = m_model.GetGrid();

  // Create a Tempest NVector for the model state
  m_Y = N_VNew_Tempest(*pGrid, m_model);

  if (m_Y == NULL) _EXCEPTIONT("ERROR: N_VNew_Tempest returned NULL");

  int iY = NV_INDEX_TEMPEST(m_Y);

  if (iY != 0) _EXCEPTION1("ERROR: iY != 0, iY = %i",iY);

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

  // Use a fixed step size
  if (m_fFixedStepSize) {

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

      if (m_iARKodeButcherTable == 2 || m_iARKodeButcherTable == 15) {

	ierr = ARKodeSetARKTableNum(ARKodeMem, 15, 2);

      } else if (m_iARKodeButcherTable == 4 || m_iARKodeButcherTable == 20) {

	ierr = ARKodeSetARKTableNum(ARKodeMem, 20, 4);

      } else if (m_iARKodeButcherTable == 9 || m_iARKodeButcherTable == 22) {

	ierr = ARKodeSetARKTableNum(ARKodeMem, 22, 9);

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
      
      // Linear Solver Settings 
      ierr = ARKSpgmr(ARKodeMem, PREC_NONE, m_iLinIters);
      
      if (ierr < 0) _EXCEPTION1("ERROR: ARKSpgmr, ierr = %i",ierr);
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

#ifdef DEBUG_PRINT_ON
  //  int iRegistryLength = GetLengthTempestNVectorRegistry();
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

#ifdef DEBUG_PRINT_ON
  AnnounceStartBlock("ARKode Step");
#endif

#ifdef TEST_STEP
  // Call step step to compare Strang and ARKode
  ARKode_Test_Step(fFirstStep, dDeltaT);
  return;
#endif

#ifdef STEP_LIKE_ARKODE
  // Debugging step function that takes and ARKode like step within Tempest
  Step_Like_ARKode(fFirstStep, fLastStep, time, dDeltaT);
  return;
#endif

  // error flag
  int ierr = 0;

  // Adjust last time step size if using fixed step sizes
  if (m_fFixedStepSize && fLastStep) {
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
  Time timeCurrentT = m_model.GetCurrentTime();
  double dCurrentT  = timeCurrentT.GetSeconds();

  // Next time in seconds
  double dNextT;
  if (m_fFixedStepSize) {
  
    Time timeNextT  = m_model.GetCurrentTime();
    Time timeDeltaT = m_model.GetDeltaT();
    timeNextT += timeDeltaT; // Possibly too large on last step because value 
                             // that GetDeltaT returns is not changed for last step

    dNextT = timeNextT.GetSeconds();

  } else {
    // should change to use next output time as dNextT

    Time timeEnd = m_model.GetEndTime();
    dNextT = timeEnd.GetSeconds();
  }
 
  // ARKode timestep
  ierr = ARKode(ARKodeMem, dNextT, m_Y, &dCurrentT, ARK_ONE_STEP);

  if (ierr < 0) _EXCEPTION1("ERROR: ARKode, ierr = %i",ierr);

  // Get dynamic step size
  if (!m_fFixedStepSize) {
    ierr = ARKodeGetLastStep(ARKodeMem, &m_dDynamicDeltaT);
      
    if (ierr < 0) _EXCEPTION1("ERROR: ARKodeGetLastStep, ierr = %i",ierr);
  }

  // Get idex of Tempest NVector in the grid
  int iY = NV_INDEX_TEMPEST(m_Y);

  // Filter negative tracers
  pHorizontalDynamicsFEM->FilterNegativeTracers(iY);
  pVerticalDynamicsFEM->FilterNegativeTracers(iY);

  // Exchange
  pGrid->PostProcessSubstage(iY, DataType_State);
  pGrid->PostProcessSubstage(iY, DataType_Tracers);

  // Adjust hyperdiffusion step size when using adaptive step sizes
  double dLastDeltaT;
  if (m_fFixedStepSize) {
    dLastDeltaT = dDeltaT;
  } else {
    dLastDeltaT = m_dDynamicDeltaT;
  }

  // Apply hyperdiffusion (initial, update, temp)
  pGrid->CopyData(iY, 1, DataType_State);
  pGrid->CopyData(iY, 1, DataType_Tracers);

  pHorizontalDynamicsFEM->StepAfterSubCycle(iY, 1, 2, time, dLastDeltaT);

  pGrid->CopyData(1, iY, DataType_State);
  pGrid->CopyData(1, iY, DataType_Tracers);

#ifdef DEBUG_PRINT_ON
  AnnounceEndBlock("Done");
  // int iRegistryLength = GetLengthTempestNVectorRegistry();
#endif
}

///////////////////////////////////////////////////////////////////////////////

static int ARKodeExplicitRHS(
	realtype time, 
	N_Vector Y, 
	N_Vector Ydot, 
	void * user_data
) {

#ifdef DEBUG_PRINT_ON
  AnnounceStartBlock("Explicit RHS");
#endif

  // model time
  Time timeT = Time(0,0,0,time,0); // only works for integers

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

#ifdef DEBUG_PRINT_ON
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
	void *user_data
) {
 
#ifdef DEBUG_PRINT_ON
  AnnounceStartBlock("Implicit RHS");
#endif

  // model time
  Time timeT = Time(0,0,0,time,0); // only works for integers

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

#ifdef DEBUG_PRINT_ON
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

#ifdef DEBUG_PRINT_ON
  AnnounceStartBlock("Full RHS");
#endif

  // model time (never used in RHS)
  Time timeT = Time(0,0,0,0,0); // only works for integers

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

#ifdef DEBUG_PRINT_ON
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
      iPorder = 1;
      
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
      iPorder = 1;
      
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
      iPorder = 1;
      
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
      iPorder = 1;
      
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

void TimestepSchemeARKode::ARKode_Test_Step(bool fFirstStep, 
					    double dDeltaT) {

  int iIdx_initial     = 20;
  int iIdx_arkode      = 21;
  int iIdx_tempest     = 22;
  int iIdx_arkode_like = 23;
  int iIdx_diff1       = 24;
  int iIdx_diff2       = 25; 
  int iIdx_diff3       = 26; 

  if (!m_fFullyExplicit) {
    _EXCEPTIONT("ARKode_Test_Step() ERROR: Not running fully explicitly");
  }

  Announce("ARKode_Test_Step()");

  // counter for function calls, controls output at end of function
  static int iCallCount = 0;

  // inrement function call counter
  iCallCount++;

  // dummy time variable
  Time time = Time(0,0,0,0,0);
  
  // Get a copy of the grid
  Grid * pGrid = m_model.GetGrid();
  
  // Take first step with ARKode and Tempest
  if (fFirstStep) {
   
    // Get a copy of the HorizontalDynamics
    HorizontalDynamics * pHorizontalDynamics =
      m_model.GetHorizontalDynamics();
    
    // Get a copy of the VerticalDynamics
    VerticalDynamics * pVerticalDynamics =
      m_model.GetVerticalDynamics();
    
    // Save Initial State
    pGrid->CopyData(0, iIdx_initial, DataType_State);
    pGrid->CopyData(0, iIdx_initial, DataType_Tracers);
    
    AnnounceStartBlock("Computing ARKode Step");
    
    // Current model time
    Time timeCurrentT = m_model.GetCurrentTime();
    double dCurrentT  = timeCurrentT.GetSeconds();
       
    Time timeNextT  = m_model.GetCurrentTime();
    Time timeDeltaT = m_model.GetDeltaT();
    timeNextT += timeDeltaT;
    
    double dNextT;
    dNextT = timeNextT.GetSeconds();
    
    int ierr;
    ierr = ARKode(ARKodeMem, dNextT, m_Y, &dCurrentT, ARK_ONE_STEP);
    
    if (ierr < 0) _EXCEPTION1("ERROR: ARKode, ierr = %i",ierr);

    // Exchange
    // pGrid->PostProcessSubstage(0, DataType_State);
    // pGrid->PostProcessSubstage(0, DataType_Tracers);
    
    // Save ARKode Solution
    pGrid->CopyData(0, iIdx_arkode, DataType_State);
    pGrid->CopyData(0, iIdx_arkode, DataType_Tracers);

    AnnounceEndBlock("Done");

    AnnounceStartBlock("Computing Tempest Step");
    
    // Reset Initial State
    pGrid->CopyData(iIdx_initial, 0, DataType_State);
    pGrid->CopyData(iIdx_initial, 0, DataType_Tracers);
    
    if (m_iSetButcherTable == 0) {

      // Compare to Forward Euler
      pGrid->CopyData(0, 4, DataType_State);
      pGrid->CopyData(0, 4, DataType_Tracers);

      pHorizontalDynamics->StepExplicit(0, 4, time, dDeltaT);
      pVerticalDynamics->StepExplicit(0, 4, time, dDeltaT);

      // pGrid->PostProcessSubstage(4, DataType_State);
      // pGrid->PostProcessSubstage(4, DataType_Tracers);

    } else if (m_iSetButcherTable == 1) { 

      // Compare to KGU(6,3) ERK method
      DataArray1D<double> dKinnmarkGrayUllrichCombination;
      
      dKinnmarkGrayUllrichCombination.Allocate(5);
      
      dKinnmarkGrayUllrichCombination[0] = - 1.0 / 4.0;
      dKinnmarkGrayUllrichCombination[1] =   5.0 / 4.0;
      dKinnmarkGrayUllrichCombination[2] =   0.0;
      dKinnmarkGrayUllrichCombination[3] =   0.0;
      dKinnmarkGrayUllrichCombination[4] =   0.0;
      
      // Z1 = Z0 = Yn
      pGrid->CopyData(0, 1, DataType_State);
      pGrid->CopyData(0, 1, DataType_Tracers);
      
      // Z1 = Z1 + 1/5 F(Z0) 
      pHorizontalDynamics->StepExplicit(0, 1, time, dDeltaT / 5.0);
      pVerticalDynamics->StepExplicit(0, 1, time, dDeltaT / 5.0);
      
      // pGrid->PostProcessSubstage(1, DataType_State);
      // pGrid->PostProcessSubstage(1, DataType_Tracers);
      
      //-------------------------------------------------------------------------
      
      // Z2 = Z0
      pGrid->CopyData(0, 2, DataType_State);
      pGrid->CopyData(0, 2, DataType_Tracers);
      
      // Z2 = Z2 + 1/5 F(Z1)
      pHorizontalDynamics->StepExplicit(1, 2, time, dDeltaT / 5.0);
      pVerticalDynamics->StepExplicit(1, 2, time, dDeltaT / 5.0);
      
      // pGrid->PostProcessSubstage(2, DataType_State);
      // pGrid->PostProcessSubstage(2, DataType_Tracers);
      
      //-------------------------------------------------------------------------
      
      // Z3 = Z0
      pGrid->CopyData(0, 3, DataType_State);
      pGrid->CopyData(0, 3, DataType_Tracers);
      
      // Z3 = Z3 + 1/3 F(Z2)
      pHorizontalDynamics->StepExplicit(2, 3, time, dDeltaT / 3.0);
      pVerticalDynamics->StepExplicit(2, 3, time, dDeltaT / 3.0);
      
      // pGrid->PostProcessSubstage(3, DataType_State);
      // pGrid->PostProcessSubstage(3, DataType_Tracers);
      
      //-------------------------------------------------------------------------
      
      // Z4 = Z0 (save over Z2)
      pGrid->CopyData(0, 2, DataType_State);
      pGrid->CopyData(0, 2, DataType_Tracers);
      
      // Z4 = Z4 + 2/3 F(Z3)
      pHorizontalDynamics->StepExplicit(3, 2, time, 2.0 * dDeltaT / 3.0);
      pVerticalDynamics->StepExplicit(3, 2, time, 2.0 * dDeltaT / 3.0);
      
      // pGrid->PostProcessSubstage(2, DataType_State);
      // pGrid->PostProcessSubstage(2, DataType_Tracers);
      
      //-------------------------------------------------------------------------
      
      // Z* = -1/4 Z0 + 5/4 Z1
      pGrid->LinearCombineData(dKinnmarkGrayUllrichCombination, 4, DataType_State);
      pGrid->LinearCombineData(dKinnmarkGrayUllrichCombination, 4, DataType_Tracers);
      
      // Yn+1 = Z5 = Z* + 3/4 F(Z4)
      pHorizontalDynamics->StepExplicit(2, 4, time, 3.0 * dDeltaT / 4.0);
      pVerticalDynamics->StepExplicit(2, 4, time, 3.0 * dDeltaT / 4.0);
      
      // pGrid->PostProcessSubstage(4, DataType_State);
      // pGrid->PostProcessSubstage(4, DataType_Tracers); 
      
      //-------------------------------------------------------------------------
      
    } else {
      _EXCEPTIONT("ARKode_Test_Step() ERROR: Invalid ERK Butcher table vaule");
    }

    // Save Tempest solution
    pGrid->CopyData(4, iIdx_tempest, DataType_State);
    pGrid->CopyData(4, iIdx_tempest, DataType_Tracers);

    AnnounceEndBlock("Done");

    AnnounceStartBlock("ARKode-Like Tempest Step");
    
    // Reset Initial State
    pGrid->CopyData(iIdx_initial, 0, DataType_State);
    pGrid->CopyData(iIdx_initial, 0, DataType_Tracers);
    
    if (m_iSetButcherTable == 0) {

      // Compare to Forward Euler

      // Zero vector to store RHS
      pGrid->ZeroData(1, DataType_State);
      pGrid->ZeroData(1, DataType_Tracers);

      // evaluate RHS F(Z0) = F(Yn)
      pHorizontalDynamics->StepExplicit(0, 1, time, 1.0);
      pVerticalDynamics->StepExplicit(0, 1, time, 1.0);

      // compute new solution 
      pGrid->LinearSumData(1.0, 0, dDeltaT, 1, 2);

      // save solution
      pGrid->CopyData(2, iIdx_arkode_like, DataType_State);
      pGrid->CopyData(2, iIdx_arkode_like, DataType_Tracers);

    } else if (m_iSetButcherTable == 1) { 

      // scaling factor for time steps
      double dC;
      
      // F(Z0)
      pGrid->ZeroData(1, DataType_State);
      pGrid->ZeroData(1, DataType_Tracers);
      
      pHorizontalDynamics->StepExplicit(0, 1, time, 1.0);
      pVerticalDynamics->StepExplicit(0, 1, time, 1.0);

      // Z1
      dC = 0.2 * dDeltaT;
      pGrid->LinearSumData(1.0, 0, dC, 1, 2);

      // F(Z1)
      pGrid->ZeroData(3, DataType_State);
      pGrid->ZeroData(3, DataType_Tracers);
      
      pHorizontalDynamics->StepExplicit(2, 3, time, 1.0);
      pVerticalDynamics->StepExplicit(2, 3, time, 1.0);

      // Z2
      dC = 0.2 * dDeltaT;
      pGrid->LinearSumData(1.0, 0, dC, 3, 2);

      // F(Z2)
      pGrid->ZeroData(3, DataType_State);
      pGrid->ZeroData(3, DataType_Tracers);
      
      pHorizontalDynamics->StepExplicit(2, 3, time, 1.0);
      pVerticalDynamics->StepExplicit(2, 3, time, 1.0);

      // Z3
      dC = dDeltaT / 3.0;
      pGrid->LinearSumData(1.0, 0, dC, 3, 2);

      // F(Z3)
      pGrid->ZeroData(3, DataType_State);
      pGrid->ZeroData(3, DataType_Tracers);
      
      pHorizontalDynamics->StepExplicit(2, 3, time, 1.0);
      pVerticalDynamics->StepExplicit(2, 3, time, 1.0);

      // Z4
      dC = 2.0 * dDeltaT / 3.0;
      pGrid->LinearSumData(1.0, 0, dC, 3, 2);

      // F(Z4)
      pGrid->ZeroData(3, DataType_State);
      pGrid->ZeroData(3, DataType_Tracers);
      
      pHorizontalDynamics->StepExplicit(2, 3, time, 1.0);
      pVerticalDynamics->StepExplicit(2, 3, time, 1.0);

      // Z5 = Yn+1
      pGrid->LinearSumData(0.25, 1, 0.75, 3, 2);
      pGrid->LinearSumData(1.0, 0, dDeltaT, 2, 3);

      // save solution
      pGrid->CopyData(3, iIdx_arkode_like, DataType_State);
      pGrid->CopyData(3, iIdx_arkode_like, DataType_Tracers);
    }

    // Compute difference between ARKode and Tempest solutions
    pGrid->LinearSumData(1.0, iIdx_arkode, -1.0, iIdx_tempest, iIdx_diff1);

    // Compute difference between ARKode-Like and Tempest solutions
    pGrid->LinearSumData(1.0, iIdx_arkode_like, -1.0, iIdx_tempest, iIdx_diff2);

    // Difference between initial and saved initial
    pGrid->LinearSumData(1.0, 0, -1.0, iIdx_initial, iIdx_diff3);
  }

  // Output 
  if (iCallCount == 1) {
    // output ARKode
    Announce("Output ARKode Step");
    pGrid->CopyData(iIdx_arkode, 0, DataType_State);
    pGrid->CopyData(iIdx_arkode, 0, DataType_Tracers);

  } else if (iCallCount == 2) {
    // output Tempest
    Announce("Output Tempest Step");
    pGrid->CopyData(iIdx_tempest, 0, DataType_State);
    pGrid->CopyData(iIdx_tempest, 0, DataType_Tracers);

  } else if (iCallCount == 3) {
    // output diff
    Announce("Output ARKode-Like");
    pGrid->CopyData(iIdx_arkode_like, 0, DataType_State);
    pGrid->CopyData(iIdx_arkode_like, 0, DataType_Tracers);   

  } else if (iCallCount == 4) {
    // output diff
    Announce("Output Diff1");
    pGrid->CopyData(iIdx_diff1, 0, DataType_State);
    pGrid->CopyData(iIdx_diff1, 0, DataType_Tracers);   

  } else if (iCallCount == 5) {
    // output diff
    Announce("Output Diff2");
    pGrid->CopyData(iIdx_diff2, 0, DataType_State);
    pGrid->CopyData(iIdx_diff2, 0, DataType_Tracers);   

  } else if (iCallCount == 6) {
    // output diff
    Announce("OutputDiff3"); 
    pGrid->CopyData(iIdx_diff3, 0, DataType_State);
    pGrid->CopyData(iIdx_diff3, 0, DataType_Tracers);
 
 } else {
    Announce("End of Test Step"); 
  }
}

///////////////////////////////////////////////////////////////////////////////

void TimestepSchemeARKode::Step_Like_ARKode(bool fFirstStep, 
					    bool fLastStep,
					    const Time & time,
					    double dDeltaT) {
 
  // Get a copy of the grid
  Grid * pGrid = m_model.GetGrid();

  // Get a copy of the HorizontalDynamics
  HorizontalDynamics * pHorizontalDynamics =
    m_model.GetHorizontalDynamics();
  
  // Get a copy of the VerticalDynamics
  VerticalDynamics * pVerticalDynamics =
    m_model.GetVerticalDynamics();
  
  // scaling factor for time steps
  double dC;
  
  // Evalutate RHS(Z0), Z0 = Yn
  pGrid->ZeroData(1, DataType_State);
  pGrid->ZeroData(1, DataType_Tracers);
  
  pHorizontalDynamics->StepExplicit(0, 1, time, 1.0);
  pVerticalDynamics->StepExplicit(0, 1, time, 1.0);
  
  // Compute Z1 = Yn + 1/5 dt RHS(Z0)
  dC = 0.2 * dDeltaT;
  pGrid->LinearSumData(1.0, 0, dC, 1, 2);

  pGrid->PostProcessSubstage(2, DataType_State);
  pGrid->PostProcessSubstage(2, DataType_Tracers); 
  
  // Evaluate RHS(Z1)
  pGrid->ZeroData(3, DataType_State);
  pGrid->ZeroData(3, DataType_Tracers);
  
  pHorizontalDynamics->StepExplicit(2, 3, time, 1.0);
  pVerticalDynamics->StepExplicit(2, 3, time, 1.0);
  
  // Compute Z2 = Yn + 1/5 dt RHS(Z1)
  dC = 0.2 * dDeltaT;
  pGrid->LinearSumData(1.0, 0, dC, 3, 2);

  pGrid->PostProcessSubstage(2, DataType_State);
  pGrid->PostProcessSubstage(2, DataType_Tracers); 
  
  // Evaluate RHS(Z2)
  pGrid->ZeroData(3, DataType_State);
  pGrid->ZeroData(3, DataType_Tracers);
  
  pHorizontalDynamics->StepExplicit(2, 3, time, 1.0);
  pVerticalDynamics->StepExplicit(2, 3, time, 1.0);
  
  // Compute Z3 = Yn + 1/3 dt RHS(Z2)
  dC = dDeltaT / 3.0;
  pGrid->LinearSumData(1.0, 0, dC, 3, 2);

  pGrid->PostProcessSubstage(2, DataType_State);
  pGrid->PostProcessSubstage(2, DataType_Tracers); 
  
  // Evaluate RHS(Z3)
  pGrid->ZeroData(3, DataType_State);
  pGrid->ZeroData(3, DataType_Tracers);
  
  pHorizontalDynamics->StepExplicit(2, 3, time, 1.0);
  pVerticalDynamics->StepExplicit(2, 3, time, 1.0);
  
  // Compute Z4 = Yn + 2/3 dt RHS(Z3)
  dC = 2.0 * dDeltaT / 3.0;
  pGrid->LinearSumData(1.0, 0, dC, 3, 2);

  pGrid->PostProcessSubstage(2, DataType_State);
  pGrid->PostProcessSubstage(2, DataType_Tracers); 
  
  // Evaluate RHS(Z4)
  pGrid->ZeroData(3, DataType_State);
  pGrid->ZeroData(3, DataType_Tracers);
  
  pHorizontalDynamics->StepExplicit(2, 3, time, 1.0);
  pVerticalDynamics->StepExplicit(2, 3, time, 1.0);
  
  // Compute Z5 = Yn+1 = Yn + 1/4 RHS(Z0) + 3/4 RHS(Z4)
  pGrid->LinearSumData(0.25, 1, 0.75, 3, 2);
  pGrid->LinearSumData(1.0, 0, dDeltaT, 2, 1);

  pGrid->PostProcessSubstage(1, DataType_State);
  pGrid->PostProcessSubstage(1, DataType_Tracers); 

  // Apply hyperdiffusion (initial, update, temp)
  pGrid->CopyData(1, 2, DataType_State);
  pGrid->CopyData(1, 2, DataType_Tracers);

  pHorizontalDynamics->StepAfterSubCycle(1, 2, 3, time, dDeltaT);

  pGrid->CopyData(2, 0, DataType_State);
  pGrid->CopyData(2, 0, DataType_Tracers);
}

#endif
