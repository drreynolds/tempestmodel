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

#define DEBUG_PRINT_ON

#include "TimestepSchemeARKode.h"
#include "Model.h"
#include "Grid.h"
#include "HorizontalDynamics.h"
#include "VerticalDynamics.h"
#include "Announce.h"

void * TimestepSchemeARKode::ARKodeMem = NULL;

///////////////////////////////////////////////////////////////////////////////

TimestepSchemeARKode::TimestepSchemeARKode(
        Model & model,
	ARKodeCommandLineVariables & ARKodeVars
) :
        TimestepScheme(model),
	m_iNVectors(ARKodeVars.nvectors),
	m_dRelTol(ARKodeVars.rtol),
	m_dAbsTol(ARKodeVars.atol)
{
        // Allocate ARKode memory
        ARKodeMem = ARKodeCreate();

	if (ARKodeMem == NULL)
	  _EXCEPTIONT("ERROR: ARKodeCreate returned NULL");		
}
					  
///////////////////////////////////////////////////////////////////////////////

void TimestepSchemeARKode::Initialize() {
  
  AnnounceStartBlock("Initializeing ARKode");

  // error flag
  int ierr = 0;

  // Get a copy of the grid
  Grid * pGrid = m_model.GetGrid();

  // Create a Tempest NVector for the model state
  m_Y = N_VNew_Tempest(*pGrid, m_model);

  if (m_Y == NULL) _EXCEPTIONT("ERROR: N_VNew_Tempest returned NULL");

  // Copy initial condition into current state vector
  int iY = NV_INDEX_TEMPEST(m_Y);

  pGrid->CopyData(0, iY, DataType_State);
  pGrid->CopyData(0, iY, DataType_Tracers);
 
  // Get current simulation time in seconds
  Time timeCurrentT = m_model.GetCurrentTime();
  double dCurrentT  = timeCurrentT.GetSeconds();

  // Initialize ARKode
  ierr = ARKodeInit(ARKodeMem, ARKodeExplicitRHS, ARKodeImplicitRHS, dCurrentT, m_Y);

  if (ierr < 0) _EXCEPTION1("ERROR: ARKodeInit, ierr = %i",ierr);		

  // Set diagnostics output file
  // FILE * pFile = stdout; 
  // ierr = ARKodeSetDiagnostics(ARKodeMem, pFile);
  // if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetDiagnostics, ierr = %i",ierr);

  // Set fixed step size in seconds
  Time timeDeltaT = m_model.GetDeltaT();
  double dDeltaT  = timeDeltaT.GetSeconds();

  ierr = ARKodeSetFixedStep(ARKodeMem, dDeltaT);

  if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetFixedStep, ierr = %i",ierr);
  
  // Specify tolerance
  ierr = ARKodeSStolerances(ARKodeMem, m_dRelTol, m_dAbsTol);

  if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSStolerances, ierr = %i",ierr);

  // Nonlinear Solver Settings

  // Linear Solver Settings
  ierr = ARKSpgmr(ARKodeMem, PREC_NONE, 0);

  if (ierr < 0) _EXCEPTION1("ERROR: ARKSpgmr, ierr = %i",ierr);

  AnnounceEndBlock("Done");
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

  // error flag
  int ierr = 0;

  // Adjust last time step size
  if (fLastStep) {
    ierr = ARKodeSetFixedStep(ARKodeMem, dDeltaT);

    if (ierr < 0) _EXCEPTION1("ERROR: ARKodeSetFixedStep, ierr = %i",ierr);
  }

  // Get a copy of the grid
  Grid * pGrid = m_model.GetGrid();

  // Copy current state in Tempest NVector
  int iY = NV_INDEX_TEMPEST(m_Y);

  pGrid->CopyData(0, iY, DataType_State);
  pGrid->CopyData(0, iY, DataType_Tracers);

  // Get current and next time values in seconds
  Time timeCurrentT = m_model.GetCurrentTime();
  double dCurrentT  = timeCurrentT.GetSeconds();

  Time timeNextT  = m_model.GetCurrentTime();
  Time timeDeltaT = m_model.GetDeltaT();
  timeNextT += timeDeltaT;

  double dNextT  = timeNextT.GetSeconds();
  
  // ARKode timestep
  ierr = ARKode(ARKodeMem, dNextT, m_Y, &dCurrentT, ARK_ONE_STEP);

  if (ierr < 0) _EXCEPTION1("ERROR: ARKode, ierr = %i",ierr);		

  pGrid->PostProcessSubstage(iY, DataType_State);
  pGrid->PostProcessSubstage(iY, DataType_Tracers);

  // Copy new state into current state position
  pGrid->CopyData(iY, 0, DataType_State);
  pGrid->CopyData(iY, 0, DataType_Tracers);

#ifdef DEBUG_PRINT_ON
  AnnounceEndBlock("Done");
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
  Time timeT = Time(0,0,0,time,0);

  // index of input data in registry
  int iY = NV_INDEX_TEMPEST(Y);

  // index of output data in registry
  int iYdot = NV_INDEX_TEMPEST(Ydot);

  // Get a copy of the grid
  Grid * pGrid = NV_GRID_TEMPEST(Y);

  // Get a copy of the model
  Model * pModel = NV_MODEL_TEMPEST(Y);

  // Get a copy of the HorizontalDynamics
  HorizontalDynamics * pHorizontalDynamics = pModel->GetHorizontalDynamics();

  // zero out iYdot
  pGrid->ZeroData(iYdot, DataType_State);
  pGrid->ZeroData(iYdot, DataType_Tracers);

  // Compute explicit RHS
  pHorizontalDynamics->StepExplicit(iY, iYdot, timeT, 1.0);

  pGrid->PostProcessSubstage(iYdot, DataType_State);
  pGrid->PostProcessSubstage(iYdot, DataType_Tracers);

#ifdef DEBUG_PRINT_ON
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
  Time timeT = Time(0,0,0,time,0);

  // index of input data in registry
  int iY = NV_INDEX_TEMPEST(Y);

  // index of output data in registry
  int iYdot = NV_INDEX_TEMPEST(Ydot);

  // Get a copy of the grid
  Grid * pGrid = NV_GRID_TEMPEST(Y);

  // Get a copy of the model
  Model * pModel = NV_MODEL_TEMPEST(Y);

  // Get a copy of the VerticalDynamics
  VerticalDynamics * pVerticalDynamics = pModel->GetVerticalDynamics();

  // zero out iYdot
  pGrid->ZeroData(iYdot, DataType_State);
  pGrid->ZeroData(iYdot, DataType_Tracers);

  // Compute implicit RHS
  pVerticalDynamics->StepImplicitTermsExplicitly(iY, iYdot, timeT, 1.0);

  pGrid->PostProcessSubstage(iYdot, DataType_State);
  pGrid->PostProcessSubstage(iYdot, DataType_Tracers);

#ifdef DEBUG_PRINT_ON
  AnnounceEndBlock("Done");
#endif

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

#endif
