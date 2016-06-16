///////////////////////////////////////////////////////////////////////////////
///
///	\file    TimestepSchemeARKode.h
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

#ifndef _TIMESTEPSCHEMEARKODE_H_
#define _TIMESTEPSCHEMEARKODE_H_
///////////////////////////////////////////////////////////////////////////////

// require SUNDIALS for compilation
#ifdef USE_SUNDIALS

#include "TimestepScheme.h"
#include "Exception.h"

#include "TempestNVector.h"
#include "arkode/arkode.h"
#include "arkode/arkode_impl.h"
#include "arkode/arkode_spgmr.h"

class Model;
class Time;

///////////////////////////////////////////////////////////////////////////////

struct ARKodeCommandLineVariables {
  int    nvectors;
  double rtol;
  double atol;
  bool   FullyExplicit;
  bool   DynamicStepSize;
  bool   AAFP;
  int    AAFPAccelVec;
  int    NonlinIters;
  int    LinIters;
  bool   UsePreconditioning;
  bool   ColumnSolver;
  int    ARKodeButcherTable;
  std::string ButcherTable;
  bool   WriteDiagnostics;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Time stepping methods provided by the ARKode library.
///	</summary>
class TimestepSchemeARKode : public TimestepScheme {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	TimestepSchemeARKode(
		Model & model,
		ARKodeCommandLineVariables & ARKodeVars
	);

private:
	///	<summary>
	///		Initializer.  Called prior to model execution.
	///	</summary>
	virtual void Initialize();

	///	<summary>
	///		Set user supplied Butcher table.
	///	</summary>
	void SetButcherTable();

public:
	///	<summary>
	///		Get the number of component data instances (default is 50).
	///	</summary>
	virtual int GetComponentDataInstances() const {
		return m_iNVectors;
	}

	///	<summary>
	///		Get the number of tracer data instances (default is 50).
	///	</summary>
	virtual int GetTracerDataInstances() const {
		return m_iNVectors;
	}

	///	<summary>
	///		Get the most recent step size used with dynamic 
	///             timestepping.
	///	</summary>
	virtual double GetDynamicDeltaT() const {
		return m_dDynamicDeltaT;
	}

protected:
	///	<summary>
	///		Perform one time step.
	///	</summary>
	virtual void Step(
		bool fFirstStep,
		bool fLastStep,
		const Time & time,
		double dDeltaT
	);

private:
	///	<summary>
	///		ARKode memory structure.
	///	</summary>
	static void * arkode_mem;

	///	<summary>
	///		Tempest NVector state vector.
	///	</summary>
	N_Vector m_Y;

	///	<summary>
	///		Number of NVectors (default is 50).
	///	</summary>
	int m_iNVectors;

	///	<summary>
	///		ARKode Butcher table number.
	///	</summary>
	int m_iARKodeButcherTable;

	///	<summary>
	///		User supplied Butcher table name.
	///	</summary>
	std::string m_strButcherTable;

	///	<summary>
	///		ARKode absolute tolerance.
	///	</summary>
	double m_dAbsTol;

	///	<summary>
	///		ARKode relative tolerance.
	///	</summary>
	double m_dRelTol;

	///	<summary>
	///		ARKode flag for fully explicit integration.
	///	</summary>
	bool m_fFullyExplicit;

	///	<summary>
	///		ARKode flag for fully implicit integration.
	///	</summary>
	bool m_fFullyImplicit;

	///	<summary>
	///		ARKode flag to used adaptive step sizes.
	///	</summary>
	bool m_fDynamicStepSize;

	///	<summary>
	///		Most recent step size in ARKode when using dynamic
	///             timestepping.
	///	</summary>
	double m_dDynamicDeltaT;

	///	<summary>
	///		ARKode flag for Anderson Accelerated Fixed Point solver.
	///	</summary>
	bool m_fAAFP;

	///	<summary>
	///		Max number of acceleration vectors.
	///	</summary>
	int m_iAAFPAccelVec;

	///	<summary>
	///		Max number of nonlinear iterations.
	///	</summary>
	int m_iNonlinIters;

	///	<summary>
	///		Max number of linear iterations.
	///	</summary>
	int m_iLinIters;

	///	<summary>
	///		ARKode flag to write diagnostics file.
	///	</summary>
	bool m_fWriteDiagnostics;

	///	<summary>
	///		ARKode flag to enable preconditioning.
	///	</summary>
	bool m_fUsePreconditioning;

	///	<summary>
	///		ARKode flag to use Tempest column solver for 
	///             linear systems (instead of GMRES).
	///	</summary>
	bool m_fColumnSolver;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Wrapper function for ARKode to compute explicit RHS.
///	</summary>
static int ARKodeExplicitRHS(
	realtype time, 
	N_Vector Y, 
	N_Vector Ydot, 
	void *user_data
);

///	<summary>
///		Wrapper function for ARKode to compute implicit RHS.
///	</summary>
static int ARKodeImplicitRHS(
	realtype time, 
	N_Vector Y, 
	N_Vector Ydot, 
	void *user_data
);

///	<summary>
///		Wrapper function for ARKode to compute the full RHS.
///	</summary>
static int ARKodeFullRHS(
	realtype time, 
	N_Vector Y, 
	N_Vector Ydot, 
	void *user_data
);

///	<summary>
///		Function to perform operations on ARKode step solution.
///	</summary>
static int ARKodePostProcessStep(
	realtype time, 
	N_Vector Y, 
	void * user_data
);

///	<summary>
///		Function to perform columnwise preconditioner solve
///	</summary>
static int ARKodePreconditionerSolve(
	realtype time, 
	N_Vector Y, 
	N_Vector FY, 
	N_Vector R,
	N_Vector Z,
	realtype gamma,
	realtype delta,
	int lr,
	void *user_data,
	N_Vector TMP
);

///	<summary>
///		Functions for replacing GMRES solver with Tempest column-wise linear solver
///	</summary>
int ARKodeColumnLInit(ARKodeMem ark_mem);
int ARKodeColumnLSetup(
        ARKodeMem ark_mem,
        int convfail,
        N_Vector ypred,
        N_Vector fpred,
        booleantype *jcurPtr,
        N_Vector vtemp1,
        N_Vector vtemp2,
        N_Vector vtemp3
);
int ARKodeColumnLSolve(
        ARKodeMem ark_mem,
        N_Vector b,
        N_Vector weight,
        N_Vector ycur,
        N_Vector fcur
);
void ARKodeColumnLFree(ARKodeMem ark_mem);

///////////////////////////////////////////////////////////////////////////////

#endif

#endif


