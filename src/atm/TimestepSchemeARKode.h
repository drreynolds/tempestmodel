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
#include "arkode/arkode_arkstep.h"
#include "arkode/arkode_ls.h"
#include "sunlinsol/sunlinsol_spgmr.h"
#include "arkode/arkode_butcher.h"
#include "arkode/arkode_butcher_dirk.h"
#include "arkode/arkode_butcher_erk.h"
#include "sundials/sundials_config.h"
#include "sundials/sundials_linearsolver.h"
#include "sunmatrix/sunmatrix_dense.h"

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
  int    ErrController;
  int    AAFPAccelVec;
  int    NonlinIters;
  int    LinIters;
  int    Predictor;
  double VAtol_vel;
  double VAtol_rho;
  double VAtol_theta;
  bool   UsePreconditioning;
  bool   ColumnSolver;
  int    ARKodeButcherTable;
  std::string ButcherTable;
  std::string StepOut;
  bool   WriteDiagnostics;
  bool   FullyImplicit;
  Time   OutputTime;
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

	///	<summary>
	///		Assign component-wise tolerances used in the function ARKodeSVTolerances
	///	</summary>
	void AssignComponentWiseTolerances();

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
	///		Tempest NVector of tolerances;
	///	</summary>
	N_Vector m_T;

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
	///		Name of file to output time step profile.
	///	</summary>
	std::string m_strStepOut;

	///	<summary>
	///		File the time stepping profile is written to.
	///	</summary>
	FILE *m_fStep_Profile;

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
	///		ARKode flag for which adapt method.
	///	</summary>
	int m_iErrController;

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
	///		Specifies method for prediciting implicit solutions.
	///	</summary>
	int m_iPredictor;

	///	<summary>
	///		Vector component of Atol for velocity.
	///	</summary>
	double m_dVAtol_vel;

	///	<summary>
	///		Vector component of Atol for theta.
	///	</summary>
	double m_dVAtol_theta;

	///	<summary>
	///		Vector component of Atol for rho (density).
	///	</summary>
	double m_dVAtol_rho;

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

	///	<summary>
	///		ARKode variable for frequency of output
	///	</summary>
	Time m_tOutT;

	double dNextOut;

	Time m_tNextOutT;

	///	<summary>
	///		Used to send dynamic dt from ARKode to
	///		hyperviscosity stepping.
	///	</summary>
	double m_iTemp_dt;
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
	void *user_data
);

///	<summary>
///		Check function return value
///	</summary>
static int check_flag(
	void *flagvalue,
	const char *funcname,
	int opt
);

///	<summary>
///		Functions for replacing GMRES solver with Tempest column-wise linear solver
///	</summary>
static int ARKodeLinSysFn(realtype t, N_Vector y,
			N_Vector fy, SUNMatrix A,
			SUNMatrix M, booleantype jok,
			booleantype *jcur, realtype gamma,
			void *user_data, N_Vector tmp1,
			N_Vector tmp2, N_Vector tmp3
);

SUNMatrix SUNMatrix_Tempest();

SUNMatrix SUNMatClone_Tempest(SUNMatrix M);

SUNMatrix_ID SUNMatGetID_Tempest(SUNMatrix M);

SUNLinearSolver SUNLinSol_Tempest(void* ark_mem);

int ARKodeColumnLSolve(
	SUNLinearSolver S,
	SUNMatrix A,
	N_Vector b,
	N_Vector ycur,
	realtype tol
);

SUNLinearSolver_Type ARKodeColumnLType(SUNLinearSolver S);

int ARKodeColumnLFree(SUNLinearSolver S);

///////////////////////////////////////////////////////////////////////////////

#endif

#endif


