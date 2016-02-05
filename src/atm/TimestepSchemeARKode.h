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
#include "arkode/arkode_spgmr.h"

class Model;
class Time;

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
		Model & model
	);

private:
	///	<summary>
	///		Initializer.  Called prior to model execution.
	///	</summary>
	virtual void Initialize();

public:
	///	<summary>
	///		Get the number of component data instances.
	///	</summary>
	virtual int GetComponentDataInstances() const {
	        _EXCEPTIONT("Not implemented");		
		//return 7;
	}

	///	<summary>
	///		Get the number of tracer data instances.
	///	</summary>
	virtual int GetTracerDataInstances() const {
		_EXCEPTIONT("Not implemented");
		//return 7;
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
	static void * ARKodeMem;

	///	<summary>
	///		Tempest NVector state vector.
	///	</summary>
	N_Vector m_Y;

};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Wrapper funciton for ARKode to compute explicit RHS.
///	</summary>
static int ARKodeExplicitRHS(
	realtype time, 
	N_Vector Y, 
	N_Vector Ydot, 
	void *user_data
);

///	<summary>
///		Wrapper funciton for ARKode to compute implicit RHS.
///	</summary>
static int ARKodeImplicitRHS(
	realtype time, 
	N_Vector Y, 
	N_Vector Ydot, 
	void *user_data
);

///////////////////////////////////////////////////////////////////////////////

#endif

#endif


