///////////////////////////////////////////////////////////////////////////////
///
///	\file    TimestepSchemeBHR553.h
///	\author  Paul Ullrich
///	\version April 22, 2014
///
///	<remarks>
///		Copyright 2000-2010 Paul Ullrich, Jorge Guerra
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _TIMESTEPSCHEMEBHR553_H_
#define _TIMESTEPSCHEMEBHR553_H_

#include "TimestepScheme.h"
#include "Exception.h"
#include "DataVector.h"

///////////////////////////////////////////////////////////////////////////////

class Model;
class Time;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Adaptive Fourth-Order Runge-Kutta time stepping.
///	</summary>
class TimestepSchemeBHR553 : public TimestepScheme {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	TimestepSchemeBHR553(
		Model & model
	);

public:
	///	<summary>
	///		Get the number of component data instances.
	///	</summary>
	virtual int GetComponentDataInstances() const {
		return 13;
	}

	///	<summary>
	///		Get the number of tracer data instances.
	///	</summary>
	virtual int GetTracerDataInstances() const {
		return 13;
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
	///		BHR553 parameter gamma
	///	</summary>
	static const double m_dgamma;
	///	<summary>
	///		BHR553 parameter delta
	///	</summary>
	static const double m_ddelta;

    ///	<summary>
	///		Coefficients for the time increment BHR553.
	///	</summary>
	static const double m_dTimeCf[4];

    ///	<summary>
	///		Coefficients for the explicit BHR553.
	///	</summary>
	static const double m_dExpCf[6][6];

	///	<summary>
	///		Coefficients for the explicit BHR553.
	///	</summary>
	static const double m_dImpCf[6][6];

	static const double m_c2;
	static const double m_c3;
	static const double m_c4;
	static const double m_b1;
	static const double m_b2;
	static const double m_b3;
	static const double m_b4;
	static const double m_at31;
	static const double m_at32;
	static const double m_at41;
	static const double m_at42;
	static const double m_at43;
	static const double m_at51;
	static const double m_at52;
	static const double m_at53;
	static const double m_at54;
	static const double m_a31;
	static const double m_a32;
	static const double m_a41;
	static const double m_a42;
	static const double m_a43;

	///	<summary>
	///		K0 vector combination (Implicit)
	///	</summary>
	DataVector<double> m_dK0Combo;

	///	<summary>
	///		K1 vector combination (Implicit)
	///	</summary>
	DataVector<double> m_dK1Combo;

	///	<summary>
	///		Kh1 vector combination (Explicit)
	///	</summary>
	DataVector<double> m_dKh1Combo;

	///	<summary>
	///		Explicit evaluation at the 2nd substage
	///	</summary>
	DataVector<double> m_du1fCombo;

	///	<summary>
	///		Explicit evaluation at the 2nd substage
	///	</summary>
	DataVector<double> m_du2fCombo;

	///	<summary>
	///		K2 vector combination (Implicit)
	///	</summary>
	DataVector<double> m_dK2Combo;

	///	<summary>
	///		Kh2 vector combination (Explicit)
	///	</summary>
	DataVector<double> m_dKh2Combo;

	///	<summary>
	///		Explicit evaluation at the 3rd substage
	///	</summary>
	DataVector<double> m_du3fCombo;

	///	<summary>
	///		K3 vector combination (Implicit)
	///	</summary>
	DataVector<double> m_dK3Combo;

	///	<summary>
	///		Kh3 vector combination (Explicit)
	///	</summary>
	DataVector<double> m_dKh3Combo;

	///	<summary>
	///		Explicit evaluation at the 4th substage
	///	</summary>
	DataVector<double> m_du4fCombo;

	///	<summary>
	///		K4 vector combination (Implicit)
	///	</summary>
	DataVector<double> m_dK4Combo;

	///	<summary>
	///		Kh4 vector combination (Explicit)
	///	</summary>
	DataVector<double> m_dKh4Combo;

	///	<summary>
	///		Explicit evaluation at the 5th substage
	///	</summary>
	DataVector<double> m_dufCombo;

	///	<summary>
	///		Kh5 vector combination (Explicit)
	///	</summary>
	DataVector<double> m_dKh5Combo;
};

///////////////////////////////////////////////////////////////////////////////

#endif


