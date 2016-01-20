///////////////////////////////////////////////////////////////////////////////
///
///	\file    RHSAnalysis.h
///	\author  Daniel R. Reynolds
///	\version January 13, 2016
///
///	<remarks>
///		Copyright 2016 Daniel R. Reynolds, Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _RHSANALYSIS_H_
#define _RHSANALYSIS_H_

#include "TimeObj.h"
#include "DataArray1D.h"
#include "Model.h"
#include "TimestepScheme.h"
#include "Exception.h"


///////////////////////////////////////////////////////////////////////////////

//class Model;
//class Time;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Interface for RHS function analysis
///	</summary>
class RHSAnalysis {

 public:
  ///	<summary>
  ///	        Implicit/explicit RHS function choices
  ///	</sumamry>
  enum ImExRHS {
    HEVE,
    HEVI
  };

 protected:
  ///	<summary>
  ///		Constructor.
  ///	</summary>
  RHSAnalysis(TimestepScheme & TS,
	      const Time StartTime,
	      Model & model,
	      ImExRHS ImEx)
    : m_TS(TS),
      m_model(model),
      m_ImEx(ImEx),
      m_StartTime(StartTime)
  { 
    // Only perform analysis for HEVE (for now)
    if (m_ImEx == HEVI) 
      _EXCEPTIONT("HEVI RHS analysis not supported (yet)");

    // Only perform analysis if TS has sufficient RHS storage
    m_nRHS = (TS.GetComponentDataInstances() > TS.GetTracerDataInstances()) ?
      TS.GetComponentDataInstances() : TS.GetTracerDataInstances();
    if (m_nRHS < 3)
      _EXCEPTIONT("TS contains insufficient RHS array storage for analysis");

    m_RHSSum.Allocate(m_nRHS);
  }

 public:
  ///	<summary>
  ///		Virtual destructor.
  ///	</summary>
  ~RHSAnalysis() { }

  ///	<summary>
  ///		Analysis routine (tests current time to determine if it should run).
  ///	</summary>
  void Analyze(const Time & time);

 private:
  ///	<summary>
  ///		Number of RHS arrays in m_TS object
  ///	</summary>
  int m_nRHS;

  ///	<summary>
  ///		Reference to the model.
  ///	</summary>
  Model & m_model;

  ///	<summary>
  ///		Reference to the time-stepping object.
  ///	</summary>
  TimestepScheme & m_TS;

  ///	<summary>
  ///		ImEx RHS functions to use.
  ///	</summary>
  ImExRHS m_ImEx;

  ///	<summary>
  ///		Start time for running analysis.
  ///	</summary>
  Time m_StartTime;

  ///	<summary>
  ///		Sum of RHS terms
  ///	</summary>
  DataArray1D<double> m_RHSSum;

};

///////////////////////////////////////////////////////////////////////////////

#endif

