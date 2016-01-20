///////////////////////////////////////////////////////////////////////////////
///
///	\file    RHSAnalysis.cpp
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

#include "RHSAnalysis.h"
#include "Grid.h"
#include "Announce.h"
#include "HorizontalDynamics.h"
#include "VerticalDynamics.h"


///////////////////////////////////////////////////////////////////////////////

void RHSAnalysis::Analyze(const Time & time) {

  // Only perform analysis if time >= m_StartTime
  if (time < m_StartTime)  return;

  // Get processor rank
  int nRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

  // Equation set
  const EquationSet & eqn = m_model.GetEquationSet();
  
  // Get a copy of the grid
  Grid * pGrid = m_model.GetGrid();

  // Get a copy of the HorizontalDynamics
  HorizontalDynamics * pHorizontalDynamics = m_model.GetHorizontalDynamics();

  // Get a copy of the VerticalDynamics
  VerticalDynamics * pVerticalDynamics = m_model.GetVerticalDynamics();

  // Compute the horizontal dynamics using field 0, current time, 
  // and time step size of 1; store in field 1
  pGrid->ZeroData(1, DataType_State);       // zero-out state and tracer 
  pGrid->ZeroData(1, DataType_Tracers);     //    data in field 1
  pHorizontalDynamics->StepExplicit(0, 1, time, 1.0);

  // Compute norms of RHS and output
  DataArray1D<double> dRHS;
  pGrid->Checksum(DataType_State, dRHS, 1, ChecksumType_Linf);
  // pGrid->Checksum(DataType_State, dRHS, 1, ChecksumType_L1);
  // pGrid->Checksum(DataType_State, dRHS, 1, ChecksumType_L2);
  if (nRank == 0) {
    for (int c=0; c<eqn.GetComponents(); c++) {
      Announce("HorizontalDynamics RHS (%s): %1.15e",
	       eqn.GetComponentShortName(c).c_str(), dRHS[c]);
    }
  }


  // Compute the vertical dynamics using field 0, current time, 
  // and time step size of 1; store in field 2
  pGrid->ZeroData(2, DataType_State);       // zero-out state and tracer 
  pGrid->ZeroData(2, DataType_Tracers);     //    data in field 2
  pVerticalDynamics->StepExplicit(0, 2, time, 1.0);

  // Compute norms of RHS and output
  pGrid->Checksum(DataType_State, dRHS, 1, ChecksumType_Linf);
  if (nRank == 0) {
    for (int c=0; c<eqn.GetComponents(); c++) {
      Announce("VerticalDynamics RHS (%s): %1.15e",
	       eqn.GetComponentShortName(c).c_str(), dRHS[c]);
    }
  }


  // Compute sum of horizontal and vertical RHS; store in field 3
  for (int i=0; i<m_nRHS; i++)  m_RHSSum[i] = 0.0;
  m_RHSSum[1] = 1.0;
  m_RHSSum[2] = 1.0;
  pGrid->LinearCombineData(m_RHSSum, 3, DataType_State);

  // Compute norms of RHS and output
  pGrid->Checksum(DataType_State, dRHS, 1, ChecksumType_Linf);
  if (nRank == 0) {
    for (int c=0; c<eqn.GetComponents(); c++) {
      Announce("Full RHS (%s): %1.15e",
	       eqn.GetComponentShortName(c).c_str(), dRHS[c]);
    }
  }
  

}

///////////////////////////////////////////////////////////////////////////////

