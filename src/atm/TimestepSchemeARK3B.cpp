///
///	\file    TimestepSchemeARK3B.cpp
///	\author  Paul Ullrich
///	\version May 29, 2014
///
///	<remarks>
///		Copyright 2000-2010 Paul Ullrich and Jorge Guerra
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "TimestepSchemeARK3B.h"
#include "Model.h"
#include "Grid.h"
#include "HorizontalDynamics.h"
#include "VerticalDynamics.h"

///////////////////////////////////////////////////////////////////////////////
// THESE COEFFICIENTS ARE COMPUTED FROM THE ORIGINAL TABLEAUX
// IMPLEMENTS Modified ARK(3,4,3) FROM ASCHER ET AL. 1997 PG. 9
// MODIFIED BY Sebantiano Boscarino
// Time step coefficients
const double TimestepSchemeARK3B::m_dTimeCf[3] = {0.4358665215, .7179332608, 1.};
// Implicit stage coefficients - Converted to U from Ascher et al. 1997 Pg 9

const double TimestepSchemeARK3B::m_dImpCf[4][4] = {
  {0.43586652150845899941, 0., 0., 0.},
  {0.28206673924577050030, 0.43586652150845899941, 0., 0.},
  {1.2084966491760100703, -0.6443631706844690697, 0.43586652150845899941, 0.},
  {1.2084966491760100703, -0.6443631706844690697, 0.43586652150845899941, 0.}};

// Explicit stage coefficients - Converted to U from Ascher et al. 1997 Pg 9
const double TimestepSchemeARK3B::m_dExpCf[5][5] = {
  {0., 0., 0., 0., 0.},
  {0.43586652150845899941, 0., 0., 0., 0.},
  {0.32127888602862775492, 0.39665437472560174482, 0., 0., 0.},
  {-0.1058582960718796471, 0.55292914803593982356, 0.55292914803593982356, 0., 0.},
  {0., 1.2084966491760100703, -0.6443631706844690697, 0.43586652150845899941, 0.}};
/*
// Modified MARS explicit
const double TimestepSchemeARK3B::m_dExpCf[5][5] = {
  {0., 0., 0., 0., 0.},
  {0.43586652150845899941, 0., 0., 0., 0.},
  {0.53539654030735431972, 0.18253672044687517995, 0., 0., 0.},
  {0.6304125581528726616, -0.8319339010630845273, 1.2015213429102118657, 0., 0.},
  {0., 1.2084966491760100703, -0.6443631706844690697, 0.43586652150845899941, 0.}};
/*
const double TimestepSchemeARK3B::m_dImpCf[4][4] = {
  {0.4358665215, 0., 0., 0.},
  {0.2820667392, 0.4358665215, 0., 0.},
  {1.2084966492, -0.644363171, 0.4358665215, 0.},
  {1.2084966492, -0.644363171, 0.4358665215, 0.}};

// Explicit stage coefficients - Converted to U from Ascher et al. 1997 Pg 9
const double TimestepSchemeARK3B::m_dExpCf[5][5] = {
  {0., 0., 0., 0., 0.},
  {0.4358665215, 0., 0., 0., 0.},
  {0.3212788861, 0.3966543747, 0., 0., 0.},
  {-0.105858296, 0.5529291479, 0.5529291479, 0., 0.},
  {0., 1.208496649, -0.644363171, 0.4358665215, 0.}};
*/
TimestepSchemeARK3B::TimestepSchemeARK3B(
	Model & model
) :
	TimestepScheme(model)
{
    m_dKh1Combo.Initialize(5);
    m_dKh2Combo.Initialize(7);
    m_dKh3Combo.Initialize(9);
	m_dKh4Combo.Initialize(10);
    m_dK1Combo.Initialize(4);
    m_dK2Combo.Initialize(7);
    m_dK3Combo.Initialize(9);
    m_du2fCombo.Initialize(5);
    m_du3fCombo.Initialize(7);
    m_dufCombo.Initialize(10);
}

///////////////////////////////////////////////////////////////////////////////

void TimestepSchemeARK3B::Step(
	bool fFirstStep,
	bool fLastStep,
	const Time & time,
	double dDeltaT
) {
	// Get a copy of the grid
	Grid * pGrid = m_model.GetGrid();

	// Get a copy of the HorizontalDynamics
	HorizontalDynamics * pHorizontalDynamics = m_model.GetHorizontalDynamics();

	// Get a copy of the VerticalDynamics
	VerticalDynamics * pVerticalDynamics = m_model.GetVerticalDynamics();

    // Kh1 combination
    m_dKh1Combo[0] = -1.0 / (dDeltaT * m_dExpCf[1][0]);
    m_dKh1Combo[1] = 1.0 / (dDeltaT * m_dExpCf[1][0]);
	m_dKh1Combo[2] = 0.0;
	m_dKh1Combo[3] = 0.0;
	m_dKh1Combo[4] = 0.0;

    // K1 combination
    m_dK1Combo[0] = 0.0;
    m_dK1Combo[1] = -1.0 / (dDeltaT * m_dImpCf[0][0]);
    m_dK1Combo[2] = 1.0 / (dDeltaT * m_dImpCf[0][0]);
	m_dK1Combo[3] = 0.0;

    // u2 from previous data
    m_du2fCombo[0] = 1.0;
    m_du2fCombo[1] = 0.0;
    m_du2fCombo[2] = 0.0;
    m_du2fCombo[3] = dDeltaT * m_dImpCf[1][0];
    m_du2fCombo[4] = dDeltaT * m_dExpCf[2][0];

    // K2 combination
    m_dK2Combo[0] = -1.0 / (dDeltaT * m_dImpCf[1][1]);
    m_dK2Combo[1] = 0.0;
    m_dK2Combo[2] = 1.0 / (dDeltaT * m_dImpCf[1][1]);
    m_dK2Combo[3] = -m_dImpCf[1][0] / m_dImpCf[1][1];
    m_dK2Combo[4] = -m_dExpCf[2][0] / m_dImpCf[1][1];
    m_dK2Combo[5] = 0.0;
    m_dK2Combo[6] = -m_dExpCf[2][1] / m_dImpCf[1][1];

    // Kh2 combination
    m_dKh2Combo[0] = -1.0 / (dDeltaT * m_dExpCf[2][1]);
    m_dKh2Combo[1] = 1.0 / (dDeltaT * m_dExpCf[2][1]);
    m_dKh2Combo[2] = 0.0;
    m_dKh2Combo[3] = -m_dImpCf[1][0] / m_dExpCf[2][1];
    m_dKh2Combo[4] = -m_dExpCf[2][0] / m_dExpCf[2][1];
	m_dKh2Combo[5] = 0.0;
	m_dKh2Combo[6] = 0.0;

    // u3 from previous data
    m_du3fCombo[0] = 1.0;
    m_du3fCombo[1] = 0.0;
    m_du3fCombo[2] = 0.0;
    m_du3fCombo[3] = dDeltaT * m_dImpCf[2][0];
    m_du3fCombo[4] = dDeltaT * m_dExpCf[3][0];
    m_du3fCombo[5] = dDeltaT * m_dImpCf[2][1];
    m_du3fCombo[6] = dDeltaT * m_dExpCf[3][1];

    // K3 combination
    m_dK3Combo[0] = -1.0 / (dDeltaT * m_dImpCf[2][2]);
    m_dK3Combo[1] = 0.0;
    m_dK3Combo[2] = 1.0 / (dDeltaT * m_dImpCf[2][2]);
    m_dK3Combo[3] = -m_dImpCf[2][0] / m_dImpCf[2][2];
    m_dK3Combo[4] = -m_dExpCf[3][0] / m_dImpCf[2][2];
    m_dK3Combo[5] = -m_dImpCf[2][1] / m_dImpCf[2][2];
    m_dK3Combo[6] = -m_dExpCf[3][1] / m_dImpCf[2][2];
    m_dK3Combo[7] = 0.0;
    m_dK3Combo[8] = -m_dExpCf[3][2] / m_dImpCf[2][2];

    // Kh3 combination
    m_dKh3Combo[0] = -1.0 / (dDeltaT * m_dExpCf[3][2]);
    m_dKh3Combo[1] = 1.0 / (dDeltaT * m_dExpCf[3][2]);
    m_dKh3Combo[2] = 0.0;
    m_dKh3Combo[3] = -m_dImpCf[2][0] / m_dExpCf[3][2];
    m_dKh3Combo[4] = -m_dExpCf[3][0] / m_dExpCf[3][2];
    m_dKh3Combo[5] = -m_dImpCf[2][1] / m_dExpCf[3][2];
    m_dKh3Combo[6] = -m_dExpCf[3][1] / m_dExpCf[3][2];
	m_dKh3Combo[7] = 0.0;
	m_dKh3Combo[8] = 0.0;

	// Kh4 combination for final evaluation
	m_dKh4Combo[0] = -1.0 / (dDeltaT * m_dExpCf[4][3]);
    m_dKh4Combo[1] = 1.0 / (dDeltaT * m_dExpCf[4][3]);
    m_dKh4Combo[2] = 0.0;
	m_dKh4Combo[3] = -m_dImpCf[2][0] / m_dExpCf[4][3];
	m_dKh4Combo[4] = -m_dExpCf[3][0] / m_dExpCf[4][3];
	m_dKh4Combo[5] = -m_dImpCf[2][1] / m_dExpCf[4][3];
	m_dKh4Combo[6] = -m_dExpCf[3][1] / m_dExpCf[4][3];
	m_dKh4Combo[7] = -m_dImpCf[2][2] / m_dExpCf[4][3];
	m_dKh4Combo[8] = -m_dExpCf[3][2] / m_dExpCf[4][3];
	m_dKh4Combo[9] = 0.0;

    // SUBSTAGE 1
    // Compute the explicit step to index 1
    Time timeSub1 = time;
    double dtSub1 = m_dExpCf[1][0] * dDeltaT;
    pGrid->CopyData(0, 1, DataType_State);
    pHorizontalDynamics->StepExplicit(0, 1, timeSub1, dtSub1);
	pVerticalDynamics->StepExplicit(0, 1, timeSub1, dtSub1);
    pGrid->PostProcessSubstage(1, DataType_State);
    pGrid->PostProcessSubstage(1, DataType_Tracers);

    // Store the evaluation Kh1 to index 4
    pGrid->LinearCombineData(m_dKh1Combo, 4, DataType_State);
    pGrid->PostProcessSubstage(4, DataType_State);
    pGrid->PostProcessSubstage(4, DataType_Tracers);

    // Compute implicit step based on known data to index 2
    pGrid->CopyData(1, 2, DataType_State);
    dtSub1 = m_dImpCf[0][0] * dDeltaT;
    pVerticalDynamics->StepImplicit(1, 2, timeSub1, dtSub1);
    pGrid->PostProcessSubstage(2, DataType_State);
    pGrid->PostProcessSubstage(2, DataType_Tracers);

    // Store the evaluation K1 to index 3
    pGrid->LinearCombineData(m_dK1Combo, 3, DataType_State);
    pGrid->PostProcessSubstage(3, DataType_State);
    pGrid->PostProcessSubstage(3, DataType_Tracers);

    //std::cout << "Substage 1 done ... \n";

    // SUBSTAGE 2
    // Compute uf2
    Time timeSub2 = time;
    double dtSub2 = m_dExpCf[2][1] * dDeltaT;
    pGrid->LinearCombineData(m_du2fCombo, 1, DataType_State);
    pHorizontalDynamics->StepExplicit(2, 1, timeSub2, dtSub2);
    pVerticalDynamics->StepExplicit(2, 1, timeSub2, dtSub2);
    pGrid->PostProcessSubstage(1, DataType_State);
    pGrid->PostProcessSubstage(1, DataType_Tracers);

    // Store the evaluation Kh2 to index 6
    pGrid->LinearCombineData(m_dKh2Combo, 6, DataType_State);
    pGrid->PostProcessSubstage(6, DataType_State);
    pGrid->PostProcessSubstage(6, DataType_Tracers);

    // Compute u2 from uf2 and store it to index 2 (over u1)
    dtSub2 = m_dImpCf[1][1] * dDeltaT;
    pGrid->CopyData(1, 2, DataType_State);
    pVerticalDynamics->StepImplicit(1, 2, timeSub2, dtSub2);
    pGrid->PostProcessSubstage(2, DataType_State);
    pGrid->PostProcessSubstage(2, DataType_Tracers);

    // Store the evaluation K2 to index 5
    pGrid->LinearCombineData(m_dK2Combo, 5, DataType_State);
    pGrid->PostProcessSubstage(5, DataType_State);
    pGrid->PostProcessSubstage(5, DataType_Tracers);

    //std::cout << "Substage 2 done ... \n";

    // SUBSTAGE 3
    // Compute uf3
    Time timeSub3 = time;
    double dtSub3 = m_dExpCf[3][2] * dDeltaT;
    pGrid->LinearCombineData(m_du3fCombo, 1, DataType_State);
    pHorizontalDynamics->StepExplicit(2, 1, timeSub2, dtSub2);
    pVerticalDynamics->StepExplicit(2, 1, timeSub2, dtSub2);
    pGrid->PostProcessSubstage(1, DataType_State);
    pGrid->PostProcessSubstage(1, DataType_Tracers);

    // Store the evaluation Kh3 to index 8
    pGrid->LinearCombineData(m_dKh3Combo, 8, DataType_State);
    pGrid->PostProcessSubstage(8, DataType_State);
    pGrid->PostProcessSubstage(8, DataType_Tracers);

    // Compute u3 from uf3 and store it to index 2 (over u2)
    dtSub3 = m_dImpCf[2][2] * dDeltaT;
    pGrid->CopyData(1, 2, DataType_State);
    pVerticalDynamics->StepImplicit(1, 2, timeSub3, dtSub3);
    pGrid->PostProcessSubstage(2, DataType_State);
    pGrid->PostProcessSubstage(2, DataType_Tracers);

    // Store the evaluation K3 to index 7
    pGrid->LinearCombineData(m_dK3Combo, 7, DataType_State);
    pGrid->PostProcessSubstage(7, DataType_State);
    pGrid->PostProcessSubstage(7, DataType_Tracers);
/*
	// full evaluation of u3 to compute the last explicit evaluation
    m_dufCombo[0] = 1.0;
    m_dufCombo[1] = 0.0;
    m_dufCombo[2] = 0.0;
    m_dufCombo[3] = dDeltaT * m_dImpCf[2][0];
    m_dufCombo[4] = dDeltaT * m_dExpCf[3][0];
    m_dufCombo[5] = dDeltaT * m_dImpCf[2][1];
    m_dufCombo[6] = dDeltaT * m_dExpCf[3][1];
    m_dufCombo[7] = dDeltaT * m_dImpCf[2][2];
    m_dufCombo[8] = dDeltaT * m_dExpCf[3][2];
	m_dufCombo[9] = 0.0;
*/
	// Compute the last explicit function evaluation
	dtSub3 = m_dExpCf[4][3] * dDeltaT;
	pGrid->CopyData(2, 1, DataType_State);
	pHorizontalDynamics->StepExplicit(2, 1, timeSub3, dtSub3);
    pVerticalDynamics->StepExplicit(2, 1, timeSub3, dtSub3);
    pGrid->PostProcessSubstage(1, DataType_State);
    pGrid->PostProcessSubstage(1, DataType_Tracers);

	// Store the evaluation Kh4 to index 9
    pGrid->LinearCombineData(m_dKh4Combo, 9, DataType_State);
    pGrid->PostProcessSubstage(9, DataType_State);
    pGrid->PostProcessSubstage(9, DataType_Tracers);

	// final stage evaluation with IMEX weights (last row in tableaux)
    m_dufCombo[0] = 1.0;
    m_dufCombo[1] = 0.0;
    m_dufCombo[2] = 0.0;
    m_dufCombo[3] = dDeltaT * m_dImpCf[3][0];
    m_dufCombo[4] = dDeltaT * m_dExpCf[4][0];
    m_dufCombo[5] = dDeltaT * m_dImpCf[3][1];
    m_dufCombo[6] = dDeltaT * m_dExpCf[4][1];
    m_dufCombo[7] = dDeltaT * m_dImpCf[3][2];
    m_dufCombo[8] = dDeltaT * m_dExpCf[4][2];
	m_dufCombo[9] = dDeltaT * m_dExpCf[4][3];

    // FINAL STAGE
    // Compute u+1 with the weights
	pGrid->LinearCombineData(m_dufCombo, 1, DataType_State);
	pGrid->PostProcessSubstage(1, DataType_State);
    pGrid->PostProcessSubstage(1, DataType_Tracers);

    // Apply hyperdiffusion at the end of the explicit substep (ask Paul)
	pGrid->CopyData(1, 2, DataType_State);
	pHorizontalDynamics->StepAfterSubCycle(2, 1, 3, time, dDeltaT);
	pGrid->CopyData(1, 0, DataType_State);
}

///////////////////////////////////////////////////////////////////////////////

