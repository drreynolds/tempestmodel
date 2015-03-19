///
///	\file    TimestepSchemeNASA_ARK3.cpp
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

#include "TimestepSchemeNASA_ARK3.h"
#include "Model.h"
#include "Grid.h"
#include "HorizontalDynamics.h"
#include "VerticalDynamics.h"

///////////////////////////////////////////////////////////////////////////////
// THESE COEFFICIENTS ARE COMPUTED FROM THE ORIGINAL TABLEAUX
// IMPLEMENTS Carpenter and Kennedy 2003 ARK3(2)4L[2]SA
// Time step coefficients
const double TimestepSchemeNASA_ARK3::m_dTimeCf[4] = {1./6., 1./6., 1./6., 1.};
// Implicit stage coefficients

const double TimestepSchemeNASA_ARK3::m_dImpCf[5][5] = {
  {0.1, 0., 0., 0., 0.},
  {0.435866521508459, 0.435866521508459, 0., 0., 0.},
  {0.257648246066427, -0.093514767574886, 0.435866521508459, 0., 0.},
  {0.187641024346724, -0.595297473576955, 0.971789927721772, 0.435866521508459, 0.},
  {0.187641024346724, -0.595297473576955, 0.971789927721772, 0.435866521508459, 0.}};

// Explicit stage coefficients
const double TimestepSchemeNASA_ARK3::m_dExpCf[5][5] = {
  {0., 0., 0., 0., 0.},
  {0.871733043016918, 0., 0., 0., 0.},
  {0.527589011976300, 0.072410988023700, 0., 0., 0.},
  {0.399096007676070, -0.437557654613519, 1.038461646937449, 0., 0.},
  {0.187641024346724, -0.595297473576955, 0.971789927721772, 0.435866521508459, 0.}};

/*
// Implicit stage coefficients MODIFIED BY BOSCARINO
const double TimestepSchemeNASA_ARK3::m_dImpCf[5][5] = {
  {0.1, 0., 0., 0., 0.},
  {0.435866521508459, 0.435866521508459, 0., 0., 0.},
  {-4.300026621769234, 2.265413383463721, 0.435866521508459, 0., 0.},
  {0.604248324588005, 0.0, -0.040114846096464, 0.435866521508459, 0.},
  {0.604248324588005, 0.0, -0.040114846096464, 0.435866521508459, 0.}};

// Explicit stage coefficients MODIFIED BY BOSCARINO
const double TimestepSchemeNASA_ARK3::m_dExpCf[5][5] = {
  {0., 0., 0., 0., 0.},
  {0.871733043016918, 0., 0., 0., 0.},
  {-3.064786741862248, 1.466040025065194, 0., 0., 0.},
  {0.214445607621330, 0.710753649652699, 0.074800742725971, 0., 0.},
  {0.604248324588005, 0.0, -0.040114846096464, 0.435866521508459, 0.}};
*/
TimestepSchemeNASA_ARK3::TimestepSchemeNASA_ARK3(
	Model & model
) :
	TimestepScheme(model)
{
    m_dKh1Combo.Initialize(5);
    m_dKh2Combo.Initialize(7);
    m_dKh3Combo.Initialize(9);
	m_dK0Combo.Initialize(4);
    m_dK1Combo.Initialize(6);
    m_dK2Combo.Initialize(8);
    m_dK3Combo.Initialize(10);
	m_du1fCombo.Initialize(4);
    m_du2fCombo.Initialize(6);
    m_du3fCombo.Initialize(8);
    m_du4fCombo.Initialize(10);
}

///////////////////////////////////////////////////////////////////////////////

void TimestepSchemeNASA_ARK3::Step(
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

	// K0 combination
	m_dK0Combo[0] = -1.0 / (dDeltaT * m_dImpCf[0][0]);;
    m_dK0Combo[1] = 1.0 / (dDeltaT * m_dImpCf[0][0]);;
    m_dK0Combo[2] = 0.0;
	m_dK0Combo[3] = 0.0;

	// u1 implicit evaluation combination
    m_du1fCombo[0] = 1.0;
    m_du1fCombo[1] = 0.0;
    m_du1fCombo[2] = 0.0;
    m_du1fCombo[3] = dDeltaT * m_dImpCf[1][0];

    // Kh1 combination
    m_dKh1Combo[0] = -1.0 / (dDeltaT * m_dExpCf[1][0]);
    m_dKh1Combo[1] = 1.0 / (dDeltaT * m_dExpCf[1][0]);
	m_dKh1Combo[2] = 0.0;
	m_dKh1Combo[3] = -m_dImpCf[1][0] / m_dExpCf[1][0];
	m_dKh1Combo[4] = 0.0;

    // K1 combination
    m_dK1Combo[0] = -1.0 / (dDeltaT * m_dImpCf[1][1]);
    m_dK1Combo[1] = 0.0;
    m_dK1Combo[2] = 1.0 / (dDeltaT * m_dImpCf[1][1]);
	m_dK1Combo[3] = -m_dImpCf[1][0] / m_dImpCf[1][1];
	m_dK1Combo[4] = -m_dExpCf[1][0] / m_dImpCf[1][1];
	m_dK1Combo[5] = 0.0;

    // u2 implicit evaluation combination
    m_du2fCombo[0] = 1.0;
    m_du2fCombo[1] = 0.0;
    m_du2fCombo[2] = 0.0;
    m_du2fCombo[3] = dDeltaT * m_dImpCf[2][0];
    m_du2fCombo[4] = dDeltaT * m_dExpCf[2][0];
	m_du2fCombo[5] = dDeltaT * m_dImpCf[2][1];

	// Kh2 combination
    m_dKh2Combo[0] = -1.0 / (dDeltaT * m_dExpCf[2][1]);
    m_dKh2Combo[1] = 1.0 / (dDeltaT * m_dExpCf[2][1]);
    m_dKh2Combo[2] = 0.0;
    m_dKh2Combo[3] = -m_dImpCf[2][0] / m_dExpCf[2][1];
    m_dKh2Combo[4] = -m_dExpCf[2][0] / m_dExpCf[2][1];
	m_dKh2Combo[5] = -m_dImpCf[2][1] / m_dExpCf[2][1];
	m_dKh2Combo[6] = 0.0;

    // K2 combination
    m_dK2Combo[0] = -1.0 / (dDeltaT * m_dImpCf[2][2]);
    m_dK2Combo[1] = 0.0;
    m_dK2Combo[2] = 1.0 / (dDeltaT * m_dImpCf[2][2]);
    m_dK2Combo[3] = -m_dImpCf[2][0] / m_dImpCf[2][2];
    m_dK2Combo[4] = -m_dExpCf[2][0] / m_dImpCf[2][2];
    m_dK2Combo[5] = -m_dImpCf[2][1] / m_dImpCf[2][2];
    m_dK2Combo[6] = -m_dExpCf[2][1] / m_dImpCf[2][2];
	m_dK2Combo[7] = 0.0;

    // u3 implicit evaluation combination
    m_du3fCombo[0] = 1.0;
    m_du3fCombo[1] = 0.0;
    m_du3fCombo[2] = 0.0;
    m_du3fCombo[3] = dDeltaT * m_dImpCf[3][0];
    m_du3fCombo[4] = dDeltaT * m_dExpCf[3][0];
    m_du3fCombo[5] = dDeltaT * m_dImpCf[3][1];
    m_du3fCombo[6] = dDeltaT * m_dExpCf[3][1];
	m_du3fCombo[7] = dDeltaT * m_dImpCf[3][2];

	// Kh3 combination
    m_dKh3Combo[0] = -1.0 / (dDeltaT * m_dExpCf[3][2]);
    m_dKh3Combo[1] = 1.0 / (dDeltaT * m_dExpCf[3][2]);
    m_dKh3Combo[2] = 0.0;
    m_dKh3Combo[3] = -m_dImpCf[3][0] / m_dExpCf[3][2];
    m_dKh3Combo[4] = -m_dExpCf[3][0] / m_dExpCf[3][2];
    m_dKh3Combo[5] = -m_dImpCf[3][1] / m_dExpCf[3][2];
    m_dKh3Combo[6] = -m_dExpCf[3][1] / m_dExpCf[3][2];
	m_dKh3Combo[7] = -m_dImpCf[3][2] / m_dExpCf[3][2];
	m_dKh3Combo[8] = 0.0;

    // K3 combination
    m_dK3Combo[0] = -1.0 / (dDeltaT * m_dImpCf[3][3]);
    m_dK3Combo[1] = 0.0;
    m_dK3Combo[2] = 1.0 / (dDeltaT * m_dImpCf[3][3]);
    m_dK3Combo[3] = -m_dImpCf[3][0] / m_dImpCf[3][3];
    m_dK3Combo[4] = -m_dExpCf[3][0] / m_dImpCf[3][3];
    m_dK3Combo[5] = -m_dImpCf[3][1] / m_dImpCf[3][3];
    m_dK3Combo[6] = -m_dExpCf[3][1] / m_dImpCf[3][3];
    m_dK3Combo[7] = -m_dImpCf[3][2] / m_dImpCf[3][3];
    m_dK3Combo[8] = -m_dExpCf[3][2] / m_dImpCf[3][3];
	m_dK3Combo[9] = 0.0;

    // uf4 explicit evaluation combination
    m_du4fCombo[0] = 1.0;
    m_du4fCombo[1] = 0.0;
    m_du4fCombo[2] = 0.0;
    m_du4fCombo[3] = dDeltaT * m_dImpCf[4][0];
    m_du4fCombo[4] = dDeltaT * m_dExpCf[4][0];
    m_du4fCombo[5] = dDeltaT * m_dImpCf[4][1];
    m_du4fCombo[6] = dDeltaT * m_dExpCf[4][1];
    m_du4fCombo[7] = dDeltaT * m_dImpCf[4][2];
    m_du4fCombo[8] = dDeltaT * m_dExpCf[4][2];
	m_du4fCombo[9] = dDeltaT * m_dImpCf[4][3];

	// PRE-STAGE 1 - Implicit solve to start
	Time timeSub0 = time;
    double dtSub0 = m_dImpCf[0][0] * dDeltaT;
    pGrid->CopyData(0, 1, DataType_State);
    pHorizontalDynamics->StepImplicit(0, 1, timeSub0, dtSub0);
	pVerticalDynamics->StepImplicit(0, 1, timeSub0, dtSub0);
    pGrid->PostProcessSubstage(1, DataType_State);
    pGrid->PostProcessSubstage(1, DataType_Tracers);

    // Store the evaluation K1 to index 3
    pGrid->LinearCombineData(m_dK0Combo, 3, DataType_State);

    // SUBSTAGE 1
    // Compute the explicit step to index 1
    Time timeSub1 = time;
    double dtSub1 = m_dExpCf[1][0] * dDeltaT;
    pGrid->LinearCombineData(m_du1fCombo, 1, DataType_State);
    pHorizontalDynamics->StepExplicit(0, 1, timeSub1, dtSub1);
	pVerticalDynamics->StepExplicit(0, 1, timeSub1, dtSub1);
    pGrid->PostProcessSubstage(1, DataType_State);
    pGrid->PostProcessSubstage(1, DataType_Tracers);

    // Store the evaluation Kh1 to index 4
    pGrid->LinearCombineData(m_dKh1Combo, 4, DataType_State);

    // Compute implicit step based on known data to index 2
    pGrid->CopyData(1, 2, DataType_State);
    dtSub1 = m_dImpCf[1][1] * dDeltaT;
    pVerticalDynamics->StepImplicit(1, 2, timeSub1, dtSub1);
    pGrid->PostProcessSubstage(2, DataType_State);
    pGrid->PostProcessSubstage(2, DataType_Tracers);

    // Store the evaluation K1 to index 3
    pGrid->LinearCombineData(m_dK1Combo, 5, DataType_State);

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

    // Compute u2 from uf2 and store it to index 2 (over u1)
    dtSub2 = m_dImpCf[2][2] * dDeltaT;
    pGrid->CopyData(1, 2, DataType_State);
    pVerticalDynamics->StepImplicit(1, 2, timeSub2, dtSub2);
    pGrid->PostProcessSubstage(2, DataType_State);
    pGrid->PostProcessSubstage(2, DataType_Tracers);

    // Store the evaluation K2 to index 7
    pGrid->LinearCombineData(m_dK2Combo, 7, DataType_State);

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

    // Compute u3 from uf3 and store it to index 2 (over u2)
    dtSub3 = m_dImpCf[3][3] * dDeltaT;
    pGrid->CopyData(1, 2, DataType_State);
    pVerticalDynamics->StepImplicit(1, 2, timeSub3, dtSub3);
    pGrid->PostProcessSubstage(2, DataType_State);
    pGrid->PostProcessSubstage(2, DataType_Tracers);

    // Store the evaluation K3 to index 9
    pGrid->LinearCombineData(m_dK3Combo, 9, DataType_State);

	//std::cout << "Substage 3 done ... \n";

    // SUBSTAGE 4 - MODIFIED DUE TO EXPLICIT ZERO EVALUATION
    // Compute uf4 from u3 and store it to index 2
    Time timeSub4 = time;
    double dtSub4 = m_dExpCf[4][3] * dDeltaT;
    pGrid->LinearCombineData(m_du4fCombo, 1, DataType_State);
    pHorizontalDynamics->StepExplicit(2, 1, timeSub4, dtSub4);
    pVerticalDynamics->StepExplicit(2, 1, timeSub4, dtSub4);
    pGrid->PostProcessSubstage(1, DataType_State);
    pGrid->PostProcessSubstage(1, DataType_Tracers);

    // Apply hyperdiffusion at the end of the explicit substep (ask Paul)
	//pGrid->CopyData(2, 1, DataType_State);
	//pGrid->CopyData(1, 2, DataType_State);
	pHorizontalDynamics->StepAfterSubCycle(0, 1, 3, time, dDeltaT);
	pGrid->CopyData(1, 0, DataType_State);
}

///////////////////////////////////////////////////////////////////////////////
