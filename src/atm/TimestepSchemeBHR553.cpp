///
///	\file    TimestepSchemeBHR553.cpp
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

#include "TimestepSchemeBHR553.h"
#include "Model.h"
#include "Grid.h"
#include "HorizontalDynamics.h"
#include "VerticalDynamics.h"

///////////////////////////////////////////////////////////////////////////////
// THESE COEFFICIENTS ARE COMPUTED FROM THE ORIGINAL TABLEAUX
// IMPLEMENTS BHR(5,5,3) FROM BOSCARINO AND RUSSO 2009 SSP 3RD ORDER METHOD
// Tableau coefficients
/*
const double TimestepSchemeBHR553::m_dgamma = 0.43586652150845899941;
const double TimestepSchemeBHR553::m_c2 = 2.0 * m_dgamma;
const double TimestepSchemeBHR553::m_c3 = 902905985686.0 / 1035759735069.0;
const double TimestepSchemeBHR553::m_c4 = 2684624.0 / 1147171.0;
const double TimestepSchemeBHR553::m_b1 = 487698502336740678603511.0 / 
                                        1181159636928185920260208.0;
const double TimestepSchemeBHR553::m_b2 = 0.0;
const double TimestepSchemeBHR553::m_b3 = 302987763081184622639300143137943089.0 /
				                        1535359944203293318639180129368156500.0;
const double TimestepSchemeBHR553::m_b4 = -105235928335100616072938218863.0 /
                                        2282554452064661756575727198000.0;
const double TimestepSchemeBHR553::m_at31 = m_dgamma;
const double TimestepSchemeBHR553::m_at32 = m_at31;
const double TimestepSchemeBHR553::m_at41 = -475883375220285986033264.0 /
                                          594112726933437845704163.0;
const double TimestepSchemeBHR553::m_at42 = 0.0;
const double TimestepSchemeBHR553::m_at43 = 1866233449822026827708736.0 / 
                                          594112726933437845704163.0;
const double TimestepSchemeBHR553::m_at51 = 62828845818073169585635881686091391737610308247.0 /
                                          176112910684412105319781630311686343715753056000.0;
const double TimestepSchemeBHR553::m_at52 = -m_b3;
const double TimestepSchemeBHR553::m_at53 = 262315887293043739337088563996093207.0 /
                                          297427554730376353252081786906492000.0;
const double TimestepSchemeBHR553::m_at54 = -987618231894176581438124717087.0 /
                                          23877337660202969319526901856000.0;
const double TimestepSchemeBHR553::m_a31 = m_at31;
const double TimestepSchemeBHR553::m_a32 = -31733082319927313.0 /
                                         455705377221960889379854647102.0;
const double TimestepSchemeBHR553::m_a41 = -3012378541084922027361996761794919360516301377809610.0 /
                                         45123394056585269977907753045030512597955897345819349.0;
const double TimestepSchemeBHR553::m_a42 = -62865589297807153294268.0 /
                                         102559673441610672305587327019095047.0;
const double TimestepSchemeBHR553::m_a43 = 418769796920855299603146267001414900945214277000.0 /
                                         212454360385257708555954598099874818603217167139.0;
*/
// DIFFERENT VALUE OF GAMMA (Boscarino 2009)
const double TimestepSchemeBHR553::m_dgamma = 0.57281606248208;
const double TimestepSchemeBHR553::m_c2 = 2.0 * m_dgamma;
const double TimestepSchemeBHR553::m_c3 = 12015769930846.0 / 24446477850549.0;
const double TimestepSchemeBHR553::m_c4 = 3532944.0 / 5360597.0;
const double TimestepSchemeBHR553::m_b1 = -2032971420760927701493589.0 / 
                                          38017147656515384190997416.0;
const double TimestepSchemeBHR553::m_b2 = 0.0;
const double TimestepSchemeBHR553::m_b3 = 2197602776651676983265261109643897073447.0 /
				                          945067123279139583549933947379097184164.0;
const double TimestepSchemeBHR553::m_b4 = -128147215194260398070666826235339.0 /
                                          69468482710687503388562952626424.0;
const double TimestepSchemeBHR553::m_at31 = 473447115440655855452482357894373.0 /
                                            1226306256343706154920072735579148.0;
const double TimestepSchemeBHR553::m_at32 = 129298766034131882323069978722019.0 /
                                            1226306256343706154920072735579148.0;
const double TimestepSchemeBHR553::m_at41 = 37498105210828143724516848.0 /
                                            172642583546398006173766007.0;
const double TimestepSchemeBHR553::m_at42 = 0.0;
const double TimestepSchemeBHR553::m_at43 = 76283359742561480140804416.0 / 172642583546398006173766007.0;
const double TimestepSchemeBHR553::m_at51 = -3409975860212064612303539855622639333030782744869519.0 /
                                             5886704102363745137792385361113084313351870216475136.0;
const double TimestepSchemeBHR553::m_at52 = -237416352433826978856941795734073.0 / 
                                             554681702576878342891447163499456.0;
const double TimestepSchemeBHR553::m_at53 = 4298159710546228783638212411650783228275.0 /
                                            2165398513352098924587211488610407046208.0;
const double TimestepSchemeBHR553::m_at54 = 6101865615855760853571922289749.0 / 
                                            272863973025878249803640374568448.0;
const double TimestepSchemeBHR553::m_a31 = 259252258169672523902708425780469319755.0 / 
                                           4392887760843243968922388674191715336228.0;
const double TimestepSchemeBHR553::m_a32 = -172074174703261986564706189586177.0 / 
                                           1226306256343706154920072735579148.0;
const double TimestepSchemeBHR553::m_a41 = 1103202061574553405285863729195740268785131739395559693754.0 /                            
                                           9879457735937277070641522414590493459028264677925767305837.0;
const double TimestepSchemeBHR553::m_a42 = -103754520567058969566542556296087324094.0 / 
                                           459050363888246734833121482275319954529.0;
const double TimestepSchemeBHR553::m_a43 = 3863207083069979654596872190377240608602701071947128.0 /
                                           19258690251287609765240683320611425745736762681950551.0;

const double TimestepSchemeBHR553::m_dTimeCf[4] = {m_c2, m_c3, m_c4, 1.};
// Implicit stage coefficients
const double TimestepSchemeBHR553::m_dImpCf[6][6] = {
  {m_dgamma, 0., 0., 0., 0., 0.},
  {m_dgamma, m_dgamma, 0., 0., 0., 0.},
  {m_a31, m_a32, m_dgamma, 0., 0., 0.},
  {m_a41, m_a42, m_a43, m_dgamma, 0., 0.},
  {m_b1, m_b2, m_b3, m_b4, m_dgamma, 0.},
  {m_b1, m_b2, m_b3, m_b4, m_dgamma, 0.}};

// Explicit stage coefficients
const double TimestepSchemeBHR553::m_dExpCf[6][6] = {
  {0., 0., 0., 0., 0., 0.},
  {2.0 * m_dgamma, 0., 0., 0., 0., 0.},
  {m_at31, m_at32, 0., 0., 0., 0.},
  {m_at41, m_at42, m_at43, 0., 0., 0.},
  {m_at51, m_at52, m_at53, m_at54, 0., 0.},
  {m_b1, m_b2, m_b3, m_b4, m_dgamma, 0.}};

TimestepSchemeBHR553::TimestepSchemeBHR553(
	Model & model
) :
	TimestepScheme(model)
{
    m_dKh1Combo.Initialize(5);
    m_dKh2Combo.Initialize(7);
    m_dKh3Combo.Initialize(9);
	m_dKh4Combo.Initialize(11);
	m_dKh5Combo.Initialize(13);
	m_dK0Combo.Initialize(4);
    m_dK1Combo.Initialize(6);
    m_dK2Combo.Initialize(8);
    m_dK3Combo.Initialize(10);
	m_dK4Combo.Initialize(12);
	m_du1fCombo.Initialize(6);
    m_du2fCombo.Initialize(6);
    m_du3fCombo.Initialize(8);
    m_du4fCombo.Initialize(10);
	m_dufCombo.Initialize(12);
}

///////////////////////////////////////////////////////////////////////////////

void TimestepSchemeBHR553::Step(
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
	m_dK0Combo[0] = -1.0 / (dDeltaT * m_dImpCf[0][0]);
    m_dK0Combo[1] = 0.0;
    m_dK0Combo[2] = 1.0 / (dDeltaT * m_dImpCf[0][0]);
	m_dK0Combo[3] = 0.0;

	// u1 implicit evaluation combination
    m_du1fCombo[0] = 1.0;
    m_du1fCombo[1] = 0.0;
    m_du1fCombo[2] = 0.0;
    m_du1fCombo[3] = dDeltaT * m_dImpCf[1][0];
	m_du1fCombo[4] = 0.0;
	m_du1fCombo[5] = 0.0;

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
/*
    // Kh1 combination
    m_dKh1Combo[0] = -1.0 / (dDeltaT * m_dExpCf[1][0]);
    m_dKh1Combo[1] = 1.0 / (dDeltaT * m_dExpCf[1][0]);
	m_dKh1Combo[2] = 0.0;
	m_dKh1Combo[3] = 0.0;
	m_dKh1Combo[4] = 0.0;

    // K1 combination
    m_dK1Combo[0] = 0.0;
    m_dK1Combo[1] = -1.0 / (dDeltaT * m_dImpCf[1][1]);
    m_dK1Combo[2] = 1.0 / (dDeltaT * m_dImpCf[1][1]);
	m_dK1Combo[3] = 0.0;
	m_dK1Combo[4] = 0.0;
	m_dK1Combo[5] = 0.0;
*/
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

    // Kh4 combination
    m_dKh4Combo[0] = -1.0 / (dDeltaT * m_dExpCf[4][3]);
    m_dKh4Combo[1] = 1.0 / (dDeltaT * m_dExpCf[4][3]);
    m_dKh4Combo[2] = 0.0;
    m_dKh4Combo[3] = -m_dImpCf[4][0] / m_dExpCf[4][3];
    m_dKh4Combo[4] = -m_dExpCf[4][0] / m_dExpCf[4][3];
    m_dKh4Combo[5] = -m_dImpCf[4][1] / m_dExpCf[4][3];
    m_dKh4Combo[6] = -m_dExpCf[4][1] / m_dExpCf[4][3];
	m_dKh4Combo[7] = -m_dImpCf[4][2] / m_dExpCf[4][3];
	m_dKh4Combo[8] = -m_dExpCf[4][2] / m_dExpCf[4][3];
	m_dKh4Combo[9] = -m_dImpCf[4][3] / m_dExpCf[4][3];
	m_dKh4Combo[10] = 0.0;

	// K4 combination
    m_dK4Combo[0] = -1.0 / (dDeltaT * m_dImpCf[4][4]);
    m_dK4Combo[1] = 0.0;
    m_dK4Combo[2] = 1.0 / (dDeltaT * m_dImpCf[4][4]);
    m_dK4Combo[3] = -m_dImpCf[4][0] / m_dImpCf[4][4];
    m_dK4Combo[4] = -m_dExpCf[4][0] / m_dImpCf[4][4];
    m_dK4Combo[5] = -m_dImpCf[4][1] / m_dImpCf[4][4];
    m_dK4Combo[6] = -m_dExpCf[4][1] / m_dImpCf[4][4];
    m_dK4Combo[7] = -m_dImpCf[4][2] / m_dImpCf[4][4];
    m_dK4Combo[8] = -m_dExpCf[4][2] / m_dImpCf[4][4];
	m_dK4Combo[9] = -m_dImpCf[4][3] / m_dImpCf[4][4];
    m_dK4Combo[10] = -m_dExpCf[4][3] / m_dImpCf[4][4];
	m_dK4Combo[11] = 0.0;

    // Kh5 combination
    m_dKh5Combo[0] = -1.0 / (dDeltaT * m_dExpCf[5][4]);
    m_dKh5Combo[1] = 1.0 / (dDeltaT * m_dExpCf[5][4]);
    m_dKh5Combo[2] = 0.0;
    m_dKh5Combo[3] = -m_dImpCf[4][0] / m_dExpCf[5][4];
    m_dKh5Combo[4] = -m_dExpCf[4][0] / m_dExpCf[5][4];
    m_dKh5Combo[5] = -m_dImpCf[4][1] / m_dExpCf[5][4];
    m_dKh5Combo[6] = -m_dExpCf[4][1] / m_dExpCf[5][4];
	m_dKh5Combo[7] = -m_dImpCf[4][2] / m_dExpCf[5][4];
	m_dKh5Combo[8] = -m_dExpCf[4][2] / m_dExpCf[5][4];
	m_dKh5Combo[9] = -m_dImpCf[4][3] / m_dExpCf[5][4];
	m_dKh5Combo[10] = -m_dExpCf[4][3] / m_dExpCf[5][4];
	m_dKh5Combo[11] = -m_dImpCf[4][4] / m_dExpCf[5][4];
	m_dKh5Combo[12] = 0.0;

	// PRE-STAGE 1 - Implicit solve to start
	Time timeSub0 = time;
    double dtSub0 = m_dImpCf[0][0] * dDeltaT;
    pGrid->CopyData(0, 2, DataType_State);
    pHorizontalDynamics->StepImplicit(0, 2, timeSub0, dtSub0);
	pVerticalDynamics->StepImplicit(0, 2, timeSub0, dtSub0);
    pGrid->PostProcessSubstage(2, DataType_State);
    pGrid->PostProcessSubstage(2, DataType_Tracers);

    // Store the evaluation K1 to index 3
    pGrid->LinearCombineData(m_dK0Combo, 3, DataType_State);
	pGrid->PostProcessSubstage(3, DataType_State);
    pGrid->PostProcessSubstage(3, DataType_Tracers);

    // SUBSTAGE 1
    // Compute the explicit step to index 1
    Time timeSub1 = time;
    double dtSub1 = m_dExpCf[1][0] * dDeltaT;
	pGrid->LinearCombineData(m_du1fCombo, 1, DataType_State);
    pHorizontalDynamics->StepExplicit(2, 1, timeSub1, dtSub1);
	pVerticalDynamics->StepExplicit(2, 1, timeSub1, dtSub1);
    pGrid->PostProcessSubstage(1, DataType_State);
    pGrid->PostProcessSubstage(1, DataType_Tracers);

    // Store the evaluation Kh1 to index 4
    pGrid->LinearCombineData(m_dKh1Combo, 4, DataType_State);
	pGrid->PostProcessSubstage(4, DataType_State);
    pGrid->PostProcessSubstage(4, DataType_Tracers);

    // Compute implicit step based on known data to index 2
    dtSub1 = m_dImpCf[1][1] * dDeltaT;
	pGrid->CopyData(1, 2, DataType_State);
    pVerticalDynamics->StepImplicit(1, 2, timeSub1, dtSub1);
    pGrid->PostProcessSubstage(2, DataType_State);
    pGrid->PostProcessSubstage(2, DataType_Tracers);

    // Store the evaluation K1 to index 3
    pGrid->LinearCombineData(m_dK1Combo, 5, DataType_State);
	pGrid->PostProcessSubstage(5, DataType_State);
    pGrid->PostProcessSubstage(5, DataType_Tracers);

	pGrid->LinearCombineData(m_du1fCombo, 2, DataType_State);

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
    dtSub2 = m_dImpCf[2][2] * dDeltaT;
    pGrid->CopyData(1, 2, DataType_State);
    pVerticalDynamics->StepImplicit(1, 2, timeSub2, dtSub2);
    pGrid->PostProcessSubstage(2, DataType_State);
    pGrid->PostProcessSubstage(2, DataType_Tracers);

    // Store the evaluation K2 to index 7
    pGrid->LinearCombineData(m_dK2Combo, 7, DataType_State);
	pGrid->PostProcessSubstage(7, DataType_State);
    pGrid->PostProcessSubstage(7, DataType_Tracers);

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
    dtSub3 = m_dImpCf[3][3] * dDeltaT;
    pGrid->CopyData(1, 2, DataType_State);
    pVerticalDynamics->StepImplicit(1, 2, timeSub3, dtSub3);
    pGrid->PostProcessSubstage(2, DataType_State);
    pGrid->PostProcessSubstage(2, DataType_Tracers);

    // Store the evaluation K3 to index 9
    pGrid->LinearCombineData(m_dK3Combo, 9, DataType_State);
	pGrid->PostProcessSubstage(9, DataType_State);
    pGrid->PostProcessSubstage(9, DataType_Tracers);

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

	// Store the evaluation Kh4 to index 10
    pGrid->LinearCombineData(m_dKh4Combo, 10, DataType_State);
	pGrid->PostProcessSubstage(10, DataType_State);
    pGrid->PostProcessSubstage(10, DataType_Tracers);

    // Compute u4 from uf4 and store it to index 2 (over u3)
    dtSub4 = m_dImpCf[4][4] * dDeltaT;
    pGrid->CopyData(1, 2, DataType_State);
    pVerticalDynamics->StepImplicit(1, 2, timeSub4, dtSub4);
    pGrid->PostProcessSubstage(2, DataType_State);
    pGrid->PostProcessSubstage(2, DataType_Tracers);

	// Store the evaluation K4 to index 11
    pGrid->LinearCombineData(m_dK4Combo, 11, DataType_State);
	pGrid->PostProcessSubstage(11, DataType_State);
    pGrid->PostProcessSubstage(11, DataType_Tracers);
/*
	// full evaluation of u4 to compute the last explicit evaluation
    m_dufCombo[0] = 1.0;
    m_dufCombo[1] = 0.0;
    m_dufCombo[2] = 0.0;
    m_dufCombo[3] = dDeltaT * m_dImpCf[4][0];
    m_dufCombo[4] = dDeltaT * m_dExpCf[4][0];
    m_dufCombo[5] = dDeltaT * m_dImpCf[4][1];
    m_dufCombo[6] = dDeltaT * m_dExpCf[4][1];
    m_dufCombo[7] = dDeltaT * m_dImpCf[4][2];
    m_dufCombo[8] = dDeltaT * m_dExpCf[4][2];
	m_dufCombo[9] = dDeltaT * m_dImpCf[4][3];
	m_dufCombo[10] = dDeltaT * m_dExpCf[4][3];
	m_dufCombo[11] = dDeltaT * m_dImpCf[4][4];
*/
	// Compute the last explicit function evaluation
	dtSub3 = m_dExpCf[5][4] * dDeltaT;
	pGrid->CopyData(2, 1, DataType_State);
	//pGrid->LinearCombineData(m_dufCombo, 1, DataType_State);
	pHorizontalDynamics->StepExplicit(2, 1, timeSub3, dtSub3);
    pVerticalDynamics->StepExplicit(2, 1, timeSub3, dtSub3);
    pGrid->PostProcessSubstage(1, DataType_State);
    pGrid->PostProcessSubstage(1, DataType_Tracers);

	// Store the evaluation Kh5 to index 12
    pGrid->LinearCombineData(m_dKh5Combo, 12, DataType_State);
	pGrid->PostProcessSubstage(12, DataType_State);
    pGrid->PostProcessSubstage(12, DataType_Tracers);

	// final stage evaluation with IMEX weights (last row in tableaux)
    m_dufCombo[0] = 1.0;
    m_dufCombo[1] = 0.0;
    m_dufCombo[2] = 0.0;
    m_dufCombo[3] = dDeltaT * m_dImpCf[5][0];
    m_dufCombo[4] = dDeltaT * m_dExpCf[5][0];
    m_dufCombo[5] = dDeltaT * m_dImpCf[5][1];
    m_dufCombo[6] = dDeltaT * m_dExpCf[5][1];
    m_dufCombo[7] = dDeltaT * m_dImpCf[5][2];
    m_dufCombo[8] = dDeltaT * m_dExpCf[5][2];
	m_dufCombo[9] = dDeltaT * m_dImpCf[5][3];
	m_dufCombo[10] = dDeltaT * m_dExpCf[5][3];
	m_dufCombo[11] = dDeltaT * m_dImpCf[5][4];

	//std::cout << "Substage 4 done ... \n";

    // FINAL STAGE
    // Compute u+1 with the weights
	pGrid->LinearCombineData(m_dufCombo, 1, DataType_State);
	pGrid->PostProcessSubstage(1, DataType_State);
    pGrid->PostProcessSubstage(1, DataType_Tracers);

    // Apply hyperdiffusion at the end of the explicit substep (ask Paul)
	pGrid->CopyData(1, 2, DataType_State);
	pHorizontalDynamics->StepAfterSubCycle(1, 2, 3, time, dDeltaT);
	pGrid->CopyData(2, 0, DataType_State);
}

///////////////////////////////////////////////////////////////////////////////

