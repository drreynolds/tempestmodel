///////////////////////////////////////////////////////////////////////////////
///
///	\file    VerticalDynamicsFEM.cpp
///	\author  Paul Ullrich
///	\version May 20, 2013
///
///	<remarks>
///		Copyright 2000-2013 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "VerticalDynamicsFEM.h"
#include "GaussQuadrature.h"
#include "GaussLobattoQuadrature.h"
#include "FluxReconstructionFunction.h"

#include "Model.h"
#include "Grid.h"
#include "GridCSGLL.h"
#include "EquationSet.h"
#include "TimeObj.h"
#include "PolynomialInterp.h"
#include "LinearAlgebra.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Type of flux reconstruction function to use (see Huynh 2007)
///	</summary>
static const int ParamFluxReconstructionType = 2;

///////////////////////////////////////////////////////////////////////////////

VerticalDynamicsFEM::VerticalDynamicsFEM(
	Model & model,
	int nHorizontalOrder,
	int nVerticalOrder,
	int nHyperdiffusionOrder,
	bool fFullyExplicit,
	bool fExnerPressureOnLevels,
	bool fMassFluxOnLevels
) :
	VerticalDynamics(model),
	m_nHorizontalOrder(nHorizontalOrder),
	m_nVerticalOrder(nVerticalOrder),
	m_fFullyExplicit(fFullyExplicit),
	m_fExnerPressureOnLevels(fExnerPressureOnLevels),
	m_fMassFluxOnLevels(fMassFluxOnLevels),
	m_nHyperdiffusionOrder(nHyperdiffusionOrder),
	m_dHyperdiffusionCoeff(0.0)
{
	if (nHyperdiffusionOrder % 2 == 1) {
		_EXCEPTIONT("Vertical hyperdiffusion order must be even.");
	}

	if (nHyperdiffusionOrder < 0) {
		_EXCEPTIONT("Vertical hyperdiffusion order must be positive.");
	}

	if (nHyperdiffusionOrder == 0) {

	} else if (nHyperdiffusionOrder == 2) {
		m_dHyperdiffusionCoeff = 0.2;

	} else if (nHyperdiffusionOrder == 4) {
		m_dHyperdiffusionCoeff = -0.02;

	} else {
		_EXCEPTIONT("UNIMPLEMENTED: Vertical hyperdiffusion order > 4");
	}

}

///////////////////////////////////////////////////////////////////////////////

VerticalDynamicsFEM::~VerticalDynamicsFEM() {
#ifdef USE_PETSC
	SNESDestroy(&m_snes);
	VecDestroy(&m_vecX);
	VecDestroy(&m_vecR);
	MatDestroy(&m_matJ);
#endif
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::Initialize() {

	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int TIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	// Pointer to grid
	Grid * pGrid = m_model.GetGrid();

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Number of radial elements
	int nRElements = pGrid->GetRElements();

	// Number of degrees of freedom per column in rho/w/theta
	m_nColumnStateSize = 3 * nRElements;

	if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
		m_nColumnStateSize ++;
	}
	if (pGrid->GetVarLocation(TIx) == DataLocation_REdge) {
		m_nColumnStateSize ++;
	}

	// First index of each variable in the array
	m_ixTBegin = 0;

	if (pGrid->GetVarLocation(TIx) == DataLocation_REdge) {
		m_ixWBegin = m_ixTBegin + nRElements + 1;
	} else {
		m_ixWBegin = m_ixTBegin + nRElements;
	}

	if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
		m_ixRBegin = m_ixWBegin + nRElements + 1;
	} else {
		m_ixRBegin = m_ixWBegin + nRElements;
	}

	// Parent grid, containing the vertical remapping information
	const GridGLL * pGLLGrid = dynamic_cast<const GridGLL*>(pGrid);
	if (pGLLGrid == NULL) {
		_EXCEPTIONT("Invalid grid -- expected GridGLL");
	}

	m_pInterpNodeToREdge = &(pGLLGrid->GetInterpNodeToREdge());
	m_pInterpREdgeToNode = &(pGLLGrid->GetInterpREdgeToNode());

#ifdef USE_PETSC
	// Initialize the PetSc solver context
	SNESCreate(PETSC_COMM_SELF, &m_snes);

	// Create vectors
	VecCreate(PETSC_COMM_SELF, &m_vecX);
	VecSetSizes(m_vecX, PETSC_DECIDE, m_nColumnStateSize);
	VecSetFromOptions(m_vecX);
	VecDuplicate(m_vecX, &m_vecR);

	// Set tolerances
	SNESSetTolerances(
		m_snes,
		1.0e-8,
		1.0e-8,
		1.0e-8,
		1,
		50);

	// Set the function
	SNESSetFunction(
		m_snes,
		m_vecR,
		VerticalDynamicsFEM_FormFunction,
		(void*)(this));

	MatCreateSNESMF(m_snes, &m_matJ);

	SNESSetJacobian(m_snes, m_matJ, m_matJ, MatMFFDComputeJacobian, NULL);

	// Set the SNES context from options
	SNESSetFromOptions(m_snes);
#endif
#ifdef USE_GMRES
	// Initialize JFNK
	InitializeJFNK(m_nColumnStateSize, m_nColumnStateSize);
#endif
#ifdef USE_DIRECTSOLVE_APPROXJ
	// Initialize Jacobian matrix
	m_matJacobianF.Initialize(m_nColumnStateSize, m_nColumnStateSize);

	// Initialize pivot vector
	m_vecIPiv.Initialize(m_nColumnStateSize);
#endif
#ifdef USE_DIRECTSOLVE
	// Initialize Jacobian matrix
	m_matJacobianF.Initialize(m_nColumnStateSize, m_nColumnStateSize);

	// Initialize pivot vector
	m_vecIPiv.Initialize(m_nColumnStateSize);
#endif

	// Allocate column for JFNK
	m_dColumnState.Initialize(m_nColumnStateSize);

	// Allocation reference column
	m_dStateRefNode.Initialize(5, nRElements);
	m_dStateRefREdge.Initialize(5, nRElements+1);

	// Solution vector from JFNK
	m_dSoln.Initialize(m_nColumnStateSize);

	// Get points for Gaussian quadrature
	DataVector<double> dG;
	DataVector<double> dGL;
	DataVector<double> dW;
	DataVector<double> dWL;

	// Reference element [0,1] model levels
	GaussQuadrature::GetPoints(m_nVerticalOrder, 0.0, 1.0, dG, dW);

	// Reference element [0,1] model interfaces
	GaussLobattoQuadrature::GetPoints(m_nVerticalOrder+1, 0.0, 1.0, dGL, dWL);

	// State vector at levels
	m_dStateNode.Initialize(
		m_model.GetEquationSet().GetComponents(),
		nRElements);

	// State vector at interfaces
	m_dStateREdge.Initialize(
		m_model.GetEquationSet().GetComponents(),
		nRElements+1);

	// Auxiliary variables at interfaces
	int nFiniteElements = nRElements / m_nVerticalOrder;
	if (nRElements % m_nVerticalOrder != 0) {
		_EXCEPTIONT("Logic error: Vertical order must divide RElements");
	}

	m_dStateAux.Initialize(nRElements+1);
	m_dStateAuxDiff.Initialize(nRElements+1);

	m_dDiffTheta.Initialize(nRElements+1);

	m_dStateFEEdge.Initialize(nFiniteElements+1, 2);

	m_dMassFluxNode.Initialize(nRElements);
	m_dDiffMassFluxNode.Initialize(nRElements);
	m_dMassFluxREdge.Initialize(nRElements+1);
	m_dDiffMassFluxREdge.Initialize(nRElements+1);

	m_dExnerPertNode.Initialize(nRElements);
	m_dExnerRefNode.Initialize(nRElements);

	m_dDiffExnerPertNode.Initialize(nRElements);
	m_dDiffExnerRefNode.Initialize(nRElements);

	m_dExnerPertREdge.Initialize(nRElements+1);
	m_dExnerRefREdge.Initialize(nRElements+1);

	m_dDiffExnerPertREdge.Initialize(nRElements+1);
	m_dDiffExnerRefREdge.Initialize(nRElements+1);

	// Grid spacing
	double dElementDeltaXi =
		  pGrid->GetREtaInterface(m_nVerticalOrder)
		- pGrid->GetREtaInterface(0);

	// Differentiation coefficients
	m_dDiffREdgeToNode.Initialize(m_nVerticalOrder, m_nVerticalOrder+1);
	m_dDiffREdgeToREdge.Initialize(m_nVerticalOrder+1, m_nVerticalOrder+1);
	m_dDiffNodeToREdge.Initialize(m_nVerticalOrder+1, m_nVerticalOrder);
	m_dDiffNodeToNode.Initialize(m_nVerticalOrder, m_nVerticalOrder);

	// Compute differentiation coefficients
	for (int n = 0; n < m_nVerticalOrder; n++) {
		PolynomialInterp::DiffLagrangianPolynomialCoeffs(
			m_nVerticalOrder+1, dGL, m_dDiffREdgeToNode[n], dG[n]);
		PolynomialInterp::DiffLagrangianPolynomialCoeffs(
			m_nVerticalOrder, dG, m_dDiffNodeToNode[n], dG[n]);

		for (int m = 0; m <= m_nVerticalOrder; m++) {
			m_dDiffREdgeToNode[n][m] /= dElementDeltaXi;
		}
		for (int m = 0; m < m_nVerticalOrder; m++) {
			m_dDiffNodeToNode[n][m] /= dElementDeltaXi;
		}
	}
	for (int n = 0; n <= m_nVerticalOrder; n++) {
		PolynomialInterp::DiffLagrangianPolynomialCoeffs(
			m_nVerticalOrder, dG, m_dDiffNodeToREdge[n], dGL[n]);
		PolynomialInterp::DiffLagrangianPolynomialCoeffs(
			m_nVerticalOrder+1, dGL, m_dDiffREdgeToREdge[n], dGL[n]);

		for (int m = 0; m <= m_nVerticalOrder; m++) {
			m_dDiffREdgeToREdge[n][m] /= dElementDeltaXi;
		}
		for (int m = 0; m < m_nVerticalOrder; m++) {
			m_dDiffNodeToREdge[n][m] /= dElementDeltaXi;
		}
	}

	// Compute second differentiation coefficients
	m_dDiffDiffREdgeToREdge.Initialize(
		m_nVerticalOrder+1, m_nVerticalOrder+1);

	for (int n = 0; n <= m_nVerticalOrder; n++) {
	for (int m = 0; m <= m_nVerticalOrder; m++) {
		m_dDiffDiffREdgeToREdge[n][m] = 0.0;
		for (int s = 0; s <= m_nVerticalOrder; s++) {
			m_dDiffDiffREdgeToREdge[n][m] -=
				  m_dDiffREdgeToREdge[s][n]
				* m_dDiffREdgeToREdge[s][m]
				* dWL[s];
		}
		m_dDiffDiffREdgeToREdge[n][m] /= dWL[n];
	}
	}

	// Get derivatives of flux reconstruction function and scale to the
	// element [0, dElementDeltaXi]
	FluxReconstructionFunction::GetDerivatives(
		ParamFluxReconstructionType,
		m_nVerticalOrder+1, dG, m_dDiffReconsPolyNode);

	FluxReconstructionFunction::GetDerivatives(
		ParamFluxReconstructionType,
		m_nVerticalOrder+1, dGL, m_dDiffReconsPolyREdge);

	for (int n = 0; n < m_dDiffReconsPolyNode.GetRows(); n++) {
		m_dDiffReconsPolyNode[n] /= dElementDeltaXi;
	}
	for (int n = 0; n < m_dDiffReconsPolyREdge.GetRows(); n++) {
		m_dDiffReconsPolyREdge[n] /= dElementDeltaXi;
	}

	// Compute amalgamated differentiation coefficients
	m_dDiffNodeToREdgeAmal.Initialize(
		m_nVerticalOrder+1, 3*m_nVerticalOrder);

	m_dDiffNodeToREdgeLeft.Initialize(
		m_nVerticalOrder+1, 2*m_nVerticalOrder);

	m_dDiffNodeToREdgeRight.Initialize(
		m_nVerticalOrder+1, 2*m_nVerticalOrder);

	for (int n = 0; n <= m_nVerticalOrder; n++) {
	for (int m = 0; m < m_nVerticalOrder; m++) {
		m_dDiffNodeToREdgeAmal[n][m_nVerticalOrder + m] =
			m_dDiffNodeToREdge[n][m];
		m_dDiffNodeToREdgeLeft[n][m] =
			m_dDiffNodeToREdge[n][m];
		m_dDiffNodeToREdgeRight[n][m_nVerticalOrder + m] =
			m_dDiffNodeToREdge[n][m];
	}
	}

	// Overlay differentiation stencils
	for (int m = 0; m < 2 * m_nVerticalOrder; m++) {
		m_dDiffNodeToREdgeAmal[0][m] = 0.5 * (
			  m_dDiffNodeToREdgeAmal[0][m]
			+ m_dDiffNodeToREdgeAmal[m_nVerticalOrder][m_nVerticalOrder+m]);

		m_dDiffNodeToREdgeAmal[m_nVerticalOrder][m_nVerticalOrder+m] =
			m_dDiffNodeToREdgeAmal[0][m];
	}

	// Contributions due to interface values
	for (int n = 0; n <= m_nVerticalOrder; n++) {
	for (int m = 0; m < m_nVerticalOrder; m++) {
		// Contribution from element on the right
		m_dDiffNodeToREdgeAmal[n][m_nVerticalOrder + m] -=
			0.5 * (*m_pInterpNodeToREdge)[m_nVerticalOrder][m]
			* m_dDiffReconsPolyREdge[n];

		m_dDiffNodeToREdgeAmal[n][2 * m_nVerticalOrder + m] +=
			0.5 * (*m_pInterpNodeToREdge)[0][m]
			* m_dDiffReconsPolyREdge[n];

		m_dDiffNodeToREdgeLeft[n][m] -=
			0.5 * (*m_pInterpNodeToREdge)[m_nVerticalOrder][m] * (
				m_dDiffReconsPolyREdge[n]
				+ m_dDiffReconsPolyREdge[m_nVerticalOrder - n]);

		m_dDiffNodeToREdgeLeft[n][m_nVerticalOrder + m] +=
			0.5 * (*m_pInterpNodeToREdge)[0][m] * (
				m_dDiffReconsPolyREdge[n]
				+ m_dDiffReconsPolyREdge[m_nVerticalOrder - n]);

		// Contribution from element on the left
		m_dDiffNodeToREdgeAmal[n][m] +=
			- 0.5 * (*m_pInterpNodeToREdge)[m_nVerticalOrder][m]
			* m_dDiffReconsPolyREdge[m_nVerticalOrder - n];

		m_dDiffNodeToREdgeAmal[n][m_nVerticalOrder + m] -=
			- 0.5* (*m_pInterpNodeToREdge)[0][m]
			* m_dDiffReconsPolyREdge[m_nVerticalOrder - n];

		m_dDiffNodeToREdgeRight[n][m] +=
			- 0.5 * (*m_pInterpNodeToREdge)[m_nVerticalOrder][m] * (
				m_dDiffReconsPolyREdge[m_nVerticalOrder - n]
				+ m_dDiffReconsPolyREdge[n]);

		m_dDiffNodeToREdgeRight[n][m_nVerticalOrder + m] -=
			- 0.5 * (*m_pInterpNodeToREdge)[0][m] * (
				m_dDiffReconsPolyREdge[m_nVerticalOrder - n]
				+ m_dDiffReconsPolyREdge[n]);
	}
	}

	// Calculate Exner pressure reference profile
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatch * pPatch = pGrid->GetActivePatch(n);

		const PatchBox & box = pPatch->GetPatchBox();

		// Data
		const GridData4D & dataRefNode =
			pPatch->GetReferenceState(DataLocation_Node);

		const GridData4D & dataRefREdge =
			pPatch->GetReferenceState(DataLocation_REdge);

		GridData3D & dataExnerNode =
			pPatch->GetVerticalDynamicsAuxData(0, DataLocation_Node);

		GridData3D & dataDiffExnerNode =
			pPatch->GetVerticalDynamicsAuxData(1, DataLocation_Node);

		GridData3D & dataExnerREdge =
			pPatch->GetVerticalDynamicsAuxData(0, DataLocation_REdge);

		GridData3D & dataDiffExnerREdge =
			pPatch->GetVerticalDynamicsAuxData(1, DataLocation_REdge);

		// Loop over all nodes
		for (int i = box.GetAInteriorBegin(); i < box.GetAInteriorEnd(); i++) {
		for (int j = box.GetBInteriorBegin(); j < box.GetBInteriorEnd(); j++) {

			// Initialize reference Exner pressure at model levels
			for (int k = 0; k < nRElements; k++) {
				dataExnerNode[k][i][j] =
					phys.GetCp() * exp(phys.GetR() / phys.GetCv()
						* log(phys.GetR() / phys.GetP0()
							* dataRefNode[TIx][k][i][j]
							* dataRefNode[RIx][k][i][j]));

				m_dExnerRefNode[k] = dataExnerNode[k][i][j];
			}

			// Initialize reference Exner pressure at model interfaces
			for (int k = 0; k <= nRElements; k++) {
				dataExnerREdge[k][i][j] =
					phys.GetCp() * exp(phys.GetR() / phys.GetCv()
						* log(phys.GetR() / phys.GetP0()
							* dataRefREdge[TIx][k][i][j]
							* dataRefREdge[RIx][k][i][j]));

				m_dExnerRefREdge[k] = dataExnerREdge[k][i][j];
			}

			// Differentiate the reference Exner pressure
			if ((pGrid->GetVarLocation(RIx) == DataLocation_Node) &&
				(pGrid->GetVarLocation(WIx) == DataLocation_Node)
			) {
				DifferentiateNodeToNode(
					m_dExnerRefNode,
					m_dDiffExnerRefNode);

				DifferentiateREdgeToREdge(
					m_dExnerRefREdge,
					m_dDiffExnerRefREdge);

			} else if (
				(pGrid->GetVarLocation(RIx) == DataLocation_Node) &&
				(pGrid->GetVarLocation(WIx) == DataLocation_REdge)
			) {
				if (m_fExnerPressureOnLevels) {
					DifferentiateNodeToNode(
						m_dExnerRefNode,
						m_dDiffExnerRefNode);

					DifferentiateNodeToREdge(
						m_dExnerRefNode,
						m_dDiffExnerRefREdge);

				} else {
					DifferentiateREdgeToNode(
						m_dExnerRefREdge,
						m_dDiffExnerRefNode);

					DifferentiateREdgeToREdge(
						m_dExnerRefREdge,
						m_dDiffExnerRefREdge);
				}

			} else {
				_EXCEPTIONT("Invalid variable staggering"
					" (possibly UNIMPLEMENTED)");
			}

			// Store derivatives at this point
			for (int k = 0; k < nRElements; k++) {
				dataDiffExnerNode[k][i][j] = m_dDiffExnerRefNode[k];
			}

			for (int k = 0; k <= nRElements; k++) {
				dataDiffExnerREdge[k][i][j] = m_dDiffExnerRefREdge[k];
			}
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::InterpolateNodeToREdge(
	const double * dDataNode,
	const double * dDataRefNode,
	double * dDataREdge,
	const double * dDataRefREdge,
	bool fZeroBoundaries
) {
	// Number of radial elements
	int nRElements = m_model.GetGrid()->GetRElements();

	int nFiniteElements = nRElements / m_nVerticalOrder;

	// Zero the memory
	memset(dDataREdge, 0, (nRElements+1) * sizeof(double));

	// Loop over all nodes and apply the value to all interfaces
	for (int k = 0; k < nRElements; k++) {
		int a = k / m_nVerticalOrder;
		int m = k % m_nVerticalOrder;
		int lBegin = a * m_nVerticalOrder;

		// Apply node value to interface
		for (int l = 0; l <= m_nVerticalOrder; l++) {
			dDataREdge[lBegin + l] +=
				(*m_pInterpNodeToREdge)[l][m]
					* (dDataNode[k] - dDataRefNode[k]);
		}
	}

	// Scale interior element edges
	for (int a = 1; a < nFiniteElements; a++) {
		dDataREdge[a * m_nVerticalOrder] *= 0.5;
	}
	for (int k = 0; k <= nRElements; k++) {
		dDataREdge[k] += dDataRefREdge[k];
	}

	// Override boundary values if zero
	if (fZeroBoundaries) {
		dDataREdge[0] = 0.0;
		dDataREdge[nRElements] = 0.0;
	}
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::InterpolateREdgeToNode(
	const double * dDataREdge,
	const double * dDataRefREdge,
	double * dDataNode,
	const double * dDataRefNode
) {
	// Number of radial elements
	int nRElements = m_model.GetGrid()->GetRElements();

	// Zero the memory
	memset(dDataNode, 0, nRElements * sizeof(double));

	// Loop over all nodes
	for (int k = 0; k < nRElements; k++) {
		int a = k / m_nVerticalOrder;
		int m = k % m_nVerticalOrder;
		int lBegin = a * m_nVerticalOrder;

		// Apply interface values to nodes
		for (int l = 0; l <= m_nVerticalOrder; l++) {
			dDataNode[k] +=
				(*m_pInterpREdgeToNode)[m][l] *
					(dDataREdge[lBegin + l] - dDataRefREdge[lBegin + l]);
		}
	}

	// Restore the reference state
	for (int k = 0; k < nRElements; k++) {
		dDataNode[k] += dDataRefNode[k];
	}
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::InterpolateNodeToFEEdges(
	const double * dDataNode,
	bool fZeroBoundaries
) {
	const int Left = 0;
	const int Right = 1;

	// Number of radial elements
	int nRElements = m_model.GetGrid()->GetRElements();

	int nFiniteElements = nRElements / m_nVerticalOrder;

	// Interpolate to interfaces (left and right)
	m_dStateFEEdge.Zero();

	for (int k = 0; k < nRElements; k++) {
		int a = k / m_nVerticalOrder;
		int m = k % m_nVerticalOrder;

		// Apply node value to left side of interface
		m_dStateFEEdge[a+1][Left] +=
			(*m_pInterpNodeToREdge)[m_nVerticalOrder][m] * dDataNode[k];

		// Apply node value to right side of interface
		m_dStateFEEdge[a][Right] +=
			(*m_pInterpNodeToREdge)[0][m] * dDataNode[k];
	}

	// Calculate average interpolant
	for (int a = 1; a < nFiniteElements; a++) {
		double dAvg = 0.5 * (m_dStateFEEdge[a][0] + m_dStateFEEdge[a][1]);

		m_dStateFEEdge[a][Left] -= dAvg;
		m_dStateFEEdge[a][Right] -= dAvg;
	}

#pragma message "Understand why this works for free boundaries"
	// Ignore contributions due to upper and lower boundary
	// NOTE: This needs to be changed for W to enforce boundary conditions
	if (fZeroBoundaries) {

	} else if (nFiniteElements == 1) {
		m_dStateFEEdge[0][Right] = 0.0;
		m_dStateFEEdge[nFiniteElements][Left] = 0.0;

	} else {
		m_dStateFEEdge[0][Right] =
			- m_dStateFEEdge[1][Left];
		m_dStateFEEdge[nFiniteElements][Left] =
			- m_dStateFEEdge[nFiniteElements-1][Right];
	}
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::DifferentiateNodeToNode(
	const double * dDataNode,
	double * dDiffNode,
	bool fZeroBoundaries
) {
	const int Left = 0;
	const int Right = 1;

	// Number of radial elements
	int nRElements = m_model.GetGrid()->GetRElements();

	// Zero the output
	memset(dDiffNode, 0, nRElements * sizeof(double));

	// Interpolate nodes to finite-element edges
	InterpolateNodeToFEEdges(dDataNode, fZeroBoundaries);

	// Calculate derivative at each node
	for (int k = 0; k < nRElements; k++) {
		int a = k / m_nVerticalOrder;
		int m = k % m_nVerticalOrder;

		int lBegin = a * m_nVerticalOrder;

		// Calculate derivatives due to internal bits
		for (int l = 0; l < m_nVerticalOrder; l++) {
			dDiffNode[k] += m_dDiffNodeToNode[m][l] * dDataNode[lBegin+l];
		}

		// Calculate derivatives due to interfaces
		dDiffNode[k] -=
			m_dDiffReconsPolyNode[m]
			* m_dStateFEEdge[a+1][Left];
		dDiffNode[k] +=
			m_dDiffReconsPolyNode[m_nVerticalOrder - m - 1]
			* m_dStateFEEdge[a][Right];
	}
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::DifferentiateNodeToREdge(
	const double * dDataNode,
	double * dDiffREdge,
	bool fZeroBoundaries
) {
	const int Left = 0;
	const int Right = 1;

	// Number of radial elements
	int nRElements = m_model.GetGrid()->GetRElements();

	int nFiniteElements = nRElements / m_nVerticalOrder;

	// Zero the output
	memset(dDiffREdge, 0, (nRElements+1) * sizeof(double));

	// Handle the single element case
	if (nFiniteElements == 1) {
		for (int k = 0; k < m_nVerticalOrder; k++) {
		for (int l = 0; l <= m_nVerticalOrder; l++) {
			dDiffREdge[l] +=
				m_dDiffNodeToREdge[l][k]
				* dDataNode[k];
		}
		}

		return;
	}

	if (nFiniteElements == 2) {
		_EXCEPTIONT("UNIMPLEMENTED: Still working on this...");
	}

	if (fZeroBoundaries) {
		_EXCEPTIONT("UNIMPLEMENTED: Still working on transitioning this...");
	}

	// Interpolate nodes to finite-element edges
	InterpolateNodeToFEEdges(dDataNode, fZeroBoundaries);

	// Calculate derivatives at interfaces
	for (int k = 0; k < nRElements; k++) {
		int a = k / m_nVerticalOrder;
		int m = k % m_nVerticalOrder;

		int lBegin = a * m_nVerticalOrder;

		// Push value of each node onto interface derivatives
		for (int l = 0; l <= m_nVerticalOrder; l++) {
			dDiffREdge[lBegin+l] +=
				m_dDiffNodeToREdge[l][m]
				* dDataNode[k];
		}
	}

	// Calculate derivatives due to interfaces
	for (int a = 0; a < nFiniteElements; a++) {
		int lBegin = a * m_nVerticalOrder;

		for (int l = 0; l <= m_nVerticalOrder; l++) {
			dDiffREdge[lBegin+l] -=
				m_dDiffReconsPolyREdge[l]
				* m_dStateFEEdge[a+1][Left];
			dDiffREdge[lBegin+l] +=
				m_dDiffReconsPolyREdge[m_nVerticalOrder - l]
				* m_dStateFEEdge[a][Right];
		}
	}

	// Scale interior finite element edges
	for (int a = 1; a < nFiniteElements; a++) {
		dDiffREdge[a * m_nVerticalOrder] *= 0.5;
	}

#pragma message "Doesn't work with ReconstructionFunctionType = 1"
/*
	// Compute interior derivatives
	{
		int kBegin = m_nVerticalOrder;
		int kLast = (nFiniteElements - 1) * m_nVerticalOrder;

		for (int k = kBegin; k <= kLast; k++) {
			int a = k / m_nVerticalOrder;
			int l = k % m_nVerticalOrder;
			if (k == kLast) {
				a--;
				l = m_nVerticalOrder;
			}

			int lPrev = (a-1) * m_nVerticalOrder;

			for (int m = 0; m < 3 * m_nVerticalOrder; m++) {
				dDiffREdge[k] +=
					m_dDiffNodeToREdgeAmal[l][m]
					* dDataNode[lPrev + m];
			}
		}
	}

	// Compute derivatives at left boundary
	for (int l = 0; l < m_nVerticalOrder; l++) {
	for (int m = 0; m < 2 * m_nVerticalOrder; m++) {
		dDiffREdge[l] +=
			m_dDiffNodeToREdgeLeft[l][m]
			* dDataNode[m];
	}
	}

	// Compute derivatives at right boundary
	{
		int lBegin = (nFiniteElements - 1) * m_nVerticalOrder;
		int lEnd = lBegin + m_nVerticalOrder;
		for (int l = 1; l <= m_nVerticalOrder; l++) {
		for (int m = 0; m < 2 * m_nVerticalOrder; m++) {
			dDiffREdge[lBegin + l] +=
				m_dDiffNodeToREdgeRight[l][m]
				* dDataNode[lBegin - m_nVerticalOrder + m];
		}
		}
	}
*/
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::DifferentiateREdgeToNode(
	const double * dDataREdge,
	double * dDiffNode
) {
	// Number of radial elements
	int nRElements = m_model.GetGrid()->GetRElements();

	int nFiniteElements = nRElements / m_nVerticalOrder;

	// Zero the data
	memset(dDiffNode, 0, nRElements * sizeof(double));

	// Loop through all nodes
	for (int k = 0; k < nRElements; k++) {

		// Differentiate from neighboring interface values
		int a = k / m_nVerticalOrder;
		int m = k % m_nVerticalOrder;
		int lBegin = a * m_nVerticalOrder;

		for (int l = 0; l <= m_nVerticalOrder; l++) {
			dDiffNode[k] +=
				  m_dDiffREdgeToNode[m][l]
				* dDataREdge[lBegin + l];
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::DifferentiateREdgeToREdge(
	const double * dDataREdge,
	double * dDiffREdge
) {
	// Number of radial elements
	int nRElements = m_model.GetGrid()->GetRElements();

	int nFiniteElements = nRElements / m_nVerticalOrder;

	// Zero the data
	memset(dDiffREdge, 0, (nRElements+1) * sizeof(double));

	// Apply all interfaces values to all interfaces within element
	for (int a = 0; a < nFiniteElements; a++) {
	for (int l = 0; l <= m_nVerticalOrder; l++) {

		int lBegin = a * m_nVerticalOrder;

		for (int m = 0; m <= m_nVerticalOrder; m++) {
			dDiffREdge[lBegin + l] +=
				  m_dDiffREdgeToREdge[l][m]
				* dDataREdge[lBegin + m];
		}
	}
	}

	// Halve interior element interface values
	for (int a = 1; a < nFiniteElements; a++) {
		dDiffREdge[a * m_nVerticalOrder] *= 0.5;
	}
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::DiffDiffREdgeToREdge(
	const double * dDataREdge,
	double * dDiffDiffREdge
) {
	// Number of radial elements
	int nRElements = m_model.GetGrid()->GetRElements();

	int nFiniteElements = nRElements / m_nVerticalOrder;

	// Zero the data
	memset(dDiffDiffREdge, 0, (nRElements+1) * sizeof(double));

	// Apply all interfaces values to all interfaces within element
	for (int a = 0; a < nFiniteElements; a++) {
	for (int l = 0; l <= m_nVerticalOrder; l++) {

		int lBegin = a * m_nVerticalOrder;

		for (int m = 0; m <= m_nVerticalOrder; m++) {
			dDiffDiffREdge[lBegin + l] +=
				  m_dDiffDiffREdgeToREdge[l][m]
				* dDataREdge[lBegin + m];
		}
	}
	}

	// Halve interior element interface values
	for (int a = 1; a < nFiniteElements; a++) {
		dDiffDiffREdge[a * m_nVerticalOrder] *= 0.5;
	}
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::SetupReferenceColumn(
	int iA,
	int iB,
	const DataMatrix<double> & dataTopography,
	const GridData4D & dataRefNode,
	const GridData4D & dataInitialNode,
	const GridData4D & dataRefREdge,
	const GridData4D & dataInitialREdge,
	const GridData3D & dataExnerNode,
	const GridData3D & dataDiffExnerNode,
	const GridData3D & dataExnerREdge,
	const GridData3D & dataDiffExnerREdge
) {

	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int TIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	// Get a copy of the grid
	Grid * pGrid = m_model.GetGrid();

	// Store domain height
	m_dDomainHeight = pGrid->GetZtop() - dataTopography[iA][iB];

	// Construct column of state variables
	int ix = 0;

	// Copy over Theta
	if (pGrid->GetVarLocation(TIx) == DataLocation_REdge) {
		for (int k = 0; k <= pGrid->GetRElements(); k++) {
			m_dColumnState[ix] = dataInitialREdge[TIx][k][iA][iB];
			ix++;
		}
	} else {
		for (int k = 0; k < pGrid->GetRElements(); k++) {
			m_dColumnState[ix] = dataInitialNode[TIx][k][iA][iB];
			ix++;
		}
	}

	// Copy over W
	if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
		for (int k = 0; k <= pGrid->GetRElements(); k++) {
			m_dColumnState[ix] = dataInitialREdge[WIx][k][iA][iB];
			ix++;
		}
	} else {
		for (int k = 0; k < pGrid->GetRElements(); k++) {
			m_dColumnState[ix] = dataInitialNode[WIx][k][iA][iB];
			ix++;
		}
	}

	// Copy over rho
	for (int k = 0; k < pGrid->GetRElements(); k++) {
		m_dColumnState[ix] = dataInitialNode[RIx][k][iA][iB];
		ix++;
	}

	// Construct reference column
	for (int k = 0; k < pGrid->GetRElements(); k++) {
		m_dStateRefNode[RIx][k] = dataRefNode[RIx][k][iA][iB];
		m_dStateRefNode[TIx][k] = dataRefNode[TIx][k][iA][iB];
	}
	for (int k = 0; k <= pGrid->GetRElements(); k++) {
		m_dStateRefREdge[RIx][k] = dataRefREdge[RIx][k][iA][iB];
		m_dStateRefREdge[TIx][k] = dataRefREdge[TIx][k][iA][iB];
	}

	// Build the Exner pressure reference
	for (int k = 0; k < pGrid->GetRElements(); k++) {
		m_dExnerRefNode[k] = dataExnerNode[k][iA][iB];
		m_dDiffExnerRefNode[k] = dataDiffExnerNode[k][iA][iB];
	}

	for (int k = 0; k <= pGrid->GetRElements(); k++) {
		m_dExnerRefREdge[k] = dataExnerREdge[k][iA][iB];
		m_dDiffExnerRefREdge[k] = dataDiffExnerREdge[k][iA][iB];
	}
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::StepExplicit(
	int iDataInitial,
	int iDataUpdate,
	const Time & time,
	double dDeltaT
) {
	// Get a copy of the grid
	Grid * pGrid = m_model.GetGrid();

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int TIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	// Number of elements
	int nRElements = pGrid->GetRElements();

	// Number of finite elements in the vertical
	int nFiniteElements = nRElements / m_nVerticalOrder;

	// Reset the reference state
	memset(m_dStateRefNode[WIx],  0,  nRElements   *sizeof(double));
	memset(m_dStateRefREdge[WIx], 0, (nRElements+1)*sizeof(double));

	// Perform local update
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatch * pPatch = pGrid->GetActivePatch(n);

		const PatchBox & box = pPatch->GetPatchBox();

		// Data
		const GridData4D & dataRefNode =
			pPatch->GetReferenceState(DataLocation_Node);

		const GridData4D & dataInitialNode =
			pPatch->GetDataState(iDataInitial, DataLocation_Node);

		GridData4D & dataUpdateNode =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		const GridData4D & dataRefREdge =
			pPatch->GetReferenceState(DataLocation_REdge);

		const GridData4D & dataInitialREdge =
			pPatch->GetDataState(iDataInitial, DataLocation_REdge);

		GridData4D & dataUpdateREdge =
			pPatch->GetDataState(iDataUpdate, DataLocation_REdge);

		// Auxiliary data storing Exner pressure
		const GridData3D & dataExnerNode =
			pPatch->GetVerticalDynamicsAuxData(0, DataLocation_Node);

		const GridData3D & dataExnerREdge =
			pPatch->GetVerticalDynamicsAuxData(0, DataLocation_REdge);

		const GridData3D & dataDiffExnerNode =
			pPatch->GetVerticalDynamicsAuxData(1, DataLocation_Node);

		const GridData3D & dataDiffExnerREdge =
			pPatch->GetVerticalDynamicsAuxData(1, DataLocation_REdge);

		// Pointwise topography
		const DataMatrix<double> & dataTopography = pPatch->GetTopography();

#pragma message "Perform update as in StepImplicit to reduce computational cost"
		// Loop over all nodes
		for (int i = box.GetAInteriorBegin(); i < box.GetAInteriorEnd(); i++) {
		for (int j = box.GetBInteriorBegin(); j < box.GetBInteriorEnd(); j++) {

			// Update thermodynamic variables
			if (m_fFullyExplicit) {

				int iA = i;
				int iB = j;

				SetupReferenceColumn(
					iA, iB,
					dataTopography,
					dataRefNode,
					dataInitialNode,
					dataRefREdge,
					dataInitialREdge,
					dataExnerNode,
					dataDiffExnerNode,
					dataExnerREdge,
					dataDiffExnerREdge);

                Evaluate(m_dColumnState, m_dSoln);

				// Apply updated state
				int ix = 0;

				// Apply update to theta
				if (pGrid->GetVarLocation(TIx) == DataLocation_REdge) {
					for (int k = 0; k <= pGrid->GetRElements(); k++) {
						dataUpdateREdge[TIx][k][iA][iB] -=
							dDeltaT * m_dSoln[ix++];
					}
				} else {
					for (int k = 0; k < pGrid->GetRElements(); k++) {
						dataUpdateNode[TIx][k][iA][iB] -=
							dDeltaT * m_dSoln[ix++];
					}
				}

				// Apply update to W
				if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
					for (int k = 0; k <= pGrid->GetRElements(); k++) {
						dataUpdateREdge[WIx][k][iA][iB] -=
							dDeltaT * m_dSoln[ix++];
					}
				} else {
					for (int k = 0; k < pGrid->GetRElements(); k++) {
						dataUpdateNode[WIx][k][iA][iB] -=
							dDeltaT * m_dSoln[ix++];
					}
				}

				// Apply update to rho
				for (int k = 0; k < pGrid->GetRElements(); k++) {
					dataUpdateNode[RIx][k][iA][iB] -=
						dDeltaT * m_dSoln[ix++];
				}
			}

			// Store W in State structure
			if (pGrid->GetVarLocation(WIx) == DataLocation_Node) {
				for (int k = 0; k < nRElements; k++) {
					m_dStateNode[WIx][k] = dataInitialNode[WIx][k][i][j];
				}
			} else {
				for (int k = 0; k <= nRElements; k++) {
					m_dStateREdge[WIx][k] = dataInitialREdge[WIx][k][i][j];
				}
			}

			// U and V on model levels
			if (pGrid->GetVarLocation(RIx) == DataLocation_Node) {

				// Obtain W on model levels
				if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
					InterpolateREdgeToNode(
						m_dStateREdge[WIx],
						m_dStateRefREdge[WIx],
						m_dStateNode[WIx],
						m_dStateRefNode[WIx]);
				}

				// Store U and V on model levels
				for (int k = 0; k < nRElements; k++) {
					m_dStateNode[UIx][k] = dataInitialNode[UIx][k][i][j];
					m_dStateNode[VIx][k] = dataInitialNode[VIx][k][i][j];
				}

				// Calculate update to U
				DifferentiateNodeToNode(
					m_dStateNode[UIx],
					m_dStateAuxDiff);

				for (int k = 0; k < nRElements; k++) {
					dataUpdateNode[UIx][k][i][j] -=
						dDeltaT * m_dStateNode[WIx][k] * m_dStateAuxDiff[k];
				}

				// Calculate update to V
				DifferentiateNodeToNode(
					m_dStateNode[VIx],
					m_dStateAuxDiff);

				for (int k = 0; k < nRElements; k++) {
					dataUpdateNode[VIx][k][i][j] -=
						dDeltaT * m_dStateNode[WIx][k] * m_dStateAuxDiff[k];
				}

			// U and V on model interfaces
			} else {

				// Store U and V on model interfaces
				for (int k = 0; k <= nRElements; k++) {
					m_dStateREdge[UIx][k] = dataInitialREdge[UIx][k][i][j];
					m_dStateREdge[VIx][k] = dataInitialREdge[VIx][k][i][j];
				}

				// Obtain W on model interfaces
				if (pGrid->GetVarLocation(WIx) == DataLocation_Node) {
					InterpolateNodeToREdge(
						m_dStateNode[WIx],
						m_dStateRefNode[WIx],
						m_dStateREdge[WIx],
						m_dStateRefREdge[WIx]);

					// Set boundary conditions
					m_dStateREdge[WIx][0] = 0.0;
					m_dStateREdge[WIx][nRElements] = 0.0;
				}

				// Calculate update to U
				DifferentiateREdgeToREdge(
					m_dStateREdge[UIx],
					m_dStateAuxDiff);

				for (int k = 1; k < nRElements; k++) {
					dataUpdateREdge[UIx][k][i][j] -=
						dDeltaT * m_dStateREdge[WIx][k] * m_dStateAuxDiff[k];
				}

				// Calculate update to U
				DifferentiateREdgeToREdge(
					m_dStateREdge[VIx],
					m_dStateAuxDiff);

				for (int k = 1; k < nRElements; k++) {
					dataUpdateREdge[VIx][k][i][j] -=
						dDeltaT * m_dStateREdge[WIx][k] * m_dStateAuxDiff[k];
				}
			}

			// Vertical transport of vertical momentum (W on model levels)
			if (pGrid->GetVarLocation(WIx) == DataLocation_Node) {

				// Differentiate W on model levels
				DifferentiateNodeToNode(
					m_dStateNode[WIx],
					m_dStateAuxDiff);

				for (int k = 0; k < nRElements; k++) {
					dataUpdateNode[WIx][k][i][j] -=
						dDeltaT * m_dStateNode[WIx][k] * m_dStateAuxDiff[k];
				}

			// Vertical transport of vertical momentum (W at interfaces)
			} else {

				// Differentiate W on model interfaces
				DifferentiateREdgeToREdge(
					m_dStateREdge[WIx],
					m_dStateAuxDiff);

				for (int k = 1; k < nRElements; k++) {
					dataUpdateREdge[WIx][k][i][j] -=
						dDeltaT * m_dStateREdge[WIx][k] * m_dStateAuxDiff[k];
				}
			}
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::BootstrapJacobian() {

	static const double Epsilon = 1.0e-5;

	int nDim = m_dColumnState.GetRows();

	DataMatrix<double> dJacobian;
	dJacobian.Initialize(nDim, nDim);

	DataVector<double> dJC;
	dJC.Initialize(nDim);

	DataVector<double> dG;
	dG.Initialize(nDim);

	DataVector<double> dJCref;
	dJCref.Initialize(nDim);

	Evaluate(m_dColumnState, dJCref);

	for (int i = 0; i < nDim; i++) {
		dG = m_dColumnState;
		dG[i] = dG[i] + Epsilon;

		Evaluate(dG, dJC);

		for (int j = 0; j < nDim; j++) {
			dJacobian[j][i] = (dJC[j] - dJCref[j]) / Epsilon;
		}
	}

	std::cout << "DeltaT: " << m_dDeltaT << std::endl;

	FILE * fp;
	fp = fopen("DGRef.txt", "w");
	for (int i = 0; i < nDim; i++) {
		for (int j = 0; j < nDim; j++) {
			fprintf(fp, "%1.15e", dJacobian[i][j]);
			if (j != nDim-1) {
				fprintf(fp, " ");
			}
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	fp = fopen("G.txt", "w");
	for (int i = 0; i < nDim; i++) {
		fprintf(fp, "%1.15e", dJCref[i]);
		if (i != nDim-1) {
			fprintf(fp, "\n");
		}
	}
	fclose(fp);

	BuildJacobianF(m_dSoln, dJacobian);

	fp = fopen("DG.txt", "w");
	for (int i = 0; i < nDim; i++) {
		for (int j = 0; j < nDim; j++) {
			fprintf(fp, "%1.15e", dJacobian[i][j]);
			if (j != nDim-1) {
				fprintf(fp, " ");
			}
		}
		fprintf(fp, "\n");
	}
	fclose(fp);


	_EXCEPTION();
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::StepImplicit(
	int iDataInitial,
	int iDataUpdate,
	const Time & time,
	double dDeltaT
) {
	// If fully explicit do nothing
	if (m_fFullyExplicit) {
		return;
	}

	// Get a copy of the grid
	Grid * pGrid = m_model.GetGrid();

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int TIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	// Store timestep size
	m_dDeltaT = dDeltaT;

	// Perform local update
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatch * pPatch = pGrid->GetActivePatch(n);

		const PatchBox & box = pPatch->GetPatchBox();

		// Contravariant metric components
		const DataMatrix4D<double> & dContraMetricA =
			pPatch->GetContraMetricA();
		const DataMatrix4D<double> & dContraMetricB =
			pPatch->GetContraMetricB();

		// Data
		const GridData4D & dataRefNode =
			pPatch->GetReferenceState(DataLocation_Node);

		const GridData4D & dataInitialNode =
			pPatch->GetDataState(iDataInitial, DataLocation_Node);

		GridData4D & dataUpdateNode =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		const GridData4D & dataRefREdge =
			pPatch->GetReferenceState(DataLocation_REdge);

		const GridData4D & dataInitialREdge =
			pPatch->GetDataState(iDataInitial, DataLocation_REdge);

		GridData4D & dataUpdateREdge =
			pPatch->GetDataState(iDataUpdate, DataLocation_REdge);

		// Auxiliary data storing Exner pressure
		const GridData3D & dataExnerNode =
			pPatch->GetVerticalDynamicsAuxData(0, DataLocation_Node);

		const GridData3D & dataExnerREdge =
			pPatch->GetVerticalDynamicsAuxData(0, DataLocation_REdge);

		const GridData3D & dataDiffExnerNode =
			pPatch->GetVerticalDynamicsAuxData(1, DataLocation_Node);

		const GridData3D & dataDiffExnerREdge =
			pPatch->GetVerticalDynamicsAuxData(1, DataLocation_REdge);

		// Pointwise topography
		const DataMatrix<double> & dataTopography = pPatch->GetTopography();

		// Number of finite elements
		int nAElements =
			box.GetAInteriorWidth() / m_nHorizontalOrder;
		int nBElements =
			box.GetBInteriorWidth() / m_nHorizontalOrder;

		// Loop over all nodes, but only perform calculation on shared
		// nodes once
		for (int a = 0; a < nAElements; a++) {
		for (int b = 0; b < nBElements; b++) {

			int iEnd;
			int jEnd;

			if (a == nAElements-1) {
				iEnd = m_nHorizontalOrder;
			} else {
				iEnd = m_nHorizontalOrder-1;
			}

			if (b == nBElements-1) {
				jEnd = m_nHorizontalOrder;
			} else {
				jEnd = m_nHorizontalOrder-1;
			}

		for (int i = 0; i < iEnd; i++) {
		for (int j = 0; j < jEnd; j++) {

			int iA = box.GetAInteriorBegin() + a * m_nHorizontalOrder + i;
			int iB = box.GetBInteriorBegin() + b * m_nHorizontalOrder + j;

			SetupReferenceColumn(
				iA, iB,
				dataTopography,
				dataRefNode,
				dataInitialNode,
				dataRefREdge,
				dataInitialREdge,
				dataExnerNode,
				dataDiffExnerNode,
				dataExnerREdge,
				dataDiffExnerREdge);

#ifdef USE_PETSC
			// Use PetSc to solve
			double * dX;
			VecGetArray(m_vecX, &dX);
			memcpy(dX, m_dColumnState, m_nColumnStateSize * sizeof(double));
			VecRestoreArray(m_vecX, &dX);

			// Solve
			PetscErrorCode ierr;
			SNESSolve(m_snes, NULL, m_vecX);

			SNESConvergedReason reason;
			SNESGetConvergedReason(m_snes, &reason);
			if ((reason < 0) && (reason != (-5))) {
				_EXCEPTION1("PetSc solver failed to converge (%i)", reason);
			}
/*
			{
				// Evaluate
				DataVector<double> dEval;
				dEval.Initialize(m_dColumnState.GetRows());
				VecGetArray(m_vecX, &dX);
				Evaluate(dX, dEval);

				nAvgValues++;
				double dEvalNorm = 0.0;
				for (int n = 0; n < dEval.GetRows(); n++) {
					dEvalNorm += dEval[n] * dEval[n];
					std::cout << dEval[n] << std::endl;
				}
				dAvgNorm += sqrt(dEvalNorm / dEval.GetRows());

				VecRestoreArray(m_vecX, &dX);
				_EXCEPTION();
			}
*/
			VecGetArray(m_vecX, &dX);
			memcpy(m_dSoln, dX, m_nColumnStateSize * sizeof(double));
			VecRestoreArray(m_vecX, &dX);
#endif
#ifdef USE_GMRES
			// Use Jacobian-Free Newton-Krylov to solve
			m_dSoln = m_dColumnState;

			//BootstrapJacobian();

			double dError =
				PerformJFNK_NewtonStep_Safe(
				//PerformBICGSTAB_NewtonStep_Safe(
					m_dSoln,
					m_dSoln.GetRows(),
					1.0e-8);

			// DEBUG (check for NANs in output)
			if (!(m_dSoln[0] == m_dSoln[0])) {
                DataVector<double> dEval;
                dEval.Initialize(m_dColumnState.GetRows());
                Evaluate(m_dSoln, dEval);

                for (int p = 0; p < dEval.GetRows(); p++) {
                    printf("%1.15e %1.15e %1.15e\n",
						dEval[p], m_dSoln[p] - m_dColumnState[p], m_dColumnState[p]);
                }
				for (int p = 0; p < m_dExnerRefREdge.GetRows(); p++) {
					printf("%1.15e %1.15e\n",
						m_dExnerRefREdge[p], dataRefREdge[RIx][p][iA][iB]);
				}
                _EXCEPTIONT("Inversion failure");
            }

#endif
#ifdef USE_DIRECTSOLVE_APPROXJ
			static const double Epsilon = 1.0e-5;

			// Prepare the column
			PrepareColumn(m_dColumnState);

			// Build the F vector
			BuildF(m_dColumnState, m_dSoln);

			DataVector<double> dJC;
			dJC.Initialize(m_dColumnState.GetRows());

			DataVector<double> dG;
			dG.Initialize(m_dColumnState.GetRows());

			DataVector<double> dJCref;
			dJCref.Initialize(m_dColumnState.GetRows());

			Evaluate(m_dColumnState, dJCref);

			for (int i = 0; i < m_dColumnState.GetRows(); i++) {
				dG = m_dColumnState;
				dG[i] = dG[i] + Epsilon;

				Evaluate(dG, dJC);

				for (int j = 0; j < m_dColumnState.GetRows(); j++) {
					m_matJacobianF[i][j] = (dJC[j] - dJCref[j]) / Epsilon;
				}
			}

			// Use direct solver
			LAPACK::DGESV(m_matJacobianF, m_dSoln, m_vecIPiv);

			for (int k = 0; k < m_dSoln.GetRows(); k++) {
				m_dSoln[k] = m_dColumnState[k] - m_dSoln[k];
			}
#endif
#ifdef USE_DIRECTSOLVE
			// Prepare the column
			PrepareColumn(m_dColumnState);

			// Build the F vector
			BuildF(m_dColumnState, m_dSoln);

			// Build the Jacobian
			BuildJacobianF(m_dColumnState, m_matJacobianF);

			// Use direct solver
			LAPACK::DGESV(m_matJacobianF, m_dSoln, m_vecIPiv);

			for (int k = 0; k < m_dSoln.GetRows(); k++) {
				m_dSoln[k] = m_dColumnState[k] - m_dSoln[k];
			}
#endif

			// Apply updated state
			int ix = 0;

			// Apply updated state to theta
			if (pGrid->GetVarLocation(TIx) == DataLocation_REdge) {
				for (int k = 0; k <= pGrid->GetRElements(); k++) {
					dataUpdateREdge[TIx][k][iA][iB] = m_dSoln[ix++];
				}
			} else {
				for (int k = 0; k < pGrid->GetRElements(); k++) {
					dataUpdateNode[TIx][k][iA][iB] = m_dSoln[ix++];
				}
			}

			// Copy over W
			if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
				for (int k = 0; k <= pGrid->GetRElements(); k++) {
					dataUpdateREdge[WIx][k][iA][iB] = m_dSoln[ix++];
				}
			} else {
				for (int k = 0; k < pGrid->GetRElements(); k++) {
					dataUpdateNode[WIx][k][iA][iB] = m_dSoln[ix++];
				}
			}

			// Copy over rho
			for (int k = 0; k < pGrid->GetRElements(); k++) {
				dataUpdateNode[RIx][k][iA][iB] = m_dSoln[ix++];
			}

#pragma message "Vertical pressure gradient influence on horizontal velocities"
/*
			// Compute vertical pressure gradient on nodes for
			// potential temperature (theta) on interfaces
			if (pGrid->GetVarLocation(TIx) == DataLocation_REdge) {
				InterpolateREdgeToNode(
					&(m_dSoln[m_ixTBegin]),
					m_dStateRefREdge[TIx],
					m_dStateNode[TIx],
					m_dStateRefNode[TIx]);

			// potential temperature (theta) on levels
			} else {
				for (int k = 0; k < pGrid->GetRElements(); k++) {
					m_dStateNode[TIx][k] = dataUpdateNode[TIx][k][iA][iB];
				}
			}

			// Exner function perturbation on nodes
			for (int k = 0; k < pGrid->GetRElements(); k++) {
				m_dExnerPertNode[k] =
					phys.GetCp() * exp(phys.GetR() / phys.GetCv()
						* log(phys.GetR() / phys.GetP0()
							* dataUpdateNode[RIx][k][iA][iB]
							* m_dStateNode[TIx][k]));

				m_dExnerPertNode[k] -= m_dExnerRefNode[k];
			}

			// Differentiate Exner function
			DifferentiateNodeToNode(
				m_dExnerPertNode,
				m_dDiffExnerPertNode);

			// Apply Exner gradient term to horizontal velocities
			for (int k = 0; k < pGrid->GetRElements(); k++) {
				double dPressureTerm = dDeltaT * m_dStateNode[TIx][k] * (
					m_dDiffExnerRefNode[k] + m_dDiffExnerPertNode[k]);

				dataUpdateNode[UIx][k][iA][iB] -=
					dContraMetricA[k][iA][iB][2] * dPressureTerm;

				dataUpdateNode[VIx][k][iA][iB] -=
					dContraMetricB[k][iA][iB][2] * dPressureTerm;
			}
*/
/*
            // Compute vertical pressure gradient on nodes for
            // potential temperature (theta) on interfaces
            if (pGrid->GetVarLocation(TIx) == DataLocation_REdge) {
                InterpolateNodeToREdge(
                    &(m_dSoln[m_ixRBegin]),
                    m_dStateRefNode[RIx],
                    m_dStateREdge[RIx],
                    m_dStateRefREdge[RIx]
                );

                for (int k = 0; k <= pGrid->GetRElements(); k++) {
                    m_dStateAux[k] =
                        phys.PressureFromRhoTheta(
                            m_dStateREdge[RIx][k] * m_dSoln[m_ixTBegin+k]);
                }

                DifferentiateREdgeToNode(
                    m_dStateAux,
                    m_dStateAuxDiff);

            // potential temperature (theta) on levels
            } else {
                for (int k = 0; k < pGrid->GetRElements(); k++) {
                    m_dStateAux[k] =
                        phys.PressureFromRhoTheta(
                            m_dSoln[m_ixRBegin+k] * m_dSoln[m_ixTBegin+k]);
                }

                DifferentiateNodeToNode(
                    m_dStateAux,
                    m_dStateAuxDiff);
            }

            // Apply vertical pressure gradient term to horizontal velocities
            for (int k = 0; k < pGrid->GetRElements(); k++) {
                double dUpdateGradP =
                    dDeltaT
                    * m_dStateAuxDiff[k]
                    / dataUpdateNode[RIx][k][iA][iB];

                dataUpdateNode[UIx][k][iA][iB] -=
                    dUpdateGradP * dContraMetricA[k][iA][iB][2];

                dataUpdateNode[VIx][k][iA][iB] -=
                    dUpdateGradP * dContraMetricB[k][iA][iB][2];
            }
*/
		}
		}

		}
		}

		// Copy over new state on shared nodes (edges of constant alpha)
		for (int a = 1; a < nAElements; a++) {
			int iA = box.GetAInteriorBegin() + a * m_nHorizontalOrder - 1;

			for (int b = 0; b < nBElements; b++) {

				// Top element contains more information
				int jEnd;
				if (b == nBElements-1) {
					jEnd = m_nHorizontalOrder;
				} else {
					jEnd = m_nHorizontalOrder-1;
				}

				// Loop along edges of constant alpha
				for (int j = 0; j < jEnd; j++) {

					int iB = box.GetBInteriorBegin() + b * m_nHorizontalOrder + j;

					for (int k = 0; k < pGrid->GetRElements(); k++) {
						dataUpdateNode[UIx][k][iA][iB]
							= dataUpdateNode[UIx][k][iA+1][iB];
						dataUpdateNode[VIx][k][iA][iB]
							= dataUpdateNode[VIx][k][iA+1][iB];
						dataUpdateNode[TIx][k][iA][iB]
							= dataUpdateNode[TIx][k][iA+1][iB];
						dataUpdateNode[WIx][k][iA][iB]
							= dataUpdateNode[WIx][k][iA+1][iB];
						dataUpdateNode[RIx][k][iA][iB]
							= dataUpdateNode[RIx][k][iA+1][iB];
					}

					for (int k = 0; k <= pGrid->GetRElements(); k++) {
						dataUpdateREdge[UIx][k][iA][iB]
							= dataUpdateREdge[UIx][k][iA+1][iB];
						dataUpdateREdge[VIx][k][iA][iB]
							= dataUpdateREdge[VIx][k][iA+1][iB];
						dataUpdateREdge[TIx][k][iA][iB]
							= dataUpdateREdge[TIx][k][iA+1][iB];
						dataUpdateREdge[WIx][k][iA][iB]
							= dataUpdateREdge[WIx][k][iA+1][iB];
						dataUpdateREdge[RIx][k][iA][iB]
							= dataUpdateREdge[RIx][k][iA+1][iB];
					}
				}
			}
		}

		// Copy over new state on shared nodes (edges of constant beta)
		for (int b = 1; b < nBElements; b++) {
		for (int i = box.GetAInteriorBegin(); i < box.GetAInteriorEnd(); i++) {
			int iB = box.GetBInteriorBegin() + b * m_nHorizontalOrder - 1;

			for (int k = 0; k < pGrid->GetRElements(); k++) {
				dataUpdateNode[UIx][k][i][iB]
					= dataUpdateNode[UIx][k][i][iB+1];
				dataUpdateNode[VIx][k][i][iB]
					= dataUpdateNode[VIx][k][i][iB+1];
				dataUpdateNode[TIx][k][i][iB]
					= dataUpdateNode[TIx][k][i][iB+1];
				dataUpdateNode[WIx][k][i][iB]
					= dataUpdateNode[WIx][k][i][iB+1];
				dataUpdateNode[RIx][k][i][iB]
					= dataUpdateNode[RIx][k][i][iB+1];
			}

			for (int k = 0; k <= pGrid->GetRElements(); k++) {
				dataUpdateREdge[UIx][k][i][iB]
					= dataUpdateREdge[UIx][k][i][iB+1];
				dataUpdateREdge[VIx][k][i][iB]
					= dataUpdateREdge[VIx][k][i][iB+1];
				dataUpdateREdge[TIx][k][i][iB]
					= dataUpdateREdge[TIx][k][i][iB+1];
				dataUpdateREdge[WIx][k][i][iB]
					= dataUpdateREdge[WIx][k][i][iB+1];
				dataUpdateREdge[RIx][k][i][iB]
					= dataUpdateREdge[RIx][k][i][iB+1];
			}
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::PrepareColumn(
	const double * dX
) {
	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int TIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	// Finite element grid spacing
	const Grid * pGrid = m_model.GetGrid();

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Number of radial elements
	const int nRElements = pGrid->GetRElements();

	// W on model interfaces
	if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
		memcpy(
			m_dStateREdge[WIx],
			&(dX[m_ixWBegin]),
			(nRElements+1) * sizeof(double));

		// W is needed on model levels
		if ((m_fMassFluxOnLevels) ||
			(pGrid->GetVarLocation(TIx) == DataLocation_Node)
		) {
			InterpolateREdgeToNode(
				&(dX[m_ixWBegin]),
				m_dStateRefREdge[WIx],
				m_dStateNode[WIx],
				m_dStateRefNode[WIx]);
		}

	// W on model levels
	} else {
		memcpy(
			m_dStateNode[WIx],
			&(dX[m_ixWBegin]),
			nRElements * sizeof(double));

		// W is needed on model interfaces
		if ((!m_fMassFluxOnLevels) ||
			(pGrid->GetVarLocation(TIx) == DataLocation_REdge)
		) {
			InterpolateNodeToREdge(
				&(dX[m_ixWBegin]),
				m_dStateRefNode[WIx],
				m_dStateREdge[WIx],
				m_dStateRefREdge[WIx],
				true);
		}
	}

	// T on model interfaces
	if (pGrid->GetVarLocation(TIx) == DataLocation_REdge) {
		memcpy(
			m_dStateREdge[TIx],
			&(dX[m_ixTBegin]),
			(nRElements+1) * sizeof(double));

		// Calculate update for T
		DifferentiateREdgeToREdge(
			&(dX[m_ixTBegin]),
			m_dDiffTheta);

		// T is needed on model levels
		if ((m_fExnerPressureOnLevels) ||
			(pGrid->GetVarLocation(WIx) == DataLocation_Node)
		) {
			InterpolateREdgeToNode(
				&(dX[m_ixTBegin]),
				m_dStateRefREdge[TIx],
				m_dStateNode[TIx],
				m_dStateRefNode[TIx]);
		}

	// T on model levels
	} else {
		memcpy(
			m_dStateNode[TIx],
			&(dX[m_ixTBegin]),
			nRElements * sizeof(double));

		// Calculate update for T
		DifferentiateNodeToNode(
			&(dX[m_ixTBegin]),
			m_dDiffTheta);

		// T is needed on model interfaces
		if ((!m_fExnerPressureOnLevels) ||
			(pGrid->GetVarLocation(WIx) == DataLocation_REdge)
		) {
			InterpolateNodeToREdge(
				&(dX[m_ixTBegin]),
				m_dStateRefNode[TIx],
				m_dStateREdge[TIx],
				m_dStateRefREdge[TIx]);
		}
	}

	// Rho on model interfaces
	if (pGrid->GetVarLocation(RIx) == DataLocation_REdge) {
		memcpy(
			m_dStateREdge[RIx],
			&(dX[m_ixRBegin]),
			(nRElements+1) * sizeof(double));

		// Rho is needed on model levels
		if ((m_fMassFluxOnLevels) ||
			(m_fExnerPressureOnLevels)
		) {
			InterpolateREdgeToNode(
				&(dX[m_ixRBegin]),
				m_dStateRefREdge[RIx],
				m_dStateNode[RIx],
				m_dStateRefNode[RIx]);
		}

	// Rho on model levels
	} else {
		memcpy(
			m_dStateNode[RIx],
			&(dX[m_ixRBegin]),
			nRElements * sizeof(double));

		// Rho is needed on model interfaces
		if ((!m_fMassFluxOnLevels) ||
			(!m_fExnerPressureOnLevels)
		) {
			InterpolateNodeToREdge(
				&(dX[m_ixRBegin]),
				m_dStateRefNode[RIx],
				m_dStateREdge[RIx],
				m_dStateRefREdge[RIx]);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::BuildF(
	const double * dX,
	double * dF
) {
	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int TIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	// Finite element grid spacing
	const Grid * pGrid = m_model.GetGrid();

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Number of radial elements
	const int nRElements = pGrid->GetRElements();

	// Grid spacing
	double dElementDeltaXi =
		  pGrid->GetREtaInterface(m_nVerticalOrder)
		- pGrid->GetREtaInterface(0);

	double dDeltaXi = dElementDeltaXi / static_cast<double>(m_nVerticalOrder);

	// Zero F
	memset(dF, 0, m_nColumnStateSize * sizeof(double));

	// T on model interfaces
	if (pGrid->GetVarLocation(TIx) == DataLocation_REdge) {
		for (int k = 0; k <= nRElements; k++) {
			dF[m_ixTBegin+k] = m_dStateREdge[WIx][k] * m_dDiffTheta[k];
		}

	// T on model levels
	} else {
		for (int k = 0; k < nRElements; k++) {
			dF[m_ixTBegin+k] = m_dStateNode[WIx][k] * m_dDiffTheta[k];
		}
	}

	// Calculate mass flux on model levels
	if (m_fMassFluxOnLevels) {
		for (int k = 0; k < nRElements; k++) {
			m_dMassFluxNode[k] =
				  m_dStateNode[RIx][k]
				* m_dStateNode[WIx][k];
		}

		// Calculate derivative of mass flux on interfaces
		if (pGrid->GetVarLocation(RIx) == DataLocation_REdge) {
			DifferentiateNodeToREdge(
				m_dMassFluxNode,
				m_dDiffMassFluxREdge,
				true);

		// Calculate derivative of mass flux on levels
		} else {
			DifferentiateNodeToNode(
				m_dMassFluxNode,
				m_dDiffMassFluxNode,
				true);
		}

	// Calculate mass flux on model interfaces
	} else {
		for (int k = 1; k < nRElements; k++) {
			m_dMassFluxREdge[k] =
				  m_dStateREdge[RIx][k] * m_dStateREdge[WIx][k];
/*
				  - 0.5 * 50.0 * fabs(m_dStateREdge[WIx][k])
				  		* (m_dStateNode[RIx][k] - m_dStateNode[RIx][k-1]);
*/
		}
		m_dMassFluxREdge[0] = 0.0;
		m_dMassFluxREdge[nRElements] = 0.0;

		// Calculate derivative of mass flux on interfaces
		if (pGrid->GetVarLocation(RIx) == DataLocation_REdge) {
			DifferentiateREdgeToREdge(
				m_dMassFluxNode,
				m_dDiffMassFluxREdge);

		// Calculate derivative of mass flux on levels
		} else {
			DifferentiateREdgeToNode(
				m_dMassFluxREdge,
				m_dDiffMassFluxNode);
		}
	}

	// Calculate Exner pressure on model levels
	if (m_fExnerPressureOnLevels) {
		for (int k = 0; k < nRElements; k++) {
			m_dExnerPertNode[k] =
				phys.GetCp() * exp(phys.GetR() / phys.GetCv()
					* log(phys.GetR() / phys.GetP0()
						* m_dStateNode[RIx][k]
						* m_dStateNode[TIx][k]));

			m_dExnerPertNode[k] -= m_dExnerRefNode[k];
		}

		// Calculate derivative of Exner pressure on interfaces
		if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
			DifferentiateNodeToREdge(
				m_dExnerPertNode,
				m_dDiffExnerPertREdge);

		// Calculate derivative of Exner pressure on levels
		} else {
			DifferentiateNodeToNode(
				m_dExnerPertNode,
				m_dDiffExnerPertNode);
		}

	// Calculate Exner pressure on model interfaces
	} else {
		for (int k = 0; k <= nRElements; k++) {
			m_dExnerPertREdge[k] =
				phys.GetCp() * exp(phys.GetR() / phys.GetCv()
					* log(phys.GetR() / phys.GetP0()
						* m_dStateREdge[RIx][k]
						* m_dStateREdge[TIx][k]));

			m_dExnerPertREdge[k] -= m_dExnerRefREdge[k];
		}

		// Calculate derivative of Exner pressure on interfaces
		if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
			DifferentiateREdgeToREdge(
				m_dExnerPertREdge,
				m_dDiffExnerPertREdge);

		// Calculate derivative of Exner pressure on levels
		} else {
			DifferentiateREdgeToNode(
				m_dExnerPertREdge,
				m_dDiffExnerPertNode);
		}
	}

	// Compute update to Rho on model interfaces
	if (pGrid->GetVarLocation(RIx) == DataLocation_REdge) {
		for (int k = 0; k <= nRElements; k++) {
			dF[m_ixRBegin+k] = m_dDiffMassFluxREdge[k];
		}

	// Compute update to Rho on model levels
	} else {
		for (int k = 0; k < nRElements; k++) {
			dF[m_ixRBegin+k] = m_dDiffMassFluxNode[k];
		}
	}

	// Compute update to W on model interfaces
	if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
		for (int k = 1; k < nRElements; k++) {
			double dTheta = m_dStateREdge[TIx][k];
			double dThetaPert = dTheta - m_dStateRefREdge[TIx][k];

			dF[m_ixWBegin+k] += 1.0 / (m_dDomainHeight * m_dDomainHeight) * (
				+ dThetaPert * m_dDiffExnerRefREdge[k]
				+ dTheta * m_dDiffExnerPertREdge[k]);
		}

		// If no vertical reference state is specified,
		// gravity must be included
		if (!pGrid->HasReferenceState()) {
			for (int k = 1; k < nRElements; k++) {
				dF[m_ixWBegin+k] += phys.GetG() / m_dDomainHeight;
			}
		}

	// Compute update to W on model levels
	} else {
		for (int k = 0; k < nRElements; k++) {
			double dTheta = m_dStateNode[TIx][k];
			double dThetaPert = dTheta - m_dStateRefNode[TIx][k];

			dF[m_ixWBegin+k] += 1.0 / (m_dDomainHeight * m_dDomainHeight) * (
				+ dThetaPert * m_dDiffExnerRefNode[k]
				+ dTheta * m_dDiffExnerPertNode[k]);
		}

		// If no vertical reference state is specified,
		// gravity must be included
		if (!pGrid->HasReferenceState()) {
			for (int k = 0; k < nRElements; k++) {
				dF[m_ixWBegin+k] += phys.GetG() / m_dDomainHeight;
			}
		}
	}

	// Apply diffusion to theta
	if ((pGrid->GetVarLocation(TIx) == DataLocation_REdge) &&
		(m_nHyperdiffusionOrder > 0)
	) {
		_EXCEPTION();
		double dScaledNu =
			m_dHyperdiffusionCoeff
			* exp(static_cast<double>(m_nHyperdiffusionOrder - 1)
				* log(dDeltaXi));

		DiffDiffREdgeToREdge(
			m_dStateREdge[TIx],
			m_dDiffTheta
		);

		for (int h = 2; h < m_nHyperdiffusionOrder; h += 2) {
			memcpy(
				m_dStateAux,
				m_dDiffTheta,
				(nRElements + 1) * sizeof(double));

			DiffDiffREdgeToREdge(
				m_dStateAux,
				m_dDiffTheta
			);
		}

		for (int k = 0; k <= nRElements; k++) {
			dF[m_ixTBegin+k] -=
				dScaledNu
				* fabs(m_dStateREdge[WIx][k])
				* m_dDiffTheta[k];
		}
	}

	// Construct the time-dependent component of the RHS
	for (int i = 0; i < m_nColumnStateSize; i++) {
		dF[i] += (dX[i] - m_dColumnState[i]) / m_dDeltaT;
	}
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::BuildJacobianF(
	const double * dX,
	double ** dDG
) {
	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int TIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	// Finite element grid
	const Grid * pGrid = m_model.GetGrid();

	if ((pGrid->GetVarLocation(UIx) != DataLocation_Node) ||
		(pGrid->GetVarLocation(VIx) != DataLocation_Node) ||
		(pGrid->GetVarLocation(TIx) != DataLocation_REdge) ||
		(pGrid->GetVarLocation(WIx) != DataLocation_REdge) ||
		(pGrid->GetVarLocation(RIx) != DataLocation_Node) ||
		(m_fMassFluxOnLevels) ||
		(!m_fExnerPressureOnLevels)
	) {
		_EXCEPTIONT("Not implemented");
	}

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Number of radial elements
	const int nRElements = pGrid->GetRElements();

	// Grid spacing
	double dElementDeltaXi =
		  pGrid->GetREtaInterface(m_nVerticalOrder)
		- pGrid->GetREtaInterface(0);

	double dDeltaXi = dElementDeltaXi / static_cast<double>(m_nVerticalOrder);

	// Zero DG
	memset(dDG[0], 0,
		m_nColumnStateSize * m_nColumnStateSize * sizeof(double));

	// dT_k/dT
	int nFiniteElements = nRElements / m_nVerticalOrder;
	for (int a = 0; a < nFiniteElements; a++) {
	for (int l = 0; l <= m_nVerticalOrder; l++) {
		int lBegin = a * m_nVerticalOrder;

		for (int k = 0; k <= m_nVerticalOrder; k++) {
			dDG[m_ixTBegin+lBegin+l][m_ixTBegin+lBegin+k] +=
				m_dDiffREdgeToREdge[k][l]
				* m_dStateREdge[WIx][lBegin+k];
		}
	}
	}

	for (int a = 1; a < nFiniteElements; a++) {
		int k = a * m_nVerticalOrder;
		int lBegin = k - m_nVerticalOrder;
		int lEnd   = k + m_nVerticalOrder + 1;
		for (int l = lBegin; l < lEnd; l++) {
			dDG[m_ixTBegin+l][m_ixTBegin+k] *= 0.5;
		}
	}

	// dT_k/dW
	for (int k = 0; k <= nRElements; k++) {
		dDG[m_ixWBegin+k][m_ixTBegin+k] = m_dDiffTheta[k];
	}

	if ((nFiniteElements == 1) || (nFiniteElements == 2)) {
		_EXCEPTIONT("UNIMPLEMENTED: At least three elements needed");
	}

	int kBegin;
	int kEnd;

	// dW_k/dT_l (bottom element)
	for (int k = 1; k < m_nVerticalOrder; k++) {
		double dRHSWCoeff = 
			1.0 / (m_dDomainHeight * m_dDomainHeight)
			* m_dStateREdge[TIx][k]
			* phys.GetR() / phys.GetCv();

		for (int m = 0; m < 2 * m_nVerticalOrder; m++) {
			double dTEntry =
				dRHSWCoeff 
				* m_dDiffNodeToREdgeLeft[k][m]
				* m_dExnerPertNode[m]
					/ m_dStateNode[TIx][m];

			int mx = m % m_nVerticalOrder;

			for (int l = 0; l <= m_nVerticalOrder; l++) {
				dDG[m_ixTBegin+l][m_ixWBegin+k] +=
					dTEntry * (*m_pInterpREdgeToNode)[mx][l];
			}

			dDG[m_ixRBegin+m][m_ixWBegin+k] +=
				dRHSWCoeff
				* m_dDiffNodeToREdgeLeft[k][m]
				* m_dExnerPertNode[m]
					/ m_dStateNode[RIx][m];
		}
	}

	// dW_k/dT_l and dW_k/dR_m (middle elements)
	int kLast = (nFiniteElements-1) * m_nVerticalOrder;
	for (int k = m_nVerticalOrder; k <= kLast; k++) {
		int kx = k % m_nVerticalOrder;
		int a = k / m_nVerticalOrder;
		if (k == kLast) {
			a--;
			kx = m_nVerticalOrder;
		}

		int lPrev = (a-1) * m_nVerticalOrder;
		int lBegin = lPrev + m_nVerticalOrder;

		double dRHSWCoeff = 
			1.0 / (m_dDomainHeight * m_dDomainHeight)
			* m_dStateREdge[TIx][k]
			* phys.GetR() / phys.GetCv();

		for (int m = 0; m < 3 * m_nVerticalOrder; m++) {

			int mx = m % m_nVerticalOrder;
			int ma = (m / m_nVerticalOrder) * m_nVerticalOrder;

			double dTEntry =
				dRHSWCoeff 
				* m_dDiffNodeToREdgeAmal[kx][m]
				* m_dExnerPertNode[lPrev + m]
					/ m_dStateNode[TIx][lPrev + m];

			for (int l = 0; l <= m_nVerticalOrder; l++) {
				dDG[m_ixTBegin+lPrev+ma+l][m_ixWBegin+k] +=
					dTEntry * (*m_pInterpREdgeToNode)[mx][l];
			}

			dDG[m_ixRBegin+lPrev+m][m_ixWBegin+k] +=
				dRHSWCoeff
				* m_dDiffNodeToREdgeAmal[kx][m]
				* m_dExnerPertNode[lPrev + m]
					/ m_dStateNode[RIx][lPrev + m];
		}
	}

	// dW_k/dT_l (top element)
	kBegin = (nFiniteElements-1) * m_nVerticalOrder;
	kEnd = kBegin + m_nVerticalOrder;
	for (int k = kBegin + 1; k < kEnd; k++) {
		int kx = k - kBegin;
		int lPrev = kBegin - m_nVerticalOrder;

		double dRHSWCoeff = 
			1.0 / (m_dDomainHeight * m_dDomainHeight)
			* m_dStateREdge[TIx][k]
			* phys.GetR() / phys.GetCv();

		for (int m = 0; m < 2 * m_nVerticalOrder; m++) {
			double dTEntry =
				dRHSWCoeff 
				* m_dDiffNodeToREdgeRight[kx][m]
				* m_dExnerPertNode[lPrev + m]
					/ m_dStateNode[TIx][lPrev + m];

			int mx = m % m_nVerticalOrder;
			int ma = (m / m_nVerticalOrder) * m_nVerticalOrder;

			for (int l = 0; l <= m_nVerticalOrder; l++) {
				dDG[m_ixTBegin+lPrev+ma+l][m_ixWBegin+k] +=
					dTEntry * (*m_pInterpREdgeToNode)[mx][l];
			}

			dDG[m_ixRBegin+lPrev+m][m_ixWBegin+k] +=
				dRHSWCoeff
				* m_dDiffNodeToREdgeRight[kx][m]
				* m_dExnerPertNode[lPrev + m]
					/ m_dStateNode[RIx][lPrev + m];
		}
	}

	// dW_k/dT
	for (int k = 1; k < nRElements; k++) {
		dDG[m_ixTBegin+k][m_ixWBegin+k] +=
			 1.0 / (m_dDomainHeight * m_dDomainHeight)
			 * m_dDiffExnerPertREdge[k];
	}

	// dRho
	for (int k = 0; k < nRElements; k++) {
		int a = k / m_nVerticalOrder;
		int l = k % m_nVerticalOrder;
		int lBegin = a * m_nVerticalOrder;

		// dRho_k/dW
		for (int m = 0; m <= m_nVerticalOrder; m++) {
			dDG[m_ixWBegin+lBegin+m][m_ixRBegin+k] +=
				m_dDiffREdgeToNode[l][m]
				* m_dStateREdge[RIx][lBegin+m];
		}

		// dRho_k/dRho
		if (a != 0) {
			int lPrev = lBegin - m_nVerticalOrder;
			for (int n = 0; n < m_nVerticalOrder; n++) {
				dDG[m_ixRBegin+lPrev+n][m_ixRBegin+k] +=
					m_dDiffREdgeToNode[l][0]
					* 0.5 * (*m_pInterpNodeToREdge)[m_nVerticalOrder][n]
					* m_dStateREdge[WIx][lBegin];
			}
		}

		for (int m = 0; m <= m_nVerticalOrder; m++) {
		for (int n = 0; n < m_nVerticalOrder; n++) {
			double dMult = 1.0;
			if ((m == 0) || (m == m_nVerticalOrder)) {
				dMult = 0.5;
			}
			dDG[m_ixRBegin+lBegin+n][m_ixRBegin+k] +=
				m_dDiffREdgeToNode[l][m]
				* dMult * (*m_pInterpNodeToREdge)[m][n]
				* m_dStateREdge[WIx][lBegin+m];
		}
		}

		if (a != nFiniteElements-1) {
			int lNext = lBegin + m_nVerticalOrder;
			for (int n = 0; n < m_nVerticalOrder; n++) {
				dDG[m_ixRBegin+lNext+n][m_ixRBegin+k] +=
					m_dDiffREdgeToNode[l][m_nVerticalOrder]
					* 0.5 * (*m_pInterpNodeToREdge)[0][n]
					* m_dStateREdge[WIx][lNext];
			}
		}
	}


	// Add the identity components
	for (int k = 0; k < m_nColumnStateSize; k++) {
		dDG[k][k] += 1.0 / m_dDeltaT;
	}
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::Evaluate(
	const double * dX,
	double * dF
) {
	// Prepare the column
	PrepareColumn(dX);

	// Evaluate the zero equations
	BuildF(dX, dF);
}

///////////////////////////////////////////////////////////////////////////////
// GLOBAL CONTEXT
///////////////////////////////////////////////////////////////////////////////

#ifdef USE_PETSC
PetscErrorCode VerticalDynamicsFEM_FormFunction(
	SNES snes,
	Vec x,
	Vec f,
	void * pDyn
) {
	// Pointers to the PetSc vector data
	const double * dX;
	double * dF;

	// Get pointers to the vector data
	VecGetArrayRead(x, &dX);
	VecGetArray(f, &dF);

	// Cast the context to VerticalDynamics and call Evaluate
	((VerticalDynamicsFEM*)(pDyn))->Evaluate(dX, dF);

	// Restore the array
	VecRestoreArrayRead(x, &dX);
	VecRestoreArray(f, &dF);

	return 0;
}
#endif
///////////////////////////////////////////////////////////////////////////////

