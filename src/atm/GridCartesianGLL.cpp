///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridCartesianGLL.cpp
///	\author  Paul Ullrich
///	\version September 19, 2013
///
///	<remarks>
///		Copyright 2000-2010 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "GridCartesianGLL.h"
#include "GridPatchCartesianGLL.h"
#include "Model.h"
#include "TestCase.h"
#include "GridSpacing.h"

#include "Direction.h"
#include "CubedSphereTrans.h"
#include "PolynomialInterp.h"
#include "GaussQuadrature.h"
#include "GaussLobattoQuadrature.h"

#include "Announce.h"
#include "MathHelper.h"

#include <cmath>

///////////////////////////////////////////////////////////////////////////////

GridCartesianGLL::GridCartesianGLL(
	const Model & model,
	int nBaseResolution,
	int nRefinementRatio,
	int nHorizontalOrder,
	int nVerticalOrder,
	int nRElements
) :
	// Call up the stack
	GridGLL::GridGLL(
		model,
		nBaseResolution,
		nRefinementRatio,
		nHorizontalOrder,
		nVerticalOrder,
		nRElements)
{
	// Set the reference length scale
	m_dReferenceLength = 0.5 * M_PI / 30.0;
}

///////////////////////////////////////////////////////////////////////////////

void GridCartesianGLL::Initialize() {

	// Call up the stack
	Grid::Initialize();

	// Initialize the vertical coordinate
	double dDeltaElement =
		static_cast<double>(m_nVerticalOrder)
		/ static_cast<double>(m_nRElements);

	InitializeVerticalCoordinate(
		GridSpacingMixedGaussLobatto(dDeltaElement, 0.0, m_nVerticalOrder)
	);

	// Determine number of usable processors
	int nCommSize;
	MPI_Comm_size(MPI_COMM_WORLD, &nCommSize);

	int nProcsPerDirection = Max((int)ISqrt(nCommSize / 6), 1);

	int nProcsPerPanel = nProcsPerDirection * nProcsPerDirection;

	int nDistributedPatches = 6 * nProcsPerPanel;

	if (nDistributedPatches < nCommSize) {
		Announce("WARNING: Patch / thread mismatch: "
			"%i threads will be unutilized",
			nCommSize - nDistributedPatches);
	}

	// Determine arrangement of elements on processors
	if (GetABaseResolution() % nProcsPerDirection != 0) {
		_EXCEPTIONT("\n(UNIMPLEMENTED) Currently elements must be "
			"equally divided among processors.");
	}

	int nElementsPerDirection = GetABaseResolution() / nProcsPerDirection;
	DataVector<int> iBoxBegin;
	iBoxBegin.Initialize(nProcsPerDirection + 1);

	iBoxBegin[0] = 0;
	for (int n = 1; n < nProcsPerDirection; n++) {
		iBoxBegin[n] = n * nElementsPerDirection;
	}
	iBoxBegin[nProcsPerDirection] = GetABaseResolution();

	// Create master patch for each panel
	std::vector<GridPatch *> pPatches;
	pPatches.reserve(nDistributedPatches);

	for (int n = 0; n < 6; n++) {
	for (int i = 0; i < nProcsPerDirection; i++) {
	for (int j = 0; j < nProcsPerDirection; j++) {
		double dDeltaA = 0.5 * M_PI / GetABaseResolution();

		GridSpacingGaussLobattoRepeated
			glspacing(dDeltaA, -0.25 * M_PI, m_nHorizontalOrder);

		PatchBox boxMaster(
			n, 0, m_model.GetHaloElements(),
			m_nHorizontalOrder * iBoxBegin[i],
			m_nHorizontalOrder * iBoxBegin[i+1],
			m_nHorizontalOrder * iBoxBegin[j],
			m_nHorizontalOrder * iBoxBegin[j+1],
			glspacing,
			glspacing);

		int ixPatch = n * nProcsPerPanel + i * nProcsPerDirection + j;

		pPatches.push_back(
			AddPatch(
				new GridPatchCartesianGLL(
					(*this),
					ixPatch,
					boxMaster,
					m_nHorizontalOrder,
					m_nVerticalOrder)));
	}
	}
	}

	if (pPatches.size() != nDistributedPatches) {
		_EXCEPTIONT("Logic error");
	}

	// Set up connectivity
	for (int n = 0; n < nDistributedPatches; n++) {
		int iDestPanel;
		int iDestI;
		int iDestJ;
		bool fSwitchAB;
		bool fSwitchPar;
		bool fSwitchPerp;

		int iSrcPanel = n / nProcsPerPanel;
		int iSrcI = (n % nProcsPerPanel) / nProcsPerDirection;
		int iSrcJ = (n % nProcsPerPanel) % nProcsPerDirection;

		int iDestN;

		// Loop through all directions
		for (int iDir = 0; iDir < 8; iDir++) {
			Direction dir = (Direction)(iDir);

			// Find the panel index in this direction
			int iSrcInew = iSrcI;
			int iSrcJnew = iSrcJ;

			DirectionIncrement(dir, iSrcInew, iSrcJnew);

			CubedSphereTrans::RelativeCoord(
				nProcsPerDirection,
				iSrcPanel, iSrcInew, iSrcJnew,
				iDestPanel, iDestI, iDestJ,
				fSwitchAB, fSwitchPar, fSwitchPerp);

			if (iDestPanel == InvalidPanel) {
				continue;
			}

			iDestN =
				iDestPanel * nProcsPerPanel
				+ iDestI * nProcsPerDirection
				+ iDestJ;

			// Find the opposing direction to return to this panel
			Direction dirOpposing =
				CubedSphereTrans::OpposingDirection(iSrcPanel, iDestPanel, dir);

			// Set of the exterior connection
			Connectivity::ExteriorConnect(
				pPatches[n],
				dir,
				pPatches[iDestN]);
		}
	}

	// Distribute patches to processors
	Grid::DistributePatches();

	// Interpolation coefficients from nodes to interfaces and vice versa
	m_dInterpNodeToREdge.Initialize(m_nVerticalOrder+1, m_nVerticalOrder);
	m_dInterpREdgeToNode.Initialize(m_nVerticalOrder, m_nVerticalOrder+1);

	for (int n = 0; n < m_nVerticalOrder+1; n++) {
		PolynomialInterp::LagrangianPolynomialCoeffs(
			m_nVerticalOrder,
			m_dREtaLevels,
			m_dInterpNodeToREdge[n],
			m_dREtaInterfaces[n]);
	}
	for (int n = 0; n < m_nVerticalOrder; n++) {
		PolynomialInterp::LagrangianPolynomialCoeffs(
			m_nVerticalOrder+1,
			m_dREtaInterfaces,
			m_dInterpREdgeToNode[n],
			m_dREtaLevels[n]);
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridCartesianGLL::InitializeVerticalCoordinate(
	const GridSpacing & aGridSpacing
) {
	// Call to Grid
	Grid::InitializeVerticalCoordinate(aGridSpacing);
}

///////////////////////////////////////////////////////////////////////////////

void GridCartesianGLL::ConvertReferenceToABP(
	const DataVector<double> & dXReference,
	const DataVector<double> & dYReference,
	DataVector<double> & dAlpha,
	DataVector<double> & dBeta,
	DataVector<int> & iPanel
) const {
	if (dXReference.GetRows() != dYReference.GetRows()) {
		_EXCEPTIONT("XReference and YReference must have same length.");
	}

	// Resize arrays
	dAlpha.Initialize(dXReference.GetRows());
	dBeta .Initialize(dXReference.GetRows());
	iPanel.Initialize(dXReference.GetRows());

	// Loop over all coordinates
	for (int i = 0; i < dXReference.GetRows(); i++) {
		CubedSphereTrans::ABPFromRLL(
			dXReference[i],
			dYReference[i],
			dAlpha[i],
			dBeta[i],
			iPanel[i]);
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridCartesianGLL::ApplyDSS(
	int iDataUpdate,
	DataType eDataType
) {
	// Exchange data between nodes
	Exchange(eDataType, iDataUpdate);

	// Post-process velocities across panel edges and
	// perform direct stiffness summation (DSS)
	for (int n = 0; n < GetActivePatchCount(); n++) {
		GridPatchCartesianGLL * pPatch =
			dynamic_cast<GridPatchCartesianGLL*>(GetActivePatch(n));

		const PatchBox & box = pPatch->GetPatchBox();

		// Patch-specific quantities
		int nAInteriorWidth = pPatch->GetPatchBox().GetAInteriorWidth();
		int nBInteriorWidth = pPatch->GetPatchBox().GetBInteriorWidth();

		int nAElements = nAInteriorWidth / m_nHorizontalOrder;
		int nBElements = nBInteriorWidth / m_nHorizontalOrder;

		if ((nAInteriorWidth % m_nHorizontalOrder) != 0) {
			_EXCEPTIONT("Logic Error: Invalid PatchBox alpha spacing");
		}
		if ((nBInteriorWidth % m_nHorizontalOrder) != 0) {
			_EXCEPTIONT("Logic Error: Invalid PatchBox beta spacing");
		}

		// Apply panel transforms to velocity data
		if (eDataType == DataType_State) {
			pPatch->TransformHaloVelocities(iDataUpdate);
		}

		// Panels in each coordinate direction
		int ixRightPanel =
			pPatch->GetNeighborPanel(Direction_Right);
		int ixTopPanel =
			pPatch->GetNeighborPanel(Direction_Top);
		int ixLeftPanel =
			pPatch->GetNeighborPanel(Direction_Left);
		int ixBottomPanel =
			pPatch->GetNeighborPanel(Direction_Bottom);

		int ixTopRightPanel =
			pPatch->GetNeighborPanel(Direction_TopRight);
		int ixTopLeftPanel =
			pPatch->GetNeighborPanel(Direction_TopLeft);
		int ixBottomLeftPanel =
			pPatch->GetNeighborPanel(Direction_BottomLeft);
		int ixBottomRightPanel =
			pPatch->GetNeighborPanel(Direction_BottomRight);

		// Loop through all components associated with this DataType
		int nComponents;
		if (eDataType == DataType_State) {
			nComponents = m_model.GetEquationSet().GetComponents();
		} else if (eDataType == DataType_Tracers) {
			nComponents = m_model.GetEquationSet().GetTracers();
		} else if (eDataType == DataType_Vorticity) {
			nComponents = 1;
		} else if (eDataType == DataType_Divergence) {
			nComponents = 1;
		} else {
			_EXCEPTIONT("Invalid DataType");
		}

		// Perform Direct Stiffness Summation (DSS)
		for (int c = 0; c < nComponents; c++) {

			// Obtain the array of working data
			int nRElements = GetRElements();
			double *** pDataUpdate;
			if (eDataType == DataType_State) {
				pDataUpdate =
					pPatch->GetDataState(iDataUpdate, GetVarLocation(c))[c];

				if (GetVarLocation(c) == DataLocation_REdge) {
					nRElements++;
				}

			} else if (eDataType == DataType_Tracers) {
				pDataUpdate =
					pPatch->GetDataTracers(iDataUpdate)[c];
			} else if (eDataType == DataType_Vorticity) {
				pDataUpdate = pPatch->GetDataVorticity();
			} else if (eDataType == DataType_Divergence) {
				pDataUpdate = pPatch->GetDataDivergence();
			}

			for (int k = 0; k < nRElements; k++) {

				// Average in the alpha direction
				for (int a = 0; a <= nAElements; a++) {
					int iA = a * m_nHorizontalOrder + box.GetHaloElements();

					// Do not average across cubed-sphere corners
					int jBegin = box.GetBInteriorBegin()-1;
					int jEnd = box.GetBInteriorEnd()+1;

					if (((a == 0) &&
							(ixTopLeftPanel == InvalidPanel)) ||
						((a == nAElements) &&
							(ixTopRightPanel == InvalidPanel))
					) {
						jEnd -= 2;
					}
					if (((a == 0) &&
							(ixBottomLeftPanel == InvalidPanel)) ||
						((a == nAElements) &&
							(ixBottomRightPanel == InvalidPanel))
					) {
						jBegin += 2;
					}

					// Perform averaging across edge
					for (int j = jBegin; j < jEnd; j++) {
						pDataUpdate[k][iA][j] = 0.5 * (
							+ pDataUpdate[k][iA  ][j]
							+ pDataUpdate[k][iA-1][j]);

						pDataUpdate[k][iA-1][j] = pDataUpdate[k][iA][j];
					}
				}

				// Average in the beta direction
				for (int b = 0; b <= nBElements; b++) {
					int iB = b * m_nHorizontalOrder + box.GetHaloElements();

					// Do not average across cubed-sphere corners
					int iBegin = box.GetAInteriorBegin()-1;
					int iEnd = box.GetAInteriorEnd()+1;

					if (((b == 0) &&
							(ixBottomLeftPanel == InvalidPanel)) ||
						((b == nAElements) &&
							(ixTopLeftPanel == InvalidPanel))
					) {
						iBegin += 2;
					}
					if (((b == 0) &&
							(ixBottomRightPanel == InvalidPanel)) ||
						((b == nAElements) &&
							(ixTopRightPanel == InvalidPanel))
					) {
						iEnd -= 2;
					}

					for (int i = iBegin; i < iEnd; i++) {
						pDataUpdate[k][i][iB] = 0.5 * (
							+ pDataUpdate[k][i][iB  ]
							+ pDataUpdate[k][i][iB-1]);

						pDataUpdate[k][i][iB-1] = pDataUpdate[k][i][iB];
					}
				}

				// Average at cubed-sphere corners (nodes of connectivity 3)
				if (ixTopRightPanel == InvalidPanel) {
					int iA = box.GetAInteriorEnd()-1;
					int iB = box.GetBInteriorEnd()-1;

					pDataUpdate[k][iA][iB] = (1.0/3.0) * (
						+ pDataUpdate[k][iA  ][iB  ]
						+ pDataUpdate[k][iA+1][iB  ]
						+ pDataUpdate[k][iA  ][iB+1]);
				}

				if (ixTopLeftPanel == InvalidPanel) {
					int iA = box.GetAInteriorBegin();
					int iB = box.GetBInteriorEnd()-1;

					pDataUpdate[k][iA][iB] = (1.0/3.0) * (
						+ pDataUpdate[k][iA  ][iB  ]
						+ pDataUpdate[k][iA-1][iB  ]
						+ pDataUpdate[k][iA  ][iB+1]);
				}

				if (ixBottomLeftPanel == InvalidPanel) {
					int iA = box.GetAInteriorBegin();
					int iB = box.GetBInteriorBegin();

					pDataUpdate[k][iA][iB] = (1.0/3.0) * (
						+ pDataUpdate[k][iA  ][iB  ]
						+ pDataUpdate[k][iA-1][iB  ]
						+ pDataUpdate[k][iA  ][iB-1]);
				}

				if (ixBottomRightPanel == InvalidPanel) {
					int iA = box.GetAInteriorEnd()-1;
					int iB = box.GetBInteriorBegin();

					pDataUpdate[k][iA][iB] = (1.0/3.0) * (
						+ pDataUpdate[k][iA  ][iB  ]
						+ pDataUpdate[k][iA+1][iB  ]
						+ pDataUpdate[k][iA  ][iB-1]);
				}
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////