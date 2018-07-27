///////////////////////////////////////////////////////////////////////////////
///
///	\file    TempestInitialize.h
///	\author  Paul Ullrich
///	\version June 4, 2014
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

#ifndef _TEMPESTINITIALIZE_H_
#define _TEMPESTINITIALIZE_H_

#include "Defines.h"
#include "Model.h"
#include "TimestepSchemeStrang.h"
#include "TimestepSchemeARS222.h"
#include "TimestepSchemeARS232.h"
#include "TimestepSchemeARK232.h"
#include "TimestepSchemeGARK2.h"
#include "TimestepSchemeARS343.h"
#include "TimestepSchemeARS443.h"
#include "TimestepSchemeSSP3332.h"
#include "TimestepSchemeSplitExp.h"
#include "TimestepSchemeARKode.h"
#include "HorizontalDynamicsStub.h"
#include "HorizontalDynamicsFEM.h"
#include "VerticalDynamicsStub.h"
#include "VerticalDynamicsFEM.h"
#include "VerticalDynamicsSchur.h"
#include "VerticalDynamicsFLL.h"
#include "OutputManagerComposite.h"
#include "OutputManagerReference.h"
#include "OutputManagerChecksum.h"
#include "GridCSGLL.h"
#include "GridCartesianGLL.h"
#include "VerticalStretch.h"

#include "TimeObj.h"
#include "Announce.h"
#include "CommandLine.h"
#include "STLStringHelper.h"

#include <mpi.h>

#include <string>

#ifdef TEMPEST_PETSC
#include <petscsnes.h>
#endif

///////////////////////////////////////////////////////////////////////////////

struct _TempestCommandLineVariables {
	bool fNoOutput;
	std::string strOutputDir;
	std::string strOutputPrefix;
	std::string strRestartFile;
	int nOutputsPerFile;
	Time timeOutputDeltaT;
	Time timeOutputRestartDeltaT;
	Time timeDeltaT;
	Time timeEndTime;
	int nOutputResX;
	int nOutputResY;
	int nOutputResZ;
	bool fOutputVorticity;
	bool fOutputDivergence;
	bool fOutputTemperature;
	bool fOutputSurfacePressure;
	bool fOutputRichardson;
	bool fNoReferenceState;
	bool fNoTracers;
	bool fNoHyperviscosity;
	int nHyperviscosityOrder;
	double dNuScalar;
	double dNuDiv;
	double dNuVort;
	double dInstepNuDiv;
	bool fExplicitVertical;
	std::string strVerticalStaggering;
	bool fForceMassFluxOnLevels;
	std::string strVerticalStretch;
	std::string strVerticalDiscretization;
	int nVerticalHyperdiffOrder;
	std::string strTimestepScheme;
	std::string strHorizontalDynamics;
	std::string strVerticalDynamics;
        int nResolutionX;
	int nResolutionY;
	int nLevels;
	int nHorizontalOrder;
	int nVerticalOrder;
    	int iARKode_nvectors;
	double dARKode_rtol;
	double dARKode_atol;
	bool fARKode_DynamicStepSize;
	bool fARKode_aafp;       
	int iARKode_ErrController;
	int iARKode_AAFPAccelVec;
	int iARKode_NonlinIters;
	int iARKode_LinIters;
	int iARKode_Predictor;
	int iARKode_VAtol_vel;
	int iARKode_VAtol_rho;
	int iARKode_VAtol_theta;
	std::string strARKode_ButcherTable;
	std::string strARKode_StepOut;
        bool fARKode_Diagnostics;
    	bool fARKode_UsePreconditioning;
    	bool fARKode_ColumnSolver;
    	bool fFullyImplicit;
};

///////////////////////////////////////////////////////////////////////////////

#define _TempestDefineCommandLineDefault(TestCaseName) \
	CommandLineBool(_tempestvars.fNoOutput, "output_none"); \
	CommandLineString(_tempestvars.strOutputDir, "output_dir", "out" TestCaseName); \
	CommandLineString(_tempestvars.strOutputPrefix, "output_prefix", "out"); \
	CommandLineString(_tempestvars.strRestartFile, "restart_file", ""); \
	CommandLineInt(_tempestvars.nOutputsPerFile, "output_perfile", -1); \
	CommandLineDeltaTime(_tempestvars.timeOutputRestartDeltaT, "output_restart_dt", ""); \
	CommandLineInt(_tempestvars.nOutputResX, "output_x", 360); \
	CommandLineInt(_tempestvars.nOutputResY, "output_y", 180); \
	CommandLineInt(_tempestvars.nOutputResZ, "output_z", 0); \
	CommandLineBool(_tempestvars.fOutputVorticity, "output_vort"); \
	CommandLineBool(_tempestvars.fOutputDivergence, "output_div"); \
	CommandLineBool(_tempestvars.fOutputTemperature, "output_temp"); \
	CommandLineBool(_tempestvars.fOutputSurfacePressure, "output_ps"); \
	CommandLineBool(_tempestvars.fOutputRichardson, "output_Ri"); \
	CommandLineBool(_tempestvars.fNoReferenceState, "norefstate"); \
	CommandLineBool(_tempestvars.fNoTracers, "notracers"); \
	CommandLineBool(_tempestvars.fNoHyperviscosity, "nohypervis"); \
	CommandLineInt(_tempestvars.nHyperviscosityOrder, "hypervisorder", 4); \
	CommandLineDouble(_tempestvars.dNuScalar, "nu", 1.0e15); \
	CommandLineDouble(_tempestvars.dNuDiv, "nud", 1.0e15); \
	CommandLineDouble(_tempestvars.dNuVort, "nuv", 1.0e15); \
	CommandLineDouble(_tempestvars.dInstepNuDiv, "inud", 0.0); \
	CommandLineBool(_tempestvars.fExplicitVertical, "explicitvertical"); \
	CommandLineStringD(_tempestvars.strVerticalStaggering, "vstagger", "CPH", "(LEV | INT | LOR | CPH)"); \
	CommandLineStringD(_tempestvars.strVerticalDiscretization, "vdisc", "FE", "(FE | FV)"); \
	CommandLineBool(_tempestvars.fForceMassFluxOnLevels, "vmassfluxlevels"); \
	CommandLineString(_tempestvars.strVerticalStretch, "vstretch", "uniform"); \
	CommandLineInt(_tempestvars.nVerticalHyperdiffOrder, "vhypervisorder", 0); \
	CommandLineString(_tempestvars.strTimestepScheme, "timescheme", "strang"); \
	CommandLineStringD(_tempestvars.strVerticalDynamics, "vmethod", "DEFAULT", "(DEFAULT | SCHUR | FLL)"); \
	CommandLineInt(_tempestvars.iARKode_nvectors, "arkode_nvectors", 50); \
	CommandLineDouble(_tempestvars.dARKode_rtol, "arkode_rtol", 1.0e-6); \
	CommandLineDouble(_tempestvars.dARKode_atol, "arkode_atol", 1.0e-11); \
	CommandLineBool(_tempestvars.fARKode_DynamicStepSize, "arkode_dynamicstepsize"); \
	CommandLineBool(_tempestvars.fARKode_aafp, "arkode_aafp"); \
	CommandLineInt(_tempestvars.iARKode_ErrController, "arkode_errcontrol", 0); \
	CommandLineInt(_tempestvars.iARKode_AAFPAccelVec, "arkode_aafpaccelvec", 0); \
	CommandLineInt(_tempestvars.iARKode_NonlinIters, "arkode_nonliniters", 0); \
	CommandLineInt(_tempestvars.iARKode_LinIters, "arkode_liniters", 0); \
	CommandLineInt(_tempestvars.iARKode_Predictor, "arkode_predictor", 0); \
	CommandLineInt(_tempestvars.iARKode_VAtol_vel, "arkode_vatol_vel", 1.0e-3); \
	CommandLineInt(_tempestvars.iARKode_VAtol_rho, "arkode_vatol_rho", 1.0e-4); \
	CommandLineInt(_tempestvars.iARKode_VAtol_theta, "arkode_vatol_theta", 1.0e-2); \
	CommandLineString(_tempestvars.strARKode_ButcherTable, "arkode_butchertable", ""); \
	CommandLineString(_tempestvars.strARKode_StepOut, "timestep_out", "timestep_profile.out"); \
	CommandLineBool(_tempestvars.fARKode_Diagnostics, "arkode_diagnostics"); \
	CommandLineBool(_tempestvars.fARKode_UsePreconditioning, "arkode_usepreconditioning"); \
	CommandLineBool(_tempestvars.fARKode_ColumnSolver, "arkode_columnsolver"); \
	CommandLineBool(_tempestvars.fFullyImplicit, "fullyimplicit");

///////////////////////////////////////////////////////////////////////////////

#define BeginTempestCommandLine(TestCaseName) \
	_TempestCommandLineVariables _tempestvars; \
	BeginCommandLine() \
	_TempestDefineCommandLineDefault(TestCaseName)

#define EndTempestCommandLine(argv) \
	EndCommandLine(argv)

#define SetDefaultResolution(default_resolution) \
	CommandLineInt(_tempestvars.nResolutionX, "resolution", default_resolution);

#define SetDefaultResolutionX(default_resolution_x) \
	CommandLineInt(_tempestvars.nResolutionX, "resx", default_resolution_x);

#define SetDefaultResolutionY(default_resolution_y) \
	CommandLineInt(_tempestvars.nResolutionY, "resy", default_resolution_y);

#define SetDefaultLevels(default_levels) \
	CommandLineInt(_tempestvars.nLevels, "levels", default_levels);

#define SetDefaultOutputDeltaT(default_output_dt) \
	CommandLineDeltaTime(_tempestvars.timeOutputDeltaT, "outputtime", default_output_dt);

#define SetDefaultDeltaT(default_dt) \
	CommandLineDeltaTime(_tempestvars.timeDeltaT, "dt", default_dt);

#define SetDefaultEndTime(default_endtime) \
	CommandLineFixedTime(_tempestvars.timeEndTime, "endtime", default_endtime);

#define SetDefaultHorizontalOrder(default_order) \
	CommandLineInt(_tempestvars.nHorizontalOrder, "order", default_order)

#define SetDefaultVerticalOrder(default_vertorder) \
	CommandLineInt(_tempestvars.nVerticalOrder, "vertorder", default_vertorder)

///////////////////////////////////////////////////////////////////////////////

void _TempestSetupMethodOfLines(
	Model & model,
	_TempestCommandLineVariables & vars
) {
	// Set the timestep scheme
	AnnounceStartBlock("Initializing time scheme");

	STLStringHelper::ToLower(vars.strTimestepScheme);
	if (vars.strTimestepScheme == "strang") {
		model.SetTimestepScheme(
			new TimestepSchemeStrang(model));

	} else if (vars.strTimestepScheme == "strang/fe") {
		model.SetTimestepScheme(
			new TimestepSchemeStrang(
				model, 0.0, TimestepSchemeStrang::ForwardEuler));

	} else if (vars.strTimestepScheme == "strang/rk4") {
		model.SetTimestepScheme(
			new TimestepSchemeStrang(
				model, 0.0, TimestepSchemeStrang::RungeKutta4));

	} else if (vars.strTimestepScheme == "strang/rk3") {
		model.SetTimestepScheme(
			new TimestepSchemeStrang(
				model, 0.0, TimestepSchemeStrang::RungeKuttaSSP3));

	} else if (vars.strTimestepScheme == "strang/kgu35") {
		model.SetTimestepScheme(
			new TimestepSchemeStrang(
				model, 0.0, TimestepSchemeStrang::KinnmarkGrayUllrich35));

	} else if (vars.strTimestepScheme == "strang/ssprk53") {
		model.SetTimestepScheme(
			new TimestepSchemeStrang(
				model, 0.0, TimestepSchemeStrang::RungeKuttaSSPRK53));

	} else if (vars.strTimestepScheme == "ars222") {
		model.SetTimestepScheme(
			new TimestepSchemeARS222(model));

	} else if (vars.strTimestepScheme == "ars232") {
		model.SetTimestepScheme(
			new TimestepSchemeARS232(model));

	} else if (vars.strTimestepScheme == "ark232") {
		model.SetTimestepScheme(
			new TimestepSchemeARK232(model));

	} else if (vars.strTimestepScheme == "gark2") {
		model.SetTimestepScheme(
			new TimestepSchemeGARK2(model));

	} else if (vars.strTimestepScheme == "ars343") {
		model.SetTimestepScheme(
			new TimestepSchemeARS343(model));

	} else if (vars.strTimestepScheme == "ars443") {
		model.SetTimestepScheme(
			new TimestepSchemeARS443(model));

	} else if (vars.strTimestepScheme == "ssp3_332") {
		model.SetTimestepScheme(
			new TimestepSchemeSSP3332(model));

	} else if (vars.strTimestepScheme == "spex") {
		model.SetTimestepScheme(
			new TimestepSchemeSplitExp(model));

#ifdef USE_SUNDIALS
	} else if (vars.strTimestepScheme == "arkode") {
 	        ARKodeCommandLineVariables ARKodeVars;

		STLStringHelper::ToLower(vars.strARKode_ButcherTable);

		if (vars.fFullyImplicit) {
		  vars.fExplicitVertical   = true;
		  ARKodeVars.FullyExplicit = false;
		  ARKodeVars.FullyImplicit = vars.fFullyImplicit;
		} else {
		  ARKodeVars.FullyExplicit = vars.fExplicitVertical;
		  ARKodeVars.FullyImplicit = false;
		}
		
		ARKodeVars.nvectors        = vars.iARKode_nvectors;
		ARKodeVars.rtol            = vars.dARKode_rtol;
		ARKodeVars.atol            = vars.dARKode_atol;
		ARKodeVars.DynamicStepSize = vars.fARKode_DynamicStepSize;
		ARKodeVars.ErrController   = vars.iARKode_ErrController;
                ARKodeVars.AAFP            = vars.fARKode_aafp;
		ARKodeVars.AAFPAccelVec    = vars.iARKode_AAFPAccelVec;
		ARKodeVars.NonlinIters     = vars.iARKode_NonlinIters;
		ARKodeVars.LinIters        = vars.iARKode_LinIters;
		ARKodeVars.Predictor       = vars.iARKode_Predictor;
		ARKodeVars.VAtol_vel	   = vars.iARKode_VAtol_vel;
		ARKodeVars.VAtol_rho	   = vars.iARKode_VAtol_rho;
		ARKodeVars.VAtol_theta	   = vars.iARKode_VAtol_theta;
		ARKodeVars.StepOut         = vars.strARKode_StepOut;
		ARKodeVars.ButcherTable       = vars.strARKode_ButcherTable;
		ARKodeVars.WriteDiagnostics   = vars.fARKode_Diagnostics;
		ARKodeVars.UsePreconditioning = vars.fARKode_UsePreconditioning,
		ARKodeVars.ColumnSolver       = vars.fARKode_ColumnSolver,

		model.SetTimestepScheme(
			new TimestepSchemeARKode(model, ARKodeVars));
#endif

	} else {
		_EXCEPTIONT("Invalid timescheme: Expected "
			"\"Strang\", \"ARS222\", \"ARS232\", \"ARK232\", "
			"\"ARS343\", \"ARS443\", \"SSP3_332\"");
	}
	AnnounceEndBlock("Done");

	// Set the horizontal dynamics
	AnnounceStartBlock("Initializing horizontal dynamics");

	if (vars.fNoHyperviscosity) {
		vars.dNuScalar = 0.0;
		vars.dNuDiv = 0.0;
		vars.dNuVort = 0.0;
	}

	model.SetHorizontalDynamics(
		new HorizontalDynamicsFEM(
			model,
			vars.nHorizontalOrder,
			vars.nHyperviscosityOrder,
			vars.dNuScalar,
			vars.dNuDiv,
			vars.dNuVort,
			vars.dInstepNuDiv));

	AnnounceEndBlock("Done");

	// Vertical staggering
	STLStringHelper::ToLower(vars.strVerticalStaggering);

	// Vertical method
	STLStringHelper::ToLower(vars.strVerticalDynamics);

	// Set the vertical dynamics
	AnnounceStartBlock("Initializing vertical dynamics");
	if (vars.nLevels == 1) {
		model.SetVerticalDynamics(
			new VerticalDynamicsStub(model));

	} else if (vars.strVerticalDynamics == "default") {
		model.SetVerticalDynamics(
			new VerticalDynamicsFEM(
				model,
				vars.nHorizontalOrder,
				vars.nVerticalOrder,
				vars.nVerticalHyperdiffOrder,
				vars.fExplicitVertical,
				!vars.fNoReferenceState,
				vars.fForceMassFluxOnLevels));

	} else if (vars.strVerticalDynamics == "schur") {
		model.SetVerticalDynamics(
			new VerticalDynamicsSchur(
				model,
				vars.nHorizontalOrder,
				vars.nVerticalOrder,
				vars.nVerticalHyperdiffOrder,
				vars.fExplicitVertical,
				!vars.fNoReferenceState,
				vars.fForceMassFluxOnLevels));

	} else if (vars.strVerticalDynamics == "fll") {
		model.SetVerticalDynamics(
			new VerticalDynamicsFLL(
				model,
				vars.nHorizontalOrder,
				vars.nVerticalOrder,
				vars.nVerticalHyperdiffOrder,
				vars.fExplicitVertical,
				!vars.fNoReferenceState,
				vars.fForceMassFluxOnLevels));

	} else {
		_EXCEPTIONT("Invalid --vmethod");
	}

	AnnounceEndBlock("Done");
}

///////////////////////////////////////////////////////////////////////////////

void _TempestSetupOutputManagers(
	Model & model,
	_TempestCommandLineVariables & vars
) {
	// Set the reference output manager for the model
	if (!vars.fNoOutput) {
		AnnounceStartBlock("Creating reference output manager");
		OutputManagerReference * pOutmanRef =
			new OutputManagerReference(
				*(model.GetGrid()),
				vars.timeOutputDeltaT,
				vars.strOutputDir,
				vars.strOutputPrefix,
				vars.nOutputsPerFile,
				vars.nOutputResX,
				vars.nOutputResY,
				vars.nOutputResZ,
				false,
				false);

		if (vars.fOutputVorticity) {
			pOutmanRef->OutputVorticity();
		}
		if (vars.fOutputDivergence) {
			pOutmanRef->OutputDivergence();
		}
		if (vars.fOutputTemperature) {
			pOutmanRef->OutputTemperature();
		}
		if (vars.fOutputSurfacePressure) {
			pOutmanRef->OutputSurfacePressure();
		}
		if (vars.fOutputRichardson) {
			pOutmanRef->OutputRichardson();
		}

		model.AttachOutputManager(pOutmanRef);
		AnnounceEndBlock("Done");
	}

	// Set the composite output manager for the model
	if (!vars.timeOutputRestartDeltaT.IsZero()) {
		AnnounceStartBlock("Creating composite output manager");
		model.AttachOutputManager(
			new OutputManagerComposite(
				*(model.GetGrid()),
				vars.timeOutputRestartDeltaT,
				vars.strOutputDir,
				vars.strOutputPrefix));
		AnnounceEndBlock("Done");
	}

	// Set the checksum output manager for the model
	AnnounceStartBlock("Creating checksum output manager");
	model.AttachOutputManager(
		new OutputManagerChecksum(
			*(model.GetGrid()),
			vars.timeOutputDeltaT));
	AnnounceEndBlock("Done");
}

///////////////////////////////////////////////////////////////////////////////

void _TempestSetupCubedSphereModel(
	Model & model,
	_TempestCommandLineVariables & vars
) {
	// Set the time step size and end time
	model.SetDeltaT(vars.timeDeltaT);
	model.SetEndTime(vars.timeEndTime);

	// Setup Method of Lines
	_TempestSetupMethodOfLines(model, vars);

	// Get the vertical discretization from --vdisc
	Grid::VerticalDiscretization eVerticalDiscretization;
	STLStringHelper::ToLower(vars.strVerticalDiscretization);
	if (vars.strVerticalDiscretization == "fe") {
		eVerticalDiscretization = Grid::VerticalDiscretization_FiniteElement;

	} else if (vars.strVerticalDiscretization == "fv") {
		eVerticalDiscretization = Grid::VerticalDiscretization_FiniteVolume;

	} else {
		_EXCEPTIONT("Invalid value for --vdisc");
	}

	// Get the vertical staggering from --vstagger
	Grid::VerticalStaggering eVerticalStaggering;
	STLStringHelper::ToLower(vars.strVerticalStaggering);
	if (vars.strVerticalStaggering == "lev") {
		eVerticalStaggering = Grid::VerticalStaggering_Levels;

	} else if (vars.strVerticalStaggering == "int") {
		eVerticalStaggering = Grid::VerticalStaggering_Interfaces;

	} else if (vars.strVerticalStaggering == "lor") {
		eVerticalStaggering = Grid::VerticalStaggering_Lorenz;

	} else if (vars.strVerticalStaggering == "cph") {
		eVerticalStaggering = Grid::VerticalStaggering_CharneyPhillips;

	} else {
		_EXCEPTIONT("Invalid value for --vstagger");
	}

	// Construct the Grid
	if (vars.strRestartFile == "") {
		AnnounceStartBlock("Constructing grid");
	
		// Maximum number of patches currently equals communicator size
		int nCommSize;
		MPI_Comm_size(MPI_COMM_WORLD, &nCommSize);

		if (nCommSize < 6) {
			nCommSize = 6;
		}

		GridCSGLL * pGrid = new GridCSGLL(model);

		pGrid->DefineParameters();

		pGrid->SetParameters(
			vars.nLevels,
			nCommSize,
			vars.nResolutionX,
			4,
			vars.nHorizontalOrder,
			vars.nVerticalOrder,
			eVerticalDiscretization,
			eVerticalStaggering);

		pGrid->InitializeDataLocal();

		// Set the vertical stretching function
		STLStringHelper::ToLower(vars.strVerticalStretch);
		if (vars.strVerticalStretch == "uniform") {

		} else if (vars.strVerticalStretch == "cubic") {
			pGrid->SetVerticalStretchFunction(
				new VerticalStretchCubic);

		} else if (vars.strVerticalStretch == "pwlinear") {
			pGrid->SetVerticalStretchFunction(
				new VerticalStretchPiecewiseLinear);

		} else {
			_EXCEPTIONT("Invalid value for --vstretch");
		}

		// Vertical type
		STLStringHelper::ToLower(vars.strVerticalDiscretization);

		// Set the Model Grid
		model.SetGrid(pGrid);

	// Set the Grid from Restart file
	} else {
		AnnounceStartBlock("Constructing Grid from restart file");

		Grid * pGrid = new GridCSGLL(model);

		pGrid->DefineParameters();

		model.SetGridFromRestartFile(pGrid, vars.strRestartFile);
	}

	AnnounceEndBlock("Done");

	// Setup OutputManagers
	_TempestSetupOutputManagers(model, vars);
}

///////////////////////////////////////////////////////////////////////////////

void _TempestSetupCartesianModel(
	Model & model,
	double dGDim[],
	double dRefLat,
	int iLatBC[],
	bool fCartesianXZ,
	_TempestCommandLineVariables & vars
) {
	// Set the time step size and end time
	model.SetDeltaT(vars.timeDeltaT);
	model.SetEndTime(vars.timeEndTime);

	// Setup Method of Lines
	_TempestSetupMethodOfLines(model, vars);

	// Get the vertical discretization from --vdisc
	Grid::VerticalDiscretization eVerticalDiscretization;
	STLStringHelper::ToLower(vars.strVerticalDiscretization);
	if (vars.strVerticalDiscretization == "fe") {
		eVerticalDiscretization = Grid::VerticalDiscretization_FiniteElement;

	} else if (vars.strVerticalDiscretization == "fv") {
		eVerticalDiscretization = Grid::VerticalDiscretization_FiniteVolume;

	} else {
		_EXCEPTIONT("Invalid value for --vdisc");
	}

	// Get the vertical staggering from --vstagger
	Grid::VerticalStaggering eVerticalStaggering;
	STLStringHelper::ToLower(vars.strVerticalStaggering);
	if (vars.strVerticalStaggering == "lev") {
		eVerticalStaggering = Grid::VerticalStaggering_Levels;

	} else if (vars.strVerticalStaggering == "int") {
		eVerticalStaggering = Grid::VerticalStaggering_Interfaces;

	} else if (vars.strVerticalStaggering == "lor") {
		eVerticalStaggering = Grid::VerticalStaggering_Lorenz;

	} else if (vars.strVerticalStaggering == "cph") {
		eVerticalStaggering = Grid::VerticalStaggering_CharneyPhillips;

	} else {
		_EXCEPTIONT("Invalid value for --vstagger");
	}

	// Construct the Grid
	if (vars.strRestartFile == "") {
		AnnounceStartBlock("Constructing grid");

		// Maximum number of patches currently equals communicator size
		int nCommSize;
		MPI_Comm_size(MPI_COMM_WORLD, &nCommSize);

		GridCartesianGLL * pGrid = new GridCartesianGLL(model);

		pGrid->DefineParameters();

		pGrid->SetParameters(
			vars.nLevels,
			nCommSize,
			vars.nResolutionX,
			vars.nResolutionY,
			4,
			vars.nHorizontalOrder,
			vars.nVerticalOrder,
			dGDim,
			dRefLat,
			iLatBC,
			fCartesianXZ,
			eVerticalDiscretization,
			eVerticalStaggering);

		pGrid->InitializeDataLocal();

		// Set the vertical stretching function
		STLStringHelper::ToLower(vars.strHorizontalDynamics);
		if (vars.strVerticalStretch == "uniform") {

		} else if (vars.strVerticalStretch == "cubic") {
			pGrid->SetVerticalStretchFunction(
				new VerticalStretchCubic);

		} else if (vars.strVerticalStretch == "pwlinear") {
			pGrid->SetVerticalStretchFunction(
				new VerticalStretchPiecewiseLinear);

		} else {
			_EXCEPTIONT("Invalid value for --vstretch");

		}

		// Set the Model Grid
		model.SetGrid(pGrid);

	// Set the Grid from Restart file
	} else {
		AnnounceStartBlock("Constructing Grid from restart file");

		Grid * pGrid = new GridCartesianGLL(model);

		pGrid->DefineParameters();

		model.SetGridFromRestartFile(pGrid, vars.strRestartFile);
	}

	AnnounceEndBlock("Done");

	// Setup OutputManagers
	_TempestSetupOutputManagers(model, vars);
}

///////////////////////////////////////////////////////////////////////////////

#define TempestSetupCubedSphereModel(model) \
	_TempestSetupCubedSphereModel(model, _tempestvars);

#define TempestSetupCartesianModel(model, dimensions, latitude, lateralBC, CartesianXZ) \
	_TempestSetupCartesianModel(model, dimensions, latitude, lateralBC, CartesianXZ, _tempestvars);

///////////////////////////////////////////////////////////////////////////////

void TempestInitialize(int * argc, char*** argv) {

#ifdef TEMPEST_PETSC
	// Initialize PetSc
	PetscInitialize(argc, argv, NULL, NULL);
#endif
#ifdef TEMPEST_MPIOMP
	// Initialize MPI
	MPI_Init(argc, argv);
#endif

}

///////////////////////////////////////////////////////////////////////////////

void TempestAbort() {

#ifdef TEMPEST_MPIOMP
	// Abort
	MPI_Abort(MPI_COMM_WORLD, 1);
#endif

}

///////////////////////////////////////////////////////////////////////////////

void TempestDeinitialize() {

#ifdef TEMPEST_PETSC
	// Finalize PetSc
	PetscFinalize();
#endif
#ifdef TEMPEST_MPIOMP
	// Finalize MPI
	MPI_Finalize();
#endif

}

///////////////////////////////////////////////////////////////////////////////

#endif

