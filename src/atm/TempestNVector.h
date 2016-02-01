///////////////////////////////////////////////////////////////////////////////
///
///	\file    TempestNVector.h
///	\author  Daniel R. Reynolds
///	\version February 1, 2016
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

#ifndef _TEMPESTNVECTOR_H_
#define _TEMPESTNVECTOR_H_

///////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////
/// This is the header file for an implementation of an NVECTOR
/// package, specifically suited to interface with the TEMPEST "Grid" 
/// data structure.
///
/// Part I of this file contains global variables required to handle the static
/// NVector registry.
///
/// Part II of this file contains the NVector data structure declaration.
///
/// Part III of this file contains the prototype for new NVector initialization.
///
/// Part IV of this file contains prototypes for the required vector operations
/// which operate on the NVector.
///
/// NOTES:
///
/// The generic N_Vector structure is defined in the SUNDIALS header file nvector.h
///
/// The definition of the type "realtype" is in the SUNDIALS header file 
/// sundialstypes.h and may be changed according to the user's needs during 
/// SUNDIALS installation -- this value should match the datatype used for 
/// TEMPEST state variables (double).
///
///////////////////////////////////////////////////////////////////////////////

#ifdef __cplusplus     // flatten namespace for use by both C and C++
extern "C" {
#endif

#include <sundials/sundials_nvector.h>

///// Part I -- global variables for NVector registry /////

#define MAX_TEMPEST_NVECTORS 100


///// Part II -- NVector data structure /////

struct _N_VectorContent_Tempest {
  int mVectorIndex = -1;      // index into registry
  Grid * mGrid = NULL;        // pointer to main Grid object
};

typedef struct _N_VectorContent_Tempest *N_VectorContent_Tempest;


///// Part III -- Functions and macros exported by TempestNVector /////

// N_VNew_Tempest
// This function creates a new TempestNVector by locking a vector from the registry
N_Vector N_VNew_Tempest(Grid & grid);

#define NV_CONTENT_TEMPEST(v) ( (N_VectorContent_Tempest)(v->content) )
#define NV_GRID_TEMPEST(v) ( NV_CONTENT_TEMPEST(v)->mGrid )
#define NV_INDEX_TEMPEST(v) ( NV_CONTENT_TEMPEST(v)->mVectorIndex )


///// Part IV -- Required vector operations on a TempestNVector /////

N_Vector N_VClone_Tempest(N_Vector);
void N_VDestroy_Tempest(N_Vector);
void N_VConst_Tempest(realtype, N_Vector);
void N_VAbs_Tempest(N_Vector, N_Vector);
void N_VScale_Tempest(realtype, N_Vector, N_Vector);
void N_VAddConst_Tempest(N_Vector, realtype, N_Vector);
void N_VLinearSum_Tempest(realtype, N_Vector, realtype, N_Vector, N_Vector);
void N_VProd_Tempest(N_Vector, N_Vector, N_Vector);
void N_VDiv_Tempest(N_Vector, N_Vector, N_Vector);
void N_VInv_Tempest(N_Vector, N_Vector);
realtype N_VDotProd_Tempest(N_Vector, N_Vector);
realtype N_VMin_Tempest(N_Vector);
realtype N_VWrmsNorm_Tempest(N_Vector, N_Vector);
realtype N_VMaxNorm_Tempest(N_Vector);


#ifdef __cplusplus
}
#endif

#endif



/*********************** END OF FILE ***********************/
