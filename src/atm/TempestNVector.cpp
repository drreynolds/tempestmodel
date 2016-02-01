///////////////////////////////////////////////////////////////////////////////
///
///	\file    TempestNVector.cpp
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

#include "TempestNVector.h"
#include "Model.h"
#include "Grid.h"
#include <sundials/sundials_math.h>

///////////////////////////////////////////////////////////////////////////////

#define ZERO RCONST(0.0)


////////////////////// TempestNVector registry /////////////////////////

bool TempestNVectorRegistry[MAX_TEMPEST_NVECTORS];   // global bool arrays initialize to false by default


////////////////////// Exported Functions /////////////////////////

// Function to create a new TempestNVector
N_Vector N_VNew_Tempest(Grid & grid) {

  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_Tempest content;

  // Initialize vector pointers to NULL
  v = NULL;
  ops = NULL;
  content = NULL;
  
  // Create vector
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);

  // Create vector operation structure
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) {free(v); return(NULL); }

  // Attach implemented vector routines to N_Vector_Ops structure
  ops->nvclone     = N_VClone_Tempest;
  ops->nvdestroy   = N_VDestroy_Tempest;
  ops->nvconst     = N_VConst_Tempest;
  ops->nvabs       = N_VAbs_Tempest;
  ops->nvscale     = N_VScale_Tempest;
  ops->nvaddconst  = N_VAddConst_Tempest;
  ops->nvlinearsum = N_VLinearSum_Tempest;
  ops->nvprod      = N_VProd_Tempest;
  ops->nvdiv       = N_VDiv_Tempest;
  ops->nvinv       = N_VInv_Tempest;
  ops->nvdotprod   = N_VDotProd_Tempest;
  ops->nvmin       = N_VMin_Tempest;
  ops->nvwrmsnorm  = N_VWrmsNorm_Tempest;
  ops->nvmaxnorm   = N_VMaxNorm_Tempest;

  // Signal that remaining vector routines are not implemented
  ops->nvspace           = NULL;
  ops->nvgetarraypointer = NULL;
  ops->nvsetarraypointer = NULL;
  ops->nvwrmsnormmask    = NULL;
  ops->nvcloneempty      = NULL;
  ops->nvwl2norm         = NULL;
  ops->nvl1norm          = NULL;
  ops->nvcompare         = NULL;
  ops->nvinvtest         = NULL;
  ops->nvconstrmask      = NULL;
  ops->nvminquotient     = NULL;

  // Create content
  content = 
    (N_VectorContent_Tempest) malloc(sizeof(struct _N_VectorContent_Tempest));
  if (content == NULL) {free(ops); free(v); return(NULL);}

  // Set mVectorIndex to first available vector from registry (if any are available)
  bool successful=false;
  for (int i=0; i<MAX_TEMPEST_NVECTORS; i++) {
    if (!TempestNVectorRegistry[i]) {
      content->mVectorIndex = i;
      TempestNVectorRegistry[i] = true;
      successful=true;
    }
  }
  if (!successful) {
    free(content);
    free(ops);
    free(v);
    return(NULL);
  }

  // store reference to Grid object
  content->mGrid = &grid;

  // Attach content and ops to generic N_Vector
  v->content = content;
  v->ops     = ops;

  return(v);
}


////////////////// Required Vector operations /////////////////////

// N_VClone_Tempest returns a new N_Vector of the same form 
// as the input N_Vector
N_Vector N_VClone_Tempest(N_Vector w) {

  // initialize output to NULL
  N_Vector v = NULL;

  // initialize Grid pointer to NULL
  Grid * grid = NULL;

  // attempt access to w's mGrid object
  grid = NV_GRID_TEMPEST(w);
  if (grid == NULL) return(NULL);

  // Create vector and return
  v = N_VNew_Tempest(*grid);

  return(v);
}


// N_VDestroy_Tempest unlocks the data associated with the NVector, 
// and frees the NVector data storage
void N_VDestroy_Tempest(N_Vector v) {
  TempestNVectorRegistry[NV_INDEX_TEMPEST(v)] = false;
  free(v->content); v->content = NULL;
  free(v->ops); v->ops = NULL;
  free(v); v = NULL;
  return;
}


// N_VConst_Tempest (or nvconst) calculates z[i] = c for i = 0,...,N-1
void N_VConst_Tempest(realtype c, N_Vector z) {
  // access Grid and vector index associated with z
  Grid * grid = NULL;
  grid = NV_GRID_TEMPEST(z);
  int zID = NV_INDEX_TEMPEST(z);

  // call Grid routine to do the work
  grid->ConstantData(c, zID);
  return;
}


// N_VAbs_Tempest or (nvabs) calculates z[i] = fabs(x[i]) for i = 0,...,N-1
void N_VAbs_Tempest(N_Vector x, N_Vector z) {
  // access Grid and vector indices associated with x and z
  Grid * grid = NULL;
  grid = NV_GRID_TEMPEST(x);
  int xID = NV_INDEX_TEMPEST(x);
  int zID = NV_INDEX_TEMPEST(z);

  // call Grid routine to do the work
  grid->AbsData(xID, zID);
  return;
}


// N_VScale_Tempest (or nvscale) calculates z[i] = c*x[i] for i = 0,...,N-1
void N_VScale_Tempest(realtype c, N_Vector x, N_Vector z) {
  // access Grid and vector indices associated with x and z
  Grid * grid = NULL;
  grid = NV_GRID_TEMPEST(x);
  int xID = NV_INDEX_TEMPEST(x);
  int zID = NV_INDEX_TEMPEST(z);

  // call Grid routine to do the work
  grid->ScaleData(c, xID, zID);
  return;
}


// N_VAddConst_Tempest (or nvaddconst) calculates z[i] = x[i] + b for i = 0,...,N-1
void N_VAddConst_Tempest(N_Vector x, realtype b, N_Vector z) {
  // access Grid and vector indices associated with x and z
  Grid * grid = NULL;
  grid = NV_GRID_TEMPEST(x);
  int xID = NV_INDEX_TEMPEST(x);
  int zID = NV_INDEX_TEMPEST(z);

  // call Grid routine to do the work
  grid->AddConstantData(xID, b, zID);
  return;
}


// N_VLinearSum_Tempest (or nvlinearsum) calculates z[i] = a*x[i] + b*y[i] for i = 0,...,N-1
void N_VLinearSum_Tempest(realtype a, N_Vector x, realtype b, 
			  N_Vector y, N_Vector z) {
  // access Grid and vector indices associated with x, y and z
  Grid * grid = NULL;
  grid = NV_GRID_TEMPEST(x);
  int xID = NV_INDEX_TEMPEST(x);
  int yID = NV_INDEX_TEMPEST(y);
  int zID = NV_INDEX_TEMPEST(z);

  // call Grid routine to do the work
  grid->LinearSumData(a, xID, b, yID, zID);
  return;
}


// N_VProd_Tempest (or nvprod) calculates z[i] = x[i]*y[i] for i = 0,...,N-1
void N_VProd_Tempest(N_Vector x, N_Vector y, N_Vector z) {
  // access Grid and vector indices associated with x, y and z
  Grid * grid = NULL;
  grid = NV_GRID_TEMPEST(x);
  int xID = NV_INDEX_TEMPEST(x);
  int yID = NV_INDEX_TEMPEST(y);
  int zID = NV_INDEX_TEMPEST(z);

  // call Grid routine to do the work
  grid->ProductData(xID, yID, zID);
  return;
}


// N_VDiv_Tempest (or nvdiv) calculates z[i] = x[i]/y[i] for i = 0,...,N-1
void N_VDiv_Tempest(N_Vector x, N_Vector y, N_Vector z) {
  // access Grid and vector indices associated with x, y and z
  Grid * grid = NULL;
  grid = NV_GRID_TEMPEST(x);
  int xID = NV_INDEX_TEMPEST(x);
  int yID = NV_INDEX_TEMPEST(y);
  int zID = NV_INDEX_TEMPEST(z);

  // call Grid routine to do the work
  grid->DivideData(xID, yID, zID);
  return;
}


// N_VInv_Tempest (or nvinv) calculates z[i] = 1.0/x[i]  for i = 0,...,N-1
// Note: it does not check for division by 0.  It should be called only 
// with an N_Vector x which is guaranteed to have all non-zero components.
void N_VInv_Tempest(N_Vector x, N_Vector z) {
  // access Grid and vector indices associated with x and z
  Grid * grid = NULL;
  grid = NV_GRID_TEMPEST(x);
  int xID = NV_INDEX_TEMPEST(x);
  int zID = NV_INDEX_TEMPEST(z);

  // call Grid routine to do the work
  grid->InvertData(xID, zID);
  return;
}


// N_VDotProd_Tempest (or nvdotprod) returns the value of the 
//  ordinary dot product of x and y, i.e. sum_{i=0,...,N-1} (x[i] * y[i])
realtype N_VDotProd_Tempest(N_Vector x, N_Vector y) {
  // access Grid and vector indices associated with x and y
  Grid * grid = NULL;
  grid = NV_GRID_TEMPEST(x);
  int xID = NV_INDEX_TEMPEST(x);
  int yID = NV_INDEX_TEMPEST(y);

  // call Grid routine to do the work
  realtype sum = grid->DotProductData(xID, yID);
  return(sum);
}


// N_VMin_Tempest (or nvmin) returns the smallest element of x
realtype N_VMin_Tempest(N_Vector x) {
  // access Grid and vector indices associated with x
  Grid * grid = NULL;
  grid = NV_GRID_TEMPEST(x);
  int xID = NV_INDEX_TEMPEST(x);

  // call Grid routine to do the work
  realtype minval = grid->MinimumData(xID);
  return(minval);
}

 
// N_VWrmsNorm_Tempest (or nvwrmsnorm) returns the weighted root 
// mean square norm of x with weight factor w, 
// i.e. sqrt [(sum_{i=0,...,N-1} (x[i] * w[i])^2) / N]
realtype N_VWrmsNorm_Tempest(N_Vector x, N_Vector w) {
  // access Grid and vector indices associated with x and w
  Grid * grid = NULL;
  grid = NV_GRID_TEMPEST(x);
  int xID = NV_INDEX_TEMPEST(x);
  int wID = NV_INDEX_TEMPEST(w);

  // call Grid routine to do the work
  realtype wrmsval = grid->WRMSNormData(xID, wID);
  return(wrmsval);
}


// N_VMaxNorm_Tempest (or nvmaxnorm) returns the maximum norm of x, 
// i.e. max_{i=0,...,N-1} abs(x[i])
realtype N_VMaxNorm_Tempest(N_Vector x) { 
  // access Grid and vector indices associated with x
  Grid * grid = NULL;
  grid = NV_GRID_TEMPEST(x);
  int xID = NV_INDEX_TEMPEST(x);

  // call Grid routine to do the work
  realtype maxval = grid->MaximumNormData(xID);
  return(maxval);
}


//////////////////////// END OF FILE ////////////////////////
