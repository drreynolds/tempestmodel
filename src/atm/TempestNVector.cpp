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

// require SUNDIALS for compilation
#ifdef USE_SUNDIALS

//#define DEBUG_OUTPUT

#include "TempestNVector.h"
#include "Model.h"
#include "Announce.h"
#include "Grid.h"
#include <sundials/sundials_math.h>

///////////////////////////////////////////////////////////////////////////////

#define ZERO RCONST(0.0)


////////////////////// TempestNVector registry /////////////////////////

bool TempestNVectorRegistry[MAX_TEMPEST_NVECTORS];   // global bool arrays initialize to false by default

int MaxRegistryLength = 0;

////////////////////// Exported Functions /////////////////////////

// Function to set the maximum number of vectors in the registry
int SetMaxTempestNVectorRegistryLength(int nvec) {

  MaxRegistryLength = nvec;

  if (MaxRegistryLength > MAX_TEMPEST_NVECTORS)
    return -1;

  return 0;
}

// Function to get the number of vectors in the registry
int GetTempestNVectorRegistryLength() {

  int RegistryLength;

  for (int i=0; i<MaxRegistryLength; i++) {
    if (!TempestNVectorRegistry[i]) {
      RegistryLength = i;
      break;
    }
  }

  return RegistryLength;
}

// Marks first available vector from registry as used (if any are available)
int ReserveNextTempestNVectorRegistryIdx() {
  
  for (int i=0; i<MaxRegistryLength; i++) {
    if (!TempestNVectorRegistry[i]) {
      TempestNVectorRegistry[i] = true;
      return i;
    }
  }
  // no places available
  return -1;
}

// Function to create a new TempestNVector
N_Vector N_VNew_Tempest(Grid & grid, Model & model) {

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
  for (int i=0; i<MaxRegistryLength; i++) {
    if (!TempestNVectorRegistry[i]) {
      content->mVectorIndex = i;
      TempestNVectorRegistry[i] = true;
      successful=true;
      break;
    }
  }
  if (!successful) {
    free(content);
    free(ops);
    free(v);
    return(NULL);
  }

#ifdef DEBUG_OUTPUT
  // Get process rank
  int iRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &iRank);

  if (iRank == 0) {
    std::cout << std::endl
	      << " Proc "
	      << iRank
	      << " N_Vector created at Registry Index: "
	      << content->mVectorIndex
	      << std::endl;
  }
#endif
  
  // store reference to Grid object
  content->mGrid = &grid;

  // store reference to Model object
  content->mModel = &model;

  // Attach content and ops to generic N_Vector
  v->content = content;
  v->ops     = ops;

  return(v);
}


// Function to attach an existing Tempest "state" to a new TempestNVector
N_Vector N_VAttach_Tempest(Grid & grid, Model & model, int VectorIndex) {

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

  // Set mVectorIndex to input vector index
  content->mVectorIndex = VectorIndex;
  
  // store reference to Grid object
  content->mGrid = &grid;

  // store reference to Model object
  content->mModel = &model;

  // Attach content and ops to generic N_Vector
  v->content = content;
  v->ops     = ops;

  return(v);
}


// Function to run tests to verify a working TempestNVector
void N_VTest_Tempest(N_Vector x) {

  Announce("Testing TempestNVector");
  
  // Create a few extra NVectors
  Announce("  Testing N_VClone");
  N_Vector w = N_VClone(x);
  N_Vector y = N_VClone(x);
  N_Vector z = N_VClone(x);

  // set NVector values to different constants
  N_VConst( 1.0, w);    // w[i] = 1
  N_VConst( 2.0, x);    // x[i] = 2
  N_VConst(-1.0, y);    // y[i] = -1
  N_VConst(-2.0, z);    // z[i] = -2
  
  // test N_VDotProd by computing number of vector entries different ways
  Announce("  Testing N_VDotProd:");
  if ((fabs(N_VDotProd(w,w) - N_VDotProd(y,y)) < 1e-15) && 
      (fabs(N_VDotProd(w,w) - 0.25*N_VDotProd(x,x)) < 1e-15) && 
      (fabs(N_VDotProd(y,y) - 0.25*N_VDotProd(z,z)) < 1e-15) &&
      (fabs(N_VDotProd(x,x) - N_VDotProd(z,z)) < 1e-15))
    Announce("    success");
  else
    Announce("    failure");

  // test N_VMaxNorm
  Announce("  Testing N_VMaxNorm:");
  if ((fabs(N_VMaxNorm(w) - 1.0) < 1e-15) && 
      (fabs(N_VMaxNorm(x) - 2.0) < 1e-15) && 
      (fabs(N_VMaxNorm(y) - 1.0) < 1e-15) &&
      (fabs(N_VMaxNorm(z) - 2.0) < 1e-15))
    Announce("    success");
  else
    Announce("    failure");

  // test N_VMin
  Announce("  Testing N_VMin:");
  if ((fabs(N_VMin(w) - 1.0) < 1e-15) && 
      (fabs(N_VMin(x) - 2.0) < 1e-15) && 
      (fabs(N_VMin(y) + 1.0) < 1e-15) &&
      (fabs(N_VMin(z) + 2.0) < 1e-15))
    Announce("    success");
  else
    Announce("    failure");

  // get number of entries in NVectors
  int n = int(N_VDotProd(w,w));
  Announce("  Number of entries in vectors = %i", n);

  // test N_VWrmsNorm
  Announce("  N_VWrmsNorm tests: w[i] = 1, x[i] = 2, y[i] = -1, z[i] = -2: ");
  Announce("    N_VWrmsNorm(x,w) = %g", N_VWrmsNorm(x,w));
  Announce("    N_VWrmsNorm(y,w) = %g", N_VWrmsNorm(y,w));
  Announce("    N_VWrmsNorm(z,w) = %g", N_VWrmsNorm(z,w));

  // test N_VAbs
  Announce("  Testing N_VAbs:");
  N_VConst(-2.0, z);
  N_VAbs(z, y);   // y[i] = 2
  if (((sqrt(N_VDotProd(y,y)/n) - fabs(N_VMin(z))) < 1e-15) &&
      (fabs(N_VMin(y) - 2.0) < 1e-15))
    Announce("    success");
  else {
    Announce("    failure (y[i] should equal %g)",fabs(N_VMin(z)));
    Announce("      sqrt(N_VDotProd(y,y)/n) = %g",sqrt(N_VDotProd(y,y)/n));
    Announce("      N_VMin(y) = %g",N_VMin(y));
  }

  // test N_VLinearSum
  Announce("  Testing N_VLinearSum:");
  N_VConst( 5.0, w);                   // w[i] = 5
  N_VConst(-1.0, x);                   // x[i] = -1
  N_VLinearSum( 2.0, w, 10.0, x, y);   // y[i] = 0
  if (N_VMaxNorm(y) < 1e-15)
    Announce("    success");
  else
    Announce("    failure");
  
  // test N_VScale
  Announce("  Testing N_VScale:");
  N_VConst(3.0, w);                    // w[i] = 3
  N_VScale(4.0, w, x);                 // x[i] = 12
  N_VConst(-12.0, y);                  // y[i] = -12
  N_VLinearSum( 1.0, x, 1.0, y, z);    // z[i] = 0
  if (N_VMaxNorm(z) < 1e-15)
    Announce("    success");
  else
    Announce("    failure");
  
  // test N_VAddConst
  Announce("  Testing N_VAddConst:");
  N_VConst(7.0, w);                    // w[i] = 7
  N_VAddConst(w, -2.0, x);             // x[i] = 5
  N_VConst(5.0, y);                    // y[i] = 5
  N_VLinearSum( 1.0, x, -1.0, y, z);   // z[i] = 0
  if (N_VMaxNorm(z) < 1e-15)
    Announce("    success");
  else {
    Announce("    failure (x[i] and y[i] should both equal 5)");
    Announce("      N_VMin(x) = %g",N_VMin(x));
    Announce("      N_VMin(y) = %g",N_VMin(y));
  }
  
  // test N_VProd
  Announce("  Testing N_VProd:");
  N_VConst(5.0, w);                    // w[i] = 5
  N_VConst(0.2, x);                    // x[i] = 1/5
  N_VProd(w, x, y);                    // y[i] = 1
  N_VConst(1.0, x);                    // x[i] = 1
  N_VLinearSum( 1.0, y, -1.0, x, z);   // z[i] = 0
  if (N_VMaxNorm(z) < 1e-15)
    Announce("    success");
  else
    Announce("    failure");
  
  // test N_VDiv
  Announce("  Testing N_VDiv:");
  N_VConst(6.0, w);                    // w[i] = 6
  N_VConst(-2.0, x);                   // x[i] = -2
  N_VDiv(w, x, y);                     // y[i] = -3
  N_VConst(3.0, x);                    // x[i] = 3
  N_VLinearSum( 1.0, y, 1.0, x, z);    // z[i] = 0
  if (N_VMaxNorm(z) < 1e-15)
    Announce("    success");
  else
    Announce("    failure");
  
  // test N_VInv
  Announce("  Testing N_VInv:");
  N_VConst(4.0, w);                    // w[i] = 4
  N_VInv(w, x);                        // x[i] = 1/4
  N_VConst(-0.25, y);                  // y[i] = -1/4
  N_VLinearSum( 1.0, x, 1.0, y, z);    // z[i] = 0
  if (N_VMaxNorm(z) < 1e-15)
    Announce("    success");
  else
    Announce("    failure");

  // Free temporary NVectors
  Announce("  Testing N_VDestroy");
  N_VDestroy(w);
  N_VDestroy(y);
  N_VDestroy(z);

  Announce("TempestNVector tests complete");
  return;
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

  // initialize Model pointer to NULL
  Model * model = NULL;

  // attempt access to w's mModel object
  model = NV_MODEL_TEMPEST(w);
  if (model == NULL) return(NULL);

  // Create vector and return
  v = N_VNew_Tempest(*grid, *model);

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

#endif


//////////////////////// END OF FILE ////////////////////////
