# Copyright (c) 2016      Bryce Adelstein-Lelbach aka wash
# Copyright (c) 2000-2016 Paul Ullrich 
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying 
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# Dan Reynolds' laptop

CXX=    mpicxx
F90=    mpif90
MPICXX= mpicxx
MPIF90= mpif90

# SUNDIALS
SUNDIALS_INCLUDEDIR= /home/dreynolds/research/ARKode/install_opt/include
SUNDIALS_LIBDIR= /home/dreynolds/research/ARKode/install_opt/lib
#SUNDIALS_INCLUDEDIR= /home/dreynolds/research/ARKode/install_dbg/include
#SUNDIALS_LIBDIR= /home/dreynolds/research/ARKode/install_dbg/lib

# NetCDF
NETCDF_ROOT=      /usr/local/netcdf-cxx-4.2
NETCDF_CXXFLAGS=  -I$(NETCDF_ROOT)/include
NETCDF_LIBRARIES= -lnetcdf -lnetcdf_c++
NETCDF_LDFLAGS=   -L$(NETCDF_ROOT)/lib

# LAPACK
LAPACK_INTERFACE= FORTRAN
LAPACK_CXXFLAGS=
LAPACK_LIBRARIES= -llapack -lblas
LAPACK_LDFLAGS=

# DO NOT DELETE
