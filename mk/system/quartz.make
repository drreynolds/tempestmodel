# Copyright (c) 2016      Bryce Adelstein-Lelbach aka wash
# Copyright (c) 2000-2016 Paul Ullrich 
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying 
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# Cab System at LLNL

CXX=    mpicxx
F90=    mpif90
MPICXX= mpicxx
MPIF90= mpif90

# SUNDIALS
SUNDIALS_ROOT=$(HOME)/local/sundials-2.7.0_quartz_gnu_opt
SUNDIALS_INCLUDEDIR=$(SUNDIALS_ROOT)/include
SUNDIALS_LIBDIR=$(SUNDIALS_ROOT)/lib

# NetCDF
NETCDF_ROOT=      $(HOME)/local/netcdf-4.1.3_quartz_gnu_opt
NETCDF_CXXFLAGS=  -I$(NETCDF_ROOT)/include
NETCDF_LIBRARIES= -lnetcdf -lnetcdf_c++
NETCDF_LDFLAGS=   -L$(NETCDF_ROOT)/lib

# LAPACK
LAPACK_INTERFACE= FORTRAN
LAPACK_CXXFLAGS= -g
LAPACK_LIBRARIES= -llapack -lblas
LAPACK_LDFLAGS=

# DO NOT DELETE
