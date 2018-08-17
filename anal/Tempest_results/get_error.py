#!/usr/bin/env python
# ===============================================================================
# Script to extract the maximum value of a given variable over time
#
# D.J. Gardner @ LLNL, April 2017
# ===============================================================================

def main():
    
    import argparse
    import os, sys
    import shlex

    from netCDF4 import Dataset
    import FileHelper

    import numpy as np
    
    parser = argparse.ArgumentParser(
        description='Extract global max for a given variable over time')

    parser.add_argument('RefFile', type=str,
                        help='path to reference solution output file')
    
    parser.add_argument('FileName', type=str,
                        help='Path to test parent directory')

    parser.add_argument('--Norm', dest='Norm',
                        choices=['RMS','L2','L1','Max'],
                        default='RMS',
                        help='norm to use for computing errors')

    # parse command line args
    args = parser.parse_args()
    
    # -------------------------------------------------------------------------------
    # Extract Data
    # -------------------------------------------------------------------------------

    refdata = Dataset(args.RefFile, mode="r")

    Uref = np.ravel(refdata.variables['U'][...])
    Vref = np.ravel(refdata.variables['V'][...])
    Wref = np.ravel(refdata.variables['W'][...])
    Rref = np.ravel(refdata.variables['Rho'][...])
    Tref = np.ravel(refdata.variables['Theta'][...])
    Sref = np.concatenate((Uref, Vref, Wref, Tref, Rref))
    Velref = np.concatenate((Uref, Vref, Wref))

    data = Dataset(args.FileName, mode="r")
     
    Uvar = np.ravel(data.variables['U'][...])
    Vvar = np.ravel(data.variables['V'][...])          
    Wvar = np.ravel(data.variables['W'][...])
    Rvar = np.ravel(data.variables['Rho'][...])
    Tvar = np.ravel(data.variables['Theta'][...])    
    Svar = np.concatenate((Uvar, Vvar, Wvar, Tvar, Rvar))
    Velvar = np.concatenate((Uvar, Vvar, Wvar))

    ErrU = err_norm(Uvar, Uref, args.Norm)
    ErrV = err_norm(Vvar, Vref, args.Norm)
    ErrW = err_norm(Wvar, Wref, args.Norm)
    ErrR = err_norm(Rvar, Rref, args.Norm)
    ErrT = err_norm(Tvar, Tref, args.Norm)
    ErrS = err_norm(Svar, Sref, args.Norm) 
    ErrVel = err_norm(Velvar, Velref, args.Norm)

#    print "U error: ", ErrU
#    print "V error: ", ErrV
#    print "W error: ", ErrW
    print "S error: ", ErrS
    print "R error: ", ErrR
    print "T error: ", ErrT
    print "V error: ", ErrVel
#    print "S error: ", ErrS

    data.close()

# ===============================================================================

def err_norm(data, ref, norm_type):

    import os
    import numpy as np

    if (norm_type == 'RMS'):
        return np.sqrt(np.average(np.square(data - ref)))

    elif (norm_type == 'L2'):
        return np.sqrt(np.sum(np.square(data - ref)))

    elif (norm_type == 'L1'):
        return np.sum(np.abs(data - ref))

    elif (norm_type == 'Max'):
        return np.amax(np.abs(data - ref))

    else:
        print "ERROR: Unknown norm type"
        sys.exit()

# ===============================================================================

if __name__ == "__main__":
    main()

# EOF

