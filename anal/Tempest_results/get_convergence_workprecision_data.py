#!/usr/bin/env python
# ===============================================================================
# Script to extract the errors and runtimes from multiple runs compared to a 
# reference solution for making convergence and work-precision plots.
#
# This script assumes a certain file structure and naming conventions. The input
# 'TestDirs' is the path to a parent directory (or multiple directories) which 
# possibly contains multiple individual test directories each of which has test
# output in a 'outTempest' directory. Example directory structure with parents
# ARS343 and ARS232 each with individual tests at different step sizes:
#
# ARS343/
#     dt100/
#         outTempest/
#     dt200/
#         outTempest/
#     dt400/
#         outTempest/
#
# ARS232/
#     dt100/
#         outTempest/
#     dt200/
#         outTempest/
#     dt400/
#         outTempest/
#
# The output from this script is a text file with the step sizes, errors, and 
# run times for each parent directory.
#
# Example call for directory structure above:
#
# tempest_get_errors.py Convergence_Tests \
#     ~/reference/outTempest/out.0000-01-01-03600.nc \
#     ~/tests/ARS343 \
#     ~/tests/ARS232 \
#     --Table Convergence_Tests.txt
#
# D.J. Gardner @ LLNL, Mar 2017
# ===============================================================================

def main():

    import argparse
    import os, sys
    import shlex

    from netCDF4 import Dataset
    import FileHelper

    import numpy as np
    
    parser = argparse.ArgumentParser(
        description='Extract runtime and error with respect to a reference'+
        ' solution and write data to a text file')

    parser.add_argument('Label', type=str,
                        help='label used in output file name')
        
    parser.add_argument('RefFile', type=str,
                        help='path to reference solution output file')

    parser.add_argument('TestDirs', type=str, nargs='+',
                        help='path to test parent directories')
    
    parser.add_argument('--Norm', dest='Norm', 
                        choices=['RMS','L2','L1','Max'],
                        default='RMS',
                        help='norm to use for computing errors')
    
    parser.add_argument('--Components', dest='Components', 
                        action='store_true', 
                        help='compute error for individual state variables')

    parser.add_argument('--RhoTheta', dest='RhoTheta', 
                        action='store_true',
                        help='flag to indicate results use rho-theta formulation')

    parser.add_argument('--Table', dest='Table', type=str,
                        help='file name for table of output')

    parser.add_argument('--Debug', dest='Debug', 
                        action='store_true',
                        help='turn on debugging output')

    # parse command line args
    args = parser.parse_args()

    # -------------------------------------------------------------------------------
    # Input Checking
    # -------------------------------------------------------------------------------

    # does reference file exist?
    if (not os.path.isfile(args.RefFile)):
        print "ERROR:",args.RefFile,"does not exist"
        sys.exit()

    # reference file name
    RefName = args.RefFile.split('/')[-1]

    # checks on test parent directories
    for TestParent in args.TestDirs:
        
        # do test directories exist?
        if (not os.path.isdir(TestParent)):
            print "ERROR:",TestParent,"does not exist"
            sys.exit()

        # test directories in test parent drectory
        Tests = FileHelper.get_immediate_subdirectories(TestParent)

        # checks on test directories
        for t in Tests:

            # do the tests directories contain an output file at the same time as
            # the reference solution file?
            fname = os.path.join(t, 'outTempest', RefName)

            if (not os.path.isfile(fname)):
                print "WARNING:",t,"does not contain",RefName
                continue 

            # is there only one .out and .err file? if not then this test was likely
            # restarted and needs some post-processing to combine .out and .err files
            outfiles = FileHelper.get_files_with_extension(t,'.out')
            
            if (len(outfiles) != 1):
                print "ERROR:",len(outfiles),".out files found in",t
                print outfiles
                sys.exit()

            errfiles = FileHelper.get_files_with_extension(t,'.err')

            if (len(errfiles) != 1):
                print "ERROR:",len(errfiles),".out files found in",t
                print errfiles
                sys.exit()
    
    # -------------------------------------------------------------------------------
    # Load Reference Data
    # -------------------------------------------------------------------------------

    print "Reading Reference:",RefName

    RefData = Dataset(args.RefFile, mode="r")

    Uref = np.ravel(RefData.variables['U'][...])
    Vref = np.ravel(RefData.variables['V'][...])
    Wref = np.ravel(RefData.variables['W'][...])
    Rref = np.ravel(RefData.variables['Rho'][...])
    if (args.RhoTheta):
        Tref = np.ravel(RefData.variables['RhoTheta'][...])
    else:
        Tref = np.ravel(RefData.variables['Theta'][...])

    Sref = np.concatenate((Uref, Vref, Wref, Tref, Rref))
                
    RefData.close()

    # -------------------------------------------------------------------------------
    # Load Test Data
    # -------------------------------------------------------------------------------
    if (args.Table):
        fout = open(args.Table,'w')
        
    for TestParent in args.TestDirs:

        print "Test Parent:",TestParent
        print >> fout, "Test Parent:",TestParent

        # create empty data lists
        integrator = []
        method     = []
        runtime    = []
        stepsize   = []        
        ErrS       = [] 
    
        if (args.Components):
            ErrU = []
            ErrV = []
            ErrW = []
            ErrR = []
            ErrT = []

        # get individual test directories
        Tests = FileHelper.get_immediate_subdirectories(TestParent)

        # loop over each test case
        for t in Tests:

            print "\t",t.split('/')[-1]

            ncfile  = os.path.join(t, 'outTempest', RefName)
            outfile = FileHelper.get_files_with_extension(t,'.out')[0]
            errfile = FileHelper.get_files_with_extension(t,'.err')[0]

            # check that run finished
            RunFinished = False

            with FileHelper.File(outfile) as fn:

                for line in fn.backward():
                    split_line = shlex.split(line)

                    if ("RESULTS" in split_line):
                        RunFinished = True
                        break
                    
                    if ("Step" in split_line):
                        break

            if (not RunFinished):
                print "WARNING: Run did not complete for",t
                continue

            # get run time
            FoundTime = False

            with open(errfile) as fn:

                for line in fn:
                    split_line = shlex.split(line)

                    if ("real" in split_line):
                        FoundTime = True

                        time = split_line[1]

                        j = 0
                        for c in time:
                            if (c == 'm'):
                                minutes = float(time[:j])
                                seconds = float(time[j+1:-1])
                                break
                            j += 1

            if (not FoundTime):
                print "WARNING: Run time not found for",t
                continue

            runtime.append(60.0 * minutes + seconds)

            # get run settings
            with open(outfile) as fn:

                for line in fn:
                    split_line = shlex.split(line)

                    if ("--timescheme" in split_line):
                        if (split_line[2] != "[arkode]"):
                            integrator.append("tempest")
                            method.append(split_line[2][1:-1])
                
                    if ("--arkode_butchertable" in split_line):
                        integrator.append("arkode")
                        method.append(split_line[2][1:-1])

                    if ("--dt" in split_line):

                        dt = split_line[2][1:-1]

                        if (dt[-1] == 'u'):
                            stepsize.append(1e-6 * float(dt[:-1]))
                        elif (dt[-1] == 's'):
                            stepsize.append(float(dt[:-1]))
                        else:
                            print "ERROR: Unknown time step size units in",t
                            sys.exit()

                    if ("MODEL SETUP" in line):
                        break

            # get run data
            OutData = Dataset(ncfile, mode="r")
               
            Uout = np.ravel(OutData.variables['U'][...])
            Vout = np.ravel(OutData.variables['V'][...])
            Wout = np.ravel(OutData.variables['W'][...])
            Rout = np.ravel(OutData.variables['Rho'][...])            
            if (args.RhoTheta):
                Tout = np.ravel(OutData.variables['RhoTheta'][...])
            else:
                Tout = np.ravel(OutData.variables['Theta'][...])

            Sout = np.concatenate((Uout, Vout, Wout, Tout, Rout))
        
            OutData.close()

            # compute error
            ErrS.append(err_norm(Sout, Sref, args.Norm))

            if (args.Components):
                ErrU.append(err_norm(Uout, Uref, args.Norm)) 
                ErrV.append(err_norm(Vout, Vref, args.Norm)) 
                ErrW.append(err_norm(Wout, Wref, args.Norm)) 
                ErrT.append(err_norm(Tout, Tref, args.Norm)) 
                ErrR.append(err_norm(Rout, Rref, args.Norm)) 

        # 
        # end individual test directory loop
        #

        # sort data by step size and write data to file(s)
        X, Y, Z = (list(d) for d in zip(*sorted(zip(stepsize, ErrS, runtime))))    
        X = np.array(X)
        Y = np.array(Y)
        Z = np.array(Z)

        ErrData = np.column_stack((X, Y, Z))

        fname = args.Label+"_"+TestParent.split('/')[-1]+"_errors.txt"

        np.savetxt(fname, ErrData)
        
        if (args.Table):

            fmt = '{0:<10s} {1:<10s} {2:<10f} {3:>22.16f} {4:>10f} {5:>14f}'            
            for i in range(0,len(X)):
                if (i < len(X)-1):
                    cord = np.log10(Y[i+1]/Y[i])/np.log10(X[i+1]/X[i])
                else:
                    cord = 0.0                          
                print >> fout, fmt.format(integrator[0], method[0], 
                                          X[i], Y[i], cord, Z[i])
            print >> fout, '-'*80

    # 
    # end parent directory loop
    #
            
    if (args.Table):
        # close output file
        fout.close()


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

