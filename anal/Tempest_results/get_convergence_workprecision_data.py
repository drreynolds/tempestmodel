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

    parser.add_argument('--SkipErr', dest='SkipErr', 
                        action='store_true',
                        help='skip check for .err files')   

    parser.add_argument('--Debug', dest='Debug', 
                        action='store_true',
                        help='turn on debugging output')

    # parse command line args
    args = parser.parse_args()

    # -------------------------------------------------------------------------------
    # Input Checking
    # -------------------------------------------------------------------------------

    if (args.Debug):
        print args

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

            if (not args.SkipErr):

                errfiles = FileHelper.get_files_with_extension(t,'.err')

                if (len(errfiles) != 1):
                    print "ERROR:",len(errfiles),".err files found in",t
                    print errfiles
                    sys.exit()
    
    # -------------------------------------------------------------------------------
    # Load Reference Data
    # -------------------------------------------------------------------------------

    print "Reading Reference:",RefName

    RefData = Dataset(args.RefFile, mode="r")

    reflev  = RefData.variables['lev'][...]
    refilev = RefData.variables['ilev'][...]

    Uref = RefData.variables['U'][...]
    Vref = RefData.variables['V'][...]
    Wref = RefData.variables['W'][...]
    Rref = RefData.variables['Rho'][...]
    if (args.RhoTheta):
        Tref = RefData.variables['RhoTheta'][...]
    else:
        Tref = RefData.variables['Theta'][...]
                
    RefData.close()
    
    # interpolate reference data if necessary
    # assumes all tests use the same discretization settings

    # ncfile  = os.path.join(Tests[0], 'outTempest', RefName)
    # print ncfile
    # OutData = Dataset(ncfile, mode="r")

    # testlev  = OutData.variables['lev'][...]
    # testilev = OutData.variables['ilev'][...]

    # OutData.close()

    # if (not np.array_equal(reflev,testlev)):

    #     print "Interpolating Reference Solution..."
        
    #     Utmp = np.empty([1,len(testlev),180,360])
    #     Vtmp = np.empty([1,len(testlev),180,360])
    #     Rtmp = np.empty([1,len(testlev),180,360])
    #     Ttmp = np.empty([1,len(testlev),180,360])

    #     Wtmp = np.empty([1,len(testilev),180,360])

    #     # level variables
    #     for i in range(len(testlev)):
    #         for j in range(len(reflev)):
    #             if (testlev[i] < reflev[j]):
    #                 print i,reflev[j-1], testlev[i], reflev[j]
    #                 Utmp[0,i,:,:] = (Uref[0,j-1,:,:] + Uref[0,j,:,:])*0.5
    #                 Vtmp[0,i,:,:] = (Vref[0,j-1,:,:] + Vref[0,j,:,:])*0.5
    #                 Rtmp[0,i,:,:] = (Rref[0,j-1,:,:] + Rref[0,j,:,:])*0.5
    #                 Ttmp[0,i,:,:] = (Tref[0,j-1,:,:] + Tref[0,j,:,:])*0.5
    #                 break
    
    #     # interface variables
    #     for i in range(len(testilev)):
    #         for j in range(len(refilev)):
    #             if (testilev[i] == refilev[j]):
    #                 print i, testilev[i], refilev[j]
    #                 Wtmp[0,i,:,:] = Wref[0,j,:,:]
    #                 break
                
    #     # rename interpolated solutions
    #     Uref = Utmp
    #     Vref = Vtmp
    #     Wref = Wtmp
    #     Rref = Rtmp
    #     Tref = Ttmp

    # unravel arrays to 1D arrays
    Uref = np.ravel(Uref)
    Vref = np.ravel(Vref)
    Wref = np.ravel(Wref)
    Rref = np.ravel(Rref)
    Tref = np.ravel(Tref)

    Sref = np.concatenate((Uref, Vref, Wref, Tref, Rref))

    # -------------------------------------------------------------------------------
    # Load Test Data
    # -------------------------------------------------------------------------------
    if (args.Table):
        if (os.path.isfile(args.Table)):
            fout = open(args.Table,'a')
        else:
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

            if (not os.path.isfile(ncfile)):
                print "WARNING:",t,"does not contain",ncfile
                continue 

            outfile = FileHelper.get_files_with_extension(t,'.out')[0]

            if (not args.SkipErr):
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

            # initialize in case skipping
            minutes = 0.0
            seconds = 0.0

            if (not args.SkipErr):

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
        
            OutData.close()

            Sout = np.concatenate((Uout, Vout, Wout, Tout, Rout))

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
        if (args.Components):

            A,B,C,D,E,F,G,H = (list(d) for d in zip(*sorted(zip(stepsize,
                                                                ErrS,
                                                                runtime,
                                                                ErrU,
                                                                ErrV,
                                                                ErrW,
                                                                ErrT,
                                                                ErrR ))))
            A = np.array(A)
            B = np.array(B)
            C = np.array(C)

            D = np.array(D)
            E = np.array(E)
            F = np.array(F)
            G = np.array(G)
            H = np.array(H)

            fname = args.Label+"_"+TestParent.split('/')[-2]\
                +"_"+TestParent.split('/')[-1]+"_"+args.Norm+"errors_components.txt"

            ErrData = np.column_stack((A,B,C,D,E,F,G,H))
            np.savetxt(fname, ErrData)

            if (args.Table):

                fmt = '{0:<10s} {1:<10s} {2:<10f} {3:>22.16f} {4:>10f} {5:>14f} {6:>22.16f} {7:>22.16f} {8:>22.16f} {9:>22.16f} {10:>22.16f}'
                for i in range(0,len(A)):
                    if (i < len(A)-1):
                        cord = np.log10(B[i+1]/B[i])/np.log10(A[i+1]/A[i])
                    else:
                        cord = 0.0
                    print >> fout, fmt.format(integrator[0], method[0],
                                              A[i], B[i], cord, C[i],
                                              D[i], E[i], F[i], G[i], H[i],)
                print >> fout, '-'*200

        else:

            A, B, C = (list(d) for d in zip(*sorted(zip(stepsize, ErrS, runtime))))
            A = np.array(A)
            B = np.array(B)
            C = np.array(C)

            fname = args.Label\
                +"_"+TestParent.split('/')[-3]\
                +"_"+TestParent.split('/')[-2]\
                +"_"+TestParent.split('/')[-1]+"_"+args.Norm+"errors.txt"

            ErrData = np.column_stack((A, B, C))
            np.savetxt(fname, ErrData)

            if (args.Table):

                fmt = '{0:<10s} {1:<10s} {2:<10f} {3:>22.16f} {4:>10f} {5:>14f}'
                for i in range(0,len(A)):
                    if (i < len(A)-1):
                        cord = np.log10(B[i+1]/B[i])/np.log10(A[i+1]/A[i])
                    else:
                        cord = 0.0
                    print >> fout, fmt.format(integrator[0], method[0],
                                              A[i], B[i], cord, C[i])
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

