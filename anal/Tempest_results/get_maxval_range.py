#!/usr/bin/env python
# ===============================================================================
# Script to compute a confidence interval over time for the maximum value of a 
# given variable using the student's t-distribution. The output index, mean max 
# variable value, and deviation from the mean max value for the requested 
# confidence level are written to a text file.
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
    import scipy.stats as st
    
    parser = argparse.ArgumentParser(
        description='Extract max vaule over time from multiple runs and compute'+
        ' a confidence intervals for the max value over time')

    parser.add_argument('var', type=str,
                        help='varibale to extract from netcdf file')
    
    parser.add_argument('alpha', type=float,
                        help='confidence level, e.g. 0.99 for 99%')   
    
    parser.add_argument('tests', type=str, nargs='+',
                        help='paths to test output files')

    parser.add_argument('--debug', dest='debug', action='store_true',
                        help='turn on debugging output')

    # parse command line args
    args = parser.parse_args()
            
    # ---------------------------------------------------------------------------
    # load test data
    # ---------------------------------------------------------------------------      
    
    # number of tests in sample
    ntests = len(args.tests)

    for itest in range(ntests):

        print "Test:",args.tests[itest]

        # check that test directory exists
        if (not os.path.isdir(args.tests[itest])):
            print "ERROR: ",args.tests[itest]," does not exist"            
            sys.exit()

        # get all netcdf output files in the test directory
        outfiles  = FileHelper.list_files(args.tests[itest],'*.nc')
        noutfiles = len(outfiles)

        # exit if no output files are found
        if (noutfiles == 0):
            print "ERROR: No output files in",args.tests[itest]
            sys.exit()
        
        outfiles = sorted(outfiles, key=FileHelper.numerical_sort)

        # allocate array for max data values over time for each test case
        # NOTE: assumed that all tests have the same number of output files
        # and that all output are at the same simulation time
        if (itest == 0):
            maxvals = np.empty([ntests, noutfiles])

        # extract max value in each output file
        for iout in range(noutfiles):
        
            if (args.debug):
                print outfiles[iout]

            data = Dataset(outfiles[iout], mode="r")
            
            maxvals[itest][iout] = np.amax(np.abs(np.ravel(data.variables[args.var][...])))

            data.close()

    # compute mean and standard deviation for each column (output time)
    mean = np.mean(maxvals, axis=0)
    sdev = np.std(maxvals, axis=0)

    if (args.debug):
        print "mean    =",mean
        print "std dev =",sdev
           
    # critical t value for confidence interval with area alpha
    tc = st.t.interval(args.alpha, ntests)[1]

    # confidence interval deviation from mean value
    dev = tc * sdev / np.sqrt(ntests)

    if (args.debug):
        print "tc  =",tc
        print "dev =",dev

    # output times (assuming output a 1 day intervals)
    times = np.arange(0, noutfiles)

    if (args.debug):
        print np.shape(times)
        print np.shape(mean)
        print np.shape(dev)

    # create matrix of data for output
    outdata = np.column_stack((times, mean, dev))

    # write data to file
    np.savetxt("max"+args.var+"_range.txt", outdata)

# ===============================================================================

if __name__ == "__main__":
    main()

# EOF

