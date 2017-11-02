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
    
    parser.add_argument('Label', type=str,
                        help='label used in output file name')

    parser.add_argument('var', type=str,
                        help='varibale to extract from netcdf file')
   
    parser.add_argument('TestDir', type=str,
                        help='Path to test parent directory')

    parser.add_argument('--CheckRange', dest='CheckRange', type=str, nargs=2,
                        help='compare max values against (0) a confidence'
                        ' interval and (1) deviation tolerance')

    parser.add_argument('--SmoothedDev', dest='SmoothedDev', action='store_true',
                        help='used smoothed deviation values')

    parser.add_argument('--S', dest='SingleTest', action='store_true',
                        help='flag to indicate running on a specific test rather a group')

    parser.add_argument('--Debug', dest='Debug', action='store_true',
                        help='turn on debugging output')

    # parse command line args
    args = parser.parse_args()
    
    if (args.Debug):
        print "Test Dir:",args.TestDir
        
    # -------------------------------------------------------------------------------
    # Input Checking
    # -------------------------------------------------------------------------------
       
    # check for test parent directory
    if (not os.path.isdir(args.TestDir)):
        print "ERROR:",args.TestDir,"does not exist"
        sys.exit()

    if (not args.SingleTest):
        # get individual test directories in test parent and sort
        tmp_Tests = FileHelper.get_immediate_subdirectories(args.TestDir)
        tmp_Tests = sorted(tmp_Tests, key=FileHelper.numerical_sort)    
    else:
        # create list with single test
        tmp_Tests = [args.TestDir]

    # iterate over individual tests and remove those withoutis a .out file or
    # with multiple .out files
    Tests = []
    for t in tmp_Tests:
    
        # is there only one .out and .err file? if not then this test was likely
        # restarted and needs some post-processing to combine .out and .err files
        outfiles = FileHelper.get_files_with_extension(t,'.out')

        if (len(outfiles) != 1):
            print "ERROR:",len(outfiles),".out files found in",t
            print outfiles
            continue

        errfiles = FileHelper.get_files_with_extension(t,'.err')

        if (len(errfiles) != 1):
            print "ERROR:",len(errfiles),".out files found in",t
            print errfiles
            continue

        Tests.append(t)

    # -------------------------------------------------------------------------------
    # Extract Data
    # -------------------------------------------------------------------------------

    fout = open(args.Label+"_results.txt",'w')
    
    # iterate over individual test cases
    t_count = 0

    for t in Tests:

        testname = t.split('/')[-1]        
        print testname

        # get output and error file names
        outfile = FileHelper.get_files_with_extension(t,'.out')[0]
        errfile = FileHelper.get_files_with_extension(t,'.err')[0]

        # check that the run finished
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

        # initialize settings
        integrator = '-'
        method     = None
        max_nli    = '-'
        max_li     = None
        stepsize   = None

        # get run settings
        with open(outfile) as fn:

            for line in fn:
                split_line = shlex.split(line)

                if ("--timescheme" in split_line):
                    if (split_line[2] != "[arkode]"):
                        integrator = "tempest"
                        method = split_line[2][1:-1]
                
                if ("--arkode_butchertable" in split_line):
                    integrator = "arkode"
                    method = split_line[2][1:-1]

                if ("--arkode_nonliniters" in split_line):
                    max_nli = split_line[2][1:-1]
                    if (max_nli < 0):
                        max_nli = '-'

                if ("--arkode_liniters" in split_line):
                    if (max_li == None):
                        max_li = split_line[2][1:-1]

                if ("--arkode_columnsolver" in split_line):
                    max_li = '-'                       

                if ("--dt" in split_line):

                    dt = split_line[2][1:-1]

                    if (dt[-1] == 'u'):
                        stepsize = (1e-6 * float(dt[:-1]))
                    elif (dt[-1] == 's'):
                        stepsize = (float(dt[:-1]))
                    else:
                        print "ERROR: Unknown time step size units in",t
                        sys.exit()

                if ("MODEL SETUP" in line):
                    break

        if (stepsize == None):
            print "WARNING: Step size not found for",t
            continue

        if (method == None):
            print "WARNING: Method not found for",t
            continue

        # get run time
        runtime = None

        with open(errfile) as fn:

            for line in fn:
                split_line = shlex.split(line)

                if ("real" in split_line):
                    time = split_line[1]

                    j = 0
                    for c in time:
                        if (c == 'm'):
                            minutes = float(time[:j])
                            seconds = float(time[j+1:-1])
                            runtime = 60.0 * minutes + seconds
                            break
                        j += 1

        if (runtime == None):
            print "WARNING: Run time not found for",t
            continue

        # list of output files
        ncfiles = FileHelper.list_files(os.path.join(t, 'outTempest'),'out.*.nc')

        if (len(ncfiles) == 0):
            print "Warning: No output files found in",t
            continue

        ncfiles = sorted(ncfiles, key=FileHelper.numerical_sort)

        # allocate arrays assuming all tests have the same number of output files
        times   = np.empty([len(ncfiles)])
        maxvals = np.empty([len(ncfiles)])

        # read netcdf files
        ncf_count = 0
        for ncf in ncfiles:
            
            data = Dataset(ncf, mode="r")

            times[ncf_count] = data.variables['time'][...]
                
            var = data.variables[args.var][...]
            maxvals[ncf_count] = np.amax(np.abs(var))
                
            data.close()

            ncf_count += 1                
                
        # combine data and write to files
        OutData = np.column_stack((times, maxvals))
            
        fname = args.Label+"_"+testname+"_max"+args.var+".txt"
        
        np.savetxt(fname, OutData)

        # ---------------------------------------------------------------
        # Check if values are in confidence interval
        # ---------------------------------------------------------------
        
        # initialize values for output in case we skip check range
        all_passed = True          
        bad_times  = np.array([-1])

        if (args.CheckRange):

            # load confidence interval data
            rangedata = np.loadtxt(args.CheckRange[0])
            
            times  = rangedata[:,0] # output times
            mean   = rangedata[:,1] # mean value

            # confidence interval deviation
            if (not args.SmoothedDev):
                dev = rangedata[:,2] # raw deviation
            else:
                dev = rangedata[:,3] # smoothed deviation  

            # deviation tolerance
            devtol = float(args.CheckRange[1]) * np.amax(dev) 
                        
            # deviation tolerance range
            lower_devtol = mean - devtol
            upper_devtol = mean + devtol

            # zero out negative values
            idx = lower_devtol < 0.0
            lower_devtol[idx] = 0.0       

            # max allowed deviation from the mean
            testdev = np.maximum(dev, devtol)

            # do maxvals lie in combined range
            diff = np.abs(maxvals - mean)
        
            outofbounds = np.greater(diff, testdev)
            badidx = np.where(outofbounds)
        
            if np.any(outofbounds):
                all_passed = False
                maxout = np.amax(diff[np.where(outofbounds)])
            else:
                all_passed = True
                maxout = 0       
                            
            # temporary arrays of bad values from min tol check
            bad_maxvals = maxvals[badidx]
            bad_times   = times[badidx]
    
            print "AllPassed: ", all_passed
            print "OutTimes:  ", np.array_str(bad_times, max_line_width=150)

        # ---------------------------------------------------------------
        # Write Run Stats to File
        # ---------------------------------------------------------------

        if (t_count == 0):
            method_old = method

            print >> fout, '='*80
            print >> fout, method,"(",integrator,")"
            fmt = '{0:>10s} {1:^11s} {2:^10s} {3:>15s} {4:^15s}'
            print >> fout, fmt.format("Step Size", "Max NLI", "Max LI", "Runtime (s)", 
                                      "Range Check")
            print >> fout, '='*80

        elif (method != method_old):        
            method_old = method
            print >> fout, '='*80
            print >> fout, method,"(",integrator,")"
            fmt = '{0:>10s} {1:^11s} {2:^10s} {3:>15s} {4:^15s}'
            print >> fout, fmt.format("Step Size", "Max NLI", "Max LI", "Runtime", 
                                      "Range Check")
            print >> fout, '='*80

        fmt = '{0:>10.0f} {1:^11s} {2:^10s} {3:>15.3f} {4:^15}'
        print >> fout, fmt.format(stepsize, str(max_nli), str(max_li), runtime, 
                                  all_passed)
        fmt = '{0:<5s} {1:<}'            
        print >> fout, fmt.format("Err:",np.array_str(bad_times, max_line_width=75))

        # update test counter
        t_count += 1
        
    fout.close()

# ===============================================================================

if __name__ == "__main__":
    main()

# EOF

