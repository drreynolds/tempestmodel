#!/usr/bin/env python
# ===============================================================================
# Script to generate convergence and work-precision plots using data files 
# created by get_convergence_workprecision_data.py script. 
#
# Example call:
#     python plot_convergence_workprecision_data.py \
#     HEVI-A_Convergence_and_Work-Precision_Plots \
#     ~/data/ARS343_errors.txt \
#     ~/data/ARS232_errors.txt \
#     --SetMarkers 0 1 \
#     --SetColors 0 1 \
#     --Legend ARS343 ARS232 \
#     --Order 2 --OrderHShift 0 --OrderVShift 10 --OrderRange 0 -5
#
# D.J. Gardner @ LLNL, Sept 2016
# ===============================================================================

def main():

    import argparse
    import os, sys
    import shlex
    
    import numpy as np
    import matplotlib.pyplot as plt  # plotting functions
    import matplotlib.cm     as cm   # colormap for graphics

    parser = argparse.ArgumentParser(
        description='Generate convergence plots')
    
    parser.add_argument('Label', type=str,
                        help='label used in plot title and figure file name')
        
    parser.add_argument('DataFiles', type=str, nargs='+',
                        help='data files to plot')

    parser.add_argument('--Components', dest='Components',
                        action='store_true',
                        help='data file also contains errors for individual errors')

    parser.add_argument('--PrecisionOff', dest='PrecisionPlot',
                        action='store_false',
                        help='do not plot work-precision')

    parser.add_argument('--Norm', dest='Norm', choices=['RMS','L2','L1','Max'],
                        default='RMS',
                        help='specify which norm was used in DataFiles')

    parser.add_argument('--Legend', type=str, nargs='+', dest='Legend',
                        help='set plot labels in legend')

    parser.add_argument('--LegendLoc', dest='LegendLoc', type=float, nargs=2,
                        default=[0.58, 0.13],
                        help='adjust legend location')

    parser.add_argument('--Order', dest='Order', type=int, nargs='+',
                        help='add reference order slope to convergence plot')

    parser.add_argument('--OrderHShift', dest='OrderHShift', type=float, nargs='+',
                        help='adjust horizontal position of reference line')

    parser.add_argument('--OrderVShift', dest='OrderVShift', type=float, nargs='+',
                        help='adjust vertical position of reference line')

    parser.add_argument('--OrderRange', dest='OrderRange', type=int, nargs='+',
                        help='limit step sizes used to plot reference line')

    parser.add_argument('--SetColors', type=int, nargs='+', dest='SetColors',
                        help='set line colors')
    
    parser.add_argument('--SetLines', type=int, nargs='+', dest='SetLines',
                        help='set line styles')
    
    parser.add_argument('--SetMarkers', type=int, nargs='+', dest='SetMarkers',
                        help='set marker styles')
    
    parser.add_argument('--Show', dest='Show', action='store_true', 
                        help='show plot, do not write to file')
    
    parser.add_argument('--Debug', dest='Debug', action='store_true',
                        help='turn on debugging output')

    # parse command line args
    args = parser.parse_args()

    # -------------------------------------------------------------------------------
    # Plot Settings
    # -------------------------------------------------------------------------------

    # number of outputs
    nplots = len(args.DataFiles) 

    # line color map, qualitative colors from colorbrewer2.org
    if (args.SetColors):    

        if (len(args.SetColors) != nplots):
            print "ERROR: len(SetColors) != len(outfiles)"
            sys.exit()

        totalcolors = max(args.SetColors)+1

    else:
        totalcolors = nplots

        args.SetColors = range(nplots)

    if (totalcolors < 10):
        colors = ['#000000','#e41a1c','#377eb8','#4daf4a','#984ea3',
                  '#ff7f00','#a65628','#f781bf','#999999']
    elif (totalcolors < 14):
        colors = ['#000000','#a6cee3','#1f78b4','#b2df8a','#33a02c',
                  '#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6',
                  '#6a3d9a','#ffff99','#b15928']
    else:
        colors = cm.Set1(np.linspace(0, 1, totalcolors))

    Cvalue = []
    for i in args.SetColors:
        Cvalue.append(colors[i])

    # line styles
    LineStyles = ['-','--','-.',':']
    if (args.SetLines):
        
        if (len(args.SetLines) != nplots):
            print "ERROR: len(SetLines) != len(outfiles)"
            sys.exit()

        if (max(args.SetLines)-1 > len(LineStyles)):
            print "Only "+len(LineStyles)+" line styles are available"
            sys.exit()

        Lstyle = []
        for i in args.SetLines:
            Lstyle.append(LineStyles[i])

    else:
        Lstyle = ['-'] * nplots

    # Line Markers
    MarkerStyles = ['.','o',
                    'v','^','<','>',
                    'd','D','s','p','h','H','8',
                    '*',
                    '+','x','1','2','3','4','|','_']
    if (args.SetMarkers):

        if (len(args.SetMarkers) != nplots):
            print "ERROR: len(SetMarkers) != len(outfiles)"
            sys.exit()

        if (max(args.SetMarkers)-1 > len(MarkerStyles)):
            print "ERROR: Only",len(MarkerStyles),"marker styles are available"
            sys.exit()

        Mstyle = []
        for i in args.SetMarkers:
            Mstyle.append(MarkerStyles[i])

    else:
        Mstyle = [None] * nplots
    
    # legend
    if (args.Legend):
        if (len(args.Legend) != len(args.DataFiles)):
            print "ERROR: len(Legend) != len(Tests)"
            sys.exit()

        Leg = []
        for i in range(0,len(args.DataFiles)):
            if (args.Legend[i] == 'None'):
                Leg.append(None)
            else:
                Leg.append(args.Legend[i])
    else:
        Leg = [None] * nplots

    # line width
    Lwidth = 1.0

    # -------------------------------------------------------------------------------
    # Make Plots
    # -------------------------------------------------------------------------------

    # intialize counter
    counter = 0

    # create plot figure
    if (args.PrecisionPlot):
        plt.figure("Convergence and Work-Precision")
        fig1, (ax1, ax2) = plt.subplots(1, 2, sharey=True)

    else:
        plt.figure("Convergence")
        fig1, ax1 = plt.subplots()

    # loop over data files, assume one plot per data file
    for dfile in args.DataFiles:

        # check if data file exists
        if (not os.path.isfile(dfile)):
            print "Warning:",dfile,"not found"
            continue

        # extract step size, errors, and run times
        RunData = np.loadtxt(dfile)
                
        stepsizes = RunData[:,0]
        errors    = RunData[:,1]
        runtimes  = RunData[:,2]

        if (args.Components):
            errorsU = RunData[:,3]
            errorsV = RunData[:,4]
            errorsW = RunData[:,5]
            errorsR = RunData[:,6]
            errorsT = RunData[:,7]

        # assumes smallest to largest dt in data file
        if (counter == 0):
            maxlen   = np.shape(stepsizes)[0]
            ordsteps = stepsizes
            orderr   = errors

            minstep = stepsizes[0]
            maxstep = stepsizes[-1]
            maxerr  = np.amax(errors)
            minerr  = np.amin(errors)
        else:
            if (np.shape(stepsizes)[0] > maxlen):
                maxlen   = np.shape(stepsizes)[0]
                ordsteps = stepsizes
                orderr   = errors

            if (stepsizes[0] < minstep):
                minstep = stepsizes[0] 

            if (stepsizes[-1] > maxstep):
                maxstep = stepsizes[-1]

            if (np.amax(errors) > maxerr):
                maxerr = np.amax(errors)

            if (np.amin(errors) < minerr):
                minerr = np.amin(errors)

        # plot convergence
        ax1.loglog(stepsizes, errors,
                   linestyle=Lstyle[counter], 
                   color=Cvalue[counter], 
                   marker=Mstyle[counter], 
                   linewidth = Lwidth, 
                   label=Leg[counter])

        if (args.Components):
            ax1.loglog(stepsizes, errorsU,
                       linestyle=Lstyle[counter],
                       color=Cvalue[counter],
                       marker=MarkerStyles[1],
                       linewidth = Lwidth,
                       label=Leg[counter]+' U')

            ax1.loglog(stepsizes, errorsV,
                       linestyle=Lstyle[counter],
                       color=Cvalue[counter],
                       marker=MarkerStyles[6],
                       linewidth = Lwidth,
                       label=Leg[counter]+' V')

            ax1.loglog(stepsizes, errorsW,
                       linestyle=Lstyle[counter],
                       color=Cvalue[counter],
                       marker=MarkerStyles[8],
                       linewidth = Lwidth,
                       label=Leg[counter]+' W')

            ax1.loglog(stepsizes, errorsR,
                       linestyle=Lstyle[counter],
                       color=Cvalue[counter],
                       marker=MarkerStyles[9],
                       linewidth = Lwidth,
                       label=Leg[counter]+' R')

            ax1.loglog(stepsizes, errorsT,
                       linestyle=Lstyle[counter],
                       color=Cvalue[counter],
                       marker=MarkerStyles[10],
                       linewidth = Lwidth,
                       label=Leg[counter]+' T')

        if (args.PrecisionPlot):

            # plot work-precision
            ax2.loglog(runtimes, errors,
                       linestyle=Lstyle[counter],
                       color=Cvalue[counter],
                       marker=Mstyle[counter],
                       linewidth = Lwidth,
                       label=Leg[counter])

        # update plot counter
        counter += 1

    # -------------------------------------------------------------------------------
    # Add Reference Line to Plots
    # -------------------------------------------------------------------------------
    if (args.Order):

        if (args.OrderHShift):
            if (len(args.OrderHShift) != len(args.Order)):
                print "ERROR: --OrderHShift and --Order are different lengths"
                sys.exit()
        else:
            args.OrderHShift = len(args.Order)*[0.0]

        if (args.OrderVShift):
            if (len(args.OrderVShift) != len(args.Order)):
                print "ERROR: --OrderVShift and --Order are different lengths"
                sys.exit()
        else:
            args.OrderVShift = len(args.Order)*[1.0]

        if (args.OrderRange):
            if (len(args.OrderRange) != 2*len(args.Order)):
                print "ERROR: --OrderRange is not twice the length of --Order"
                sys.exit()

        LineStyles = ['--','-.',':']

        # initialize counter
        counter = 0
        
        for p in range(0,len(args.Order)):

            if (ordsteps[-1] != maxstep):
                np.append(ordsteps, [maxstep])

            if (ordsteps[0] != minstep):
                ordsteps = np.insert(ordsteps, 0, minstep)

            X = ordsteps
            Y = minerr * (X / X[0])**args.Order[p]

            # horizontal shift
            if (args.OrderHShift[p] < 0):
                # left shift
                X = X/np.log10(10 + abs(args.OrderHShift))
            elif (args.OrderHShift[p] > 0):
                # right shift
                X = X*np.log10(10 + args.OrderHShift)
                
            # vertical shift
            Y = args.OrderVShift[p] * Y

            # truncate order line
            if (args.OrderRange):
                a = args.OrderRange[2*p]
                b = np.shape(ordsteps)[0] + args.OrderRange[2*p+1] + 1
                X = X[a:b]
                Y = Y[a:b]

            # plot order line
            ax1.loglog(X, Y,
                       linestyle=LineStyles[p],
                       color='black', 
                       linewidth = Lwidth,
                       label='Order '+str(args.Order[p]))
            
            # update counter
            counter += 1

    # -------------------------------------------------------------------------------
    # Finalize Plots
    # -------------------------------------------------------------------------------

    title = fig1.suptitle(args.Label.replace("_"," "))

    ax1.set_xlabel('Step Size (s)')
    ax1.set_ylabel(args.Norm+' Error ')
    ax1.grid(True)

    handles, labels = ax1.get_legend_handles_labels()

    if (args.PrecisionPlot):
        ax2.set_xlabel('Run Time (s)')
        ax2.grid(True)

    # shift upper right corner (0 = no shift)
    x = args.LegendLoc[0] 
    y = args.LegendLoc[1]

    lgd = fig1.legend( handles, labels, loc = 'center', 
                       bbox_to_anchor = (x, y, 1, 1),
                       bbox_transform = plt.gcf().transFigure )
    
    if (args.Show):
        # show plot on screen
        plt.show()
    else:
        # save plot to file
        figname='State_Convergence_and_Work_Precision_'+args.Label+'.pdf'
        fig1.savefig(figname, dpi='300', format='pdf', 
                    bbox_extra_artists=(lgd,title,), bbox_inches='tight')

    plt.close('all')

# ===============================================================================

if __name__ == "__main__":
    main()

# EOF
