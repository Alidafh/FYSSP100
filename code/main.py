import ROOT as ro
import ctypes
import statistics as stat
import numpy as np
import analysis_functions as analysis
import sys, os
ro.gROOT.SetBatch(ro.kTRUE)


################################################################################
# Create mass plot
analysis.MassPlot(20)
#analysis.MassPlot(10)

# a), b) Mass window that optimizes the expected and observed significance and figure.
max_exp1, bestwidth_exp1 = analysis.Significance_Optimization(LumiScale=1, plot="plot")

# c) Mass window that optimizes the expected and observed significance for 5 times Luminosity
max_exp2, bestwidth_exp2 = analysis.Significance_Optimization(LumiScale=5, plot="plot")

# d) Luminosity for discovery
print("Finding the lowest Luminosity required for discovery")
best_width = 0
best_lumi = 0
LumiScales = np.linspace(0, 10, 50)
for i in range(len(LumiScales)):
    sign, width = analysis.Significance_Optimization(LumiScales[i])
    if sign >= 5:
        print("      Expected significance: {:.2f} at width {:.2f} and luminosity {:.2f}\n".format(sign, width, LumiScales[i]))
        best_width = width
        best_lumi = LumiScales[i]
        break


################################################################################
# Part 2

# Preform sideband fit to find alpha and uncertainties
#bestalpha, sigmalow, sigmaup = analysis.sideband_fit()

# Takes time to do sideband_fit, so nubers are entered here for the time being
bestalpha = 1.11
sigmalow = 0.06
sigmaup = 0.07

hist_bkg = analysis.GetMassDistribution(2)
low = hist_bkg.FindBin(125.-0.5*(bestwidth_exp1))
high = hist_bkg.FindBin(125. + 0.5*(bestwidth_exp1))
nbkg = hist_bkg.Integral(low, high)

print("\nUsing optimal signal region width {:.2f} and bkg. scale factor {:.2f}:".format(bestwidth_exp1, bestalpha))
print("      unscaled bkg. events: {:.3f}".format(nbkg))
print("      scaled bkg. events:   {:.3f} -{:.3f} +{:.3f}".format(nbkg*bestalpha, sigmalow*nbkg, sigmaup*nbkg))

# Find number of signal events in the signal region
hist_signal = analysis.GetMassDistribution(0)
low = hist_signal.FindBin(125.-0.5*(bestwidth_exp1))
high = hist_signal.FindBin(125. + 0.5*(bestwidth_exp1))
nsig = hist_signal.Integral(low,high)
print("      unscaled sig. events: {:.3f}\n".format(nsig))

# Find expected significance
#Ntoys = 1e6
Ntoys = 100
print("Calculating the expected significance using {} toys:".format(Ntoys))
pvalue = analysis.ExpectedSignificance_ToyMC(nbkg, nsig, sigmalow*nbkg, Ntoys)

################################################################################

print("\nCalculating test-statisitic for real data:")
hist_data = analysis.GetMassDistribution(3)
hist_signal = analysis.GetMassDistribution(0)
hist_bkg = analysis.GetMassDistribution(2)

t_data = analysis.Get_TestStatistic(hist_data, hist_bkg, hist_signal)
print("      t = {:.3f}\n".format(t_data))

###############################################################################
# Execute toy generator before moving on

def calc_mc():
    #check if the MC has already been calculated:
    path = "output/histograms/test_statistic_distribution.root"

    if os.path.isfile(path) == True:
        answer = input("Do you want to use existing distributions for test-statistics? [y]")
        if answer != "y":
            analysis.Generate_toys()

    if os.path.isfile(path) == False:
        print("Need to calculate distributions for test-statistics")
        answer = input("This may take some time. Continue? [y]")
        if answer == "y":
            analysis.Generate_toys()
        else:
            print("Quitting..")
            quit()

if int(sys.argv[1]) == 1:
    calc_mc()

###############################################################################

# Get the test-statistic distributions from file
histograms = ["test_statistic_bkg", "test_statistic_sb", "test_statistic_data"]
infile = ro.TFile("output/histograms/test_statistic_distribution.root", "READ")
hist_distribution_b = infile.Get(histograms[0]).Clone("h_tdist_b")
hist_distribution_sb = infile.Get(histograms[1]).Clone("h_tdist_sb")
#hist_distribution_data = infile.Get(histograms[2]).Clone("h_tdist_data")
infile.Close()

analysis.analyze_distributions(hist_distribution_b, hist_distribution_sb, t_data, "plot")

###############################################################################
