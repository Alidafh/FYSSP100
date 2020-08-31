import ROOT as ro
import ctypes
import statistics as stat
import numpy as np
import analysis_functions as analysis
import sys, os
import matplotlib.pyplot as plt
ro.gROOT.SetBatch(ro.kTRUE)


################################################################################
# Create mass plot
analysis.MassPlot(20)
#analysis.MassPlot(10)

# Mass window that optimizes the expected and observed significance
exp_significance, signal_width = analysis.Significance_Optimization(LumiScale=1, plot="plot")

# for 5 times the luminosity
exp_significance_lumi5, width_lumi5 = analysis.Significance_Optimization(LumiScale=5, plot="plot")


print("Finding the lowest Luminosity required for discovery")
best_width = 0
best_lumi = 0
LumiScales = np.linspace(0, 10, 50)
z_array = np.zeros(50)
for i in range(len(LumiScales)):
	# Loop over luminosity scale factors
    significance, width = analysis.Significance_Optimization(LumiScales[i])
    z_array[i] = significance
    if significance >= 5:
        print("      Expected significance: {:.2f} at width {:.2f} and luminosity {:.2f}\n".format(significance, width, LumiScales[i]))
        best_width = width
        best_lumi = LumiScales[i]
        break

"""
plt.figure(1)
plt.plot(LumiScales, z_array)
plt.xlabel("Luminosity scale factro")
plt.ylabel("Significance, Z")
plt.show()
quit()
"""

# Plot this.
#exp_significance_lumiX, width_lumiX = analysis.Significance_Optimization(best_lumi, "plot")

################################################################################
# Preform sideband fit to find bkg-scale factor
bestalpha, sigmalow, sigmaup = analysis.sideband_fit()

# Plot the scaled sideband region for comparison
analysis.plotMassSideband(20, bestalpha, 0.10)

# Count number of expected signal and background events in signal region
hist_bkg = analysis.GetMassDistribution(2)
hist_signal = analysis.GetMassDistribution(0)
hist_data = analysis.GetMassDistribution(3)

nbkg = analysis.count_events(hist_bkg, signal_width)        # events under H0
nsignal = analysis.count_events(hist_signal, signal_width)  # signal events
ndata = analysis.count_events(hist_data, signal_width)      # observed events

print("\nIn the signal region of width {:.2f} GeV:".format(signal_width))
print("      expected bkg. events: {:.3f}".format(nbkg))
print("      expected sig. events: {:.3f}".format(nsignal))
print("      observed events:      {:.3f}".format(ndata))

nbkg_scaled = nbkg*bestalpha        # scaled bkg events in signal region
d_nbkg_scaled = sigmaup*nbkg        # Unc. on the scaled bkg estimate

print("      scaled bkg. events:   {:.3f} -{:.3f} +{:.3f}\n".format(nbkg_scaled, sigmalow*nbkg, sigmaup*nbkg))

# Calculate the expected and observed significance with new background estimate

if int(sys.argv[1]) == 0: Ntoys = 100
if int(sys.argv[1]) != 0: Ntoys = 1e6

d_nbkg = 0.10                       # Error on the bkg estimate
d_nbkg_scaled = sigmaup*nbkg        # Error on the scaled bkg estimate

pvalue_exp = analysis.ExpectedSignificance_ToyMC(nbkg_scaled, nsignal, d_nbkg_scaled, Ntoys)
pvalue_obs = stat.IntegratePoissonFromRight(ndata, nbkg_scaled)

sigma_exp = ro.Math.gaussian_quantile_c(pvalue_exp, 1)
sigma_obs = ro.Math.gaussian_quantile_c(pvalue_obs, 1)

print("      expected p-value = {:.3f}, sigificance = {:.3f}".format(pvalue_exp, sigma_exp))
print("      observed p-value = {:.3f}, sigificance = {:.3f}".format(pvalue_obs, sigma_obs))

print("------------------------------------------------------------------------")
print("\nUsing the unscaled MC estimate and a relative error of 10 percent:\n")

pvalue_exp_unsc = analysis.ExpectedSignificance_ToyMC(nbkg, nsignal, d_nbkg, Ntoys)
pvalue_obs_unsc = stat.IntegratePoissonFromRight(ndata, nbkg)

sigma_exp_unsc = ro.Math.gaussian_quantile_c(pvalue_exp_unsc, 1)
sigma_obs_unsc = ro.Math.gaussian_quantile_c(pvalue_obs_unsc, 1)

print("      expected p-value = {:.3f}, sigificance = {:.3f}".format(pvalue_exp_unsc, sigma_exp_unsc))
print("      observed p-value = {:.3f}, sigificance = {:.3f}".format(pvalue_obs_unsc, sigma_obs_unsc))
print("------------------------------------------------------------------------")

################################################################################

print("\nCalculating observed test-statisitic:")
hist_data = analysis.GetMassDistribution(3)
hist_signal = analysis.GetMassDistribution(0)
hist_bkg = analysis.GetMassDistribution(2)

t_data = analysis.Get_TestStatistic(hist_data, hist_bkg, hist_signal)
print("      t = {:.3f}\n".format(t_data))

###############################################################################
# Execute toy generator before moving on, if this is already calculated there
# is no need to do it again.
#(ntoys, sf_bkg= 1, sf_sig=1, rebinN=1)

if int(sys.argv[1]) != 0:
    analysis.Generate_toys(10000, 1, 1, 1)

###############################################################################

# Get the test-statistic distributions from file
histograms = ["test_statistic_bkg", "test_statistic_sb"]
infile = ro.TFile("output/histograms/test_statistic_distribution_sfb1_sfs1_toys10000_bin1.root", "READ")
hist_distribution_b = infile.Get(histograms[0]).Clone("h_tdist_b")
hist_distribution_sb = infile.Get(histograms[1]).Clone("h_tdist_sb")
infile.Close()

analysis.analyze_distributions(hist_distribution_b, hist_distribution_sb, t_data, "plot")


###############################################################################
# different luminosity scale factor

def new_scales(sfb, sfs, rebinN, ntoys):
    print("-----------------------------------------------------")
    print("Using luminosity scale factor:", sfb)

    hist_dat_scaled = analysis.GetMassDistribution(3, sfs, rebinN)
    hist_sig_scaled = analysis.GetMassDistribution(0, sfs, rebinN)
    hist_bkg_scaled = analysis.GetMassDistribution(2, sfb, rebinN)

    t_obs_scaled = analysis.Get_TestStatistic(hist_dat_scaled, hist_bkg_scaled, hist_sig_scaled)
    print("      t = {:.3f}\n".format(t_obs_scaled))

    analysis.Generate_toys(ntoys, sfb, sfs, rebinN)

    histograms = ["test_statistic_bkg", "test_statistic_sb"]
    infile1 = ro.TFile("output/histograms/test_statistic_distribution_sfb{}_sfs{}_toys{}_bin{}.root".format(sfb, sfs, ntoys,rebinN), "READ")
    hist_dist_b_scaled = infile1.Get(histograms[0]).Clone("h_tdist_b_scaled")
    hist_dist_sb_scaled = infile1.Get(histograms[1]).Clone("h_tdist_sb_scaled")
    infile1.Close()

    analysis.analyze_distributions(hist_dist_b_scaled, hist_dist_sb_scaled, t_obs_scaled, "plot")

    print(" ")

#new_scales(1.5, 1.5, 10, 10000)
#new_scales(2.0, 2.0, 10, 10000)
#new_scales(5.0, 5.0, 10, 10000)

########################################################################

hist_loglik, sf_sig, sf_bkg  = analysis.muFit(1, 100)
analysis.plot_mu(hist_loglik, sf_sig, sf_bkg)


# Last minute function, values in the function are found running new_scales above
# and manually entered in the extrapolate function in analysis_functions.py...
c11 = analysis.extrapolate("clb")
c22 = analysis.extrapolate("clsb")
