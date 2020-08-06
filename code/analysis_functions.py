import ROOT as ro
import numpy as np
from statistics import IntegratePoisson, IntegrateFromRight
import statistics as stat

def Quiet(func,level = ro.kInfo + 1):
    def qfunc(*args,**kwargs):
        oldlevel = ro.gErrorIgnoreLevel
        ro.gErrorIgnoreLevel = level
        try:
            return func(*args,**kwargs)
        finally:
            ro.gErrorIgnoreLevel = oldlevel
    return qfunc

def text(x, y, text_string, size=0.022):
    t = ro.TLatex()
    t.SetNDC(1)
    t.SetTextAlign(13);
    t.SetTextSize(size);
    t.SetTextFont(42);
    t.DrawLatex(x, y, text_string)
    #t.Draw("same")

def GetMassDistribution(type, scaleFactor = 1):
    """
    Returns the scaled histogram from file.
        type: 0 = Higgs125, 1 = Higgs200, 2 = ZZ SM background, 3 = data
        scaleFactor: 1 if no scaling is neccesary
    """
    ro.TH1.AddDirectory(ro.kFALSE)
    histograms = ["h_m4l_Higgs125_fake", "h_m4l_Higgs200_fake", "h_m4l_ZZ_fake", "h_m4l_data_fake"]
    infile = ro.TFile("input/Histograms_fake.root", "READ")
    h_mass = infile.Get(histograms[type]).Clone("h_mass")
    h_mass.SetTitle("; m_{4l} GeV"+";Events/ {:.1f}".format(h_mass.GetBinWidth(1)))
    nbins = h_mass.GetNbinsX()
    for i in range(nbins):
        bin_cont = h_mass.GetBinContent(i)
        h_mass.SetBinContent(i,scaleFactor*bin_cont)
    infile.Close()
    return h_mass

def MassPlot(rebinN):
	"""
	plots the mass distribution
	rebin: rebin histograms for better plotting
	"""
	print("Plotting the invariant mass distribution with rebin = {:.2f}".format(rebinN))
	hist_signal = GetMassDistribution(0)
	hist_bkg = GetMassDistribution(2)
	hist_data = GetMassDistribution(3)

	# Rebin histograms
	hist_signal.Rebin(rebinN)
	hist_bkg.Rebin(rebinN)
	hist_data.Rebin(rebinN)

	# make signal+bkg histogram
	hist_sb = hist_bkg.Clone("hist_sb")
	hist_sb.Add(hist_signal)

	hist_sb.SetFillColor(7)
	hist_sb.SetAxisRange(0, 22, "Y")
	hist_sb.SetAxisRange(0,400, "X")
	hist_bkg.SetFillColor(2)

	legend = ro.TLegend(0.75, 0.75, 0.90, 0.85)
	legend.AddEntry(hist_sb, "Higgs signal", "f")
	legend.AddEntry(hist_bkg, "Background", "f")
	legend.AddEntry(hist_data, "Data", "p")

	sr1 = ro.TLine(150, 0, 150, 20)
	sr1.SetLineStyle(7)

	asr1 = ro.TArrow(150, 12, 170, 12, 0.01, "-|>")

	c = ro.TCanvas("c", "c", 1000, 600)
	hist_sb.Draw("hist")
	hist_sb.GetYaxis().SetTitle("Events/ {:.1f} GeV".format(hist_sb.GetBinWidth(1)))
	hist_bkg.Draw("same")
	hist_data.Draw("e same")
	sr1.Draw()
	asr1.Draw()
	text(0.38, 0.65, "Sideband", size=0.022)
	text(0.385, 0.63, "Region", size=0.022)
	legend.Draw()

	c.Update()
	Quiet(c.SaveAs)("output/figures/mass_plot_125_{}.png".format(rebinN))
	print("      figure saved in: output/figures/mass_plot_125_{}.png\n".format(rebinN))

def Significance_Optimization(LumiScale, plot="no"):
	"""
	Finds the mass window that gives the highest expected significance
		LumiScale - If you want to scale the luminosity,
		plot="no" - If a figure is desired, do plot="plot"
	"""
	# Get histograms
	if plot != "no": print("For a Luminosity scale: {:.1f}".format(LumiScale))
	hist_signal = GetMassDistribution(0, LumiScale)
	hist_bkg = GetMassDistribution(2, LumiScale)
	hist_data = GetMassDistribution(3, LumiScale)

	hist_window = ro.TH1D("h_masswindow", ";m_{4l} GeV; ",250,0.,25.); # make a mass window - full width between 0 and 25 GeV
	hist_window_exp = ro.TH1D("h_masswindow_expected",";width [GeV]; Significance Z",250,0.,25.);
	hist_window_obs = ro.TH1D("h_masswindow_observed",";width [GeV]; Significance Z",250,0.,25.);

	nbins = hist_window.GetNbinsX()
	for i in range(nbins):
		window_width = hist_window.GetBinCenter(i)
		#print("Trying mass window: {:.2f}\n".format(window_width))
		x1 = 125-0.5*window_width
		x2 = 125+0.5*window_width
		ndata = hist_data.Integral(hist_data.FindBin(x1), hist_data.FindBin(x2))
		nbkg  = hist_bkg.Integral(hist_bkg.FindBin(x1), hist_bkg.FindBin(x2))
		nsig  = hist_signal.Integral(hist_signal.FindBin(x1), hist_signal.FindBin(x2))

		if( (nbkg+nsig)<1 ): continue

		#Calculating the expected significance
		exp_p = 1 - IntegratePoisson(nsig+nbkg, nbkg)
		exp_z = ro.Math.normal_quantile_c(exp_p, 1)
		hist_window_exp.SetBinContent(i, exp_z)

		# Calculating observed significance
		obs_p = 1 - IntegratePoisson(ndata, nbkg)
		obs_z = ro.Math.normal_quantile_c(obs_p, 1)
		hist_window_obs.SetBinContent(i, obs_z)

	if plot=="plot":
		# plot the results and save figure
		c = ro.TCanvas("c{}".format(LumiScale), "c", 1000, 600)
		#legend = ro.TLegend(0.75, 0.75, 0.90, 0.85)
		legend = ro.TLegend(0.75, 0.79, 0.93, 0.93)
		legend.AddEntry(hist_window_exp, "Expected", "l")
		hist_window_obs.SetLineColor(1)
		hist_window_exp.SetLineColor(4)

		if( abs(LumiScale-1.00) < 0.01 ):
		    hist_window_obs.Draw("L same")
		    legend.AddEntry(hist_window_obs, "Observed", "l")

		hist_window_exp.Draw("L same")
		legend.Draw()
		c.Update()
		Quiet(c.SaveAs)("output/figures/expected_mass_window_{}.png".format(LumiScale))
		print("      figure saved in: output/figures/expected_mass_window_{}.png".format(LumiScale))

	# Find the mass window that maximises the significance:
	max_exp = hist_window_exp.GetMaximum()
	bin_exp = hist_window_exp.GetMaximumBin()
	bestwidth_exp = hist_window_exp.GetXaxis().GetBinCenter(bin_exp)

	max_obs = hist_window_obs.GetMaximum()
	bin_obs = hist_window_obs.GetMaximumBin()
	bestwidth_obs = hist_window_obs.GetXaxis().GetBinCenter(bin_obs)

	if plot != "no" :print("      maximum expected significance: {:.2f} at mass window {:.2f}".format(max_exp, bestwidth_exp))
	if plot != "no" :print("      maximum observed significance: {:.2f} at mass window {:.2f}\n".format(max_obs, bestwidth_obs))

	return max_exp, bestwidth_exp

def sideband_fit(rebinN=1, LumiScale=1):
    #Sideband fit
    hist_bkg = GetMassDistribution(2, LumiScale)
    hist_data = GetMassDistribution(3, LumiScale)

    hist_bkg.Rebin(rebinN)
    hist_data.Rebin(rebinN)

    nbins = hist_data.GetNbinsX()
    scale_bkg = ro.TH1F("scale_bkg", ";#alpha; ln L", nbins, 0.1, 3.1)
    hist_likelihood = ro.TH1F("likelihood", ";#alpha; ln L", nbins, 0.1, 3.1)

    print("Preforming sideband fit:")
    for i in range(1, nbins+1):
        # Loop over scale factors
        alpha = scale_bkg.GetBinCenter(i)
        logL = 0.
        for j in range(1, nbins+1):
            mass = hist_data.GetBinCenter(j)
            n = hist_data.GetBinContent(j)
            if (mass >= 150. and mass <= 400.):
                func = alpha*hist_bkg.GetBinContent(j)
                logL += np.log(ro.TMath.Poisson(n,func))

        hist_likelihood.SetBinContent(i, logL)

    # Find the maximum likelihood value
    maximum_bin = hist_likelihood.GetMaximumBin()
    maximum_value = hist_likelihood.GetBinContent(maximum_bin)
    bestalpha = hist_likelihood.GetBinCenter(maximum_bin)

    # rescale the likelihood function
    h_likelihood_rescaled = hist_likelihood.Clone("h_likelihood_rescaled")
    h_likelihood_rescaled.Reset()
    h_likelihood_rescaled.SetTitle(";#alpha; ln L - ln L_{max};")

    for i in range(1, nbins+1):
        content = hist_likelihood.GetBinContent(i) - maximum_value
        h_likelihood_rescaled.SetBinContent(i, content)

    low_bin = h_likelihood_rescaled.FindFirstBinAbove(-0.5)-1
    up_bin = h_likelihood_rescaled.FindLastBinAbove(-0.5)+1
    low_value = hist_likelihood.GetBinCenter(low_bin)
    up_value = hist_likelihood.GetBinCenter(up_bin)

    low = bestalpha - low_value
    up = up_value - bestalpha

    print("      Optimal bkg. scale factor: {:.2f}, with unc.: -{:.2f}, +{:.2f} ".format(bestalpha, low, up))

    if rebinN==1:
        c = ro.TCanvas("c{}".format(rebinN), "c", 1000, 600)
        c.cd()
        h_likelihood_rescaled.SetAxisRange(-1, 0, "Y")
        h_likelihood_rescaled.SetAxisRange(1, 1.25, "X")
        h_likelihood_rescaled.Draw("L")

        line = ro.TLine(bestalpha, 0, bestalpha, -1)
        line.SetLineStyle(7); line.SetLineColor(ro.kGray+2)
        line.Draw()

        line2 = ro.TLine(low_value, -0.5, low_value, -1)
        line2.SetLineStyle(7); line2.SetLineColor(ro.kGray+2)
        line2.Draw()

        line3 = ro.TLine(up_value, -0.5, up_value, -1)
        line3.SetLineStyle(7); line3.SetLineColor(ro.kGray+2)
        line3.Draw()

    # arrows showing the error
        arrow1 = ro.TArrow(low_value, -0.5, bestalpha, -0.5, 0.01,"<|>")
        arrow2 = ro.TArrow(bestalpha, -0.5, up_value, -0.5, 0.01,"<|>")
        text(0.38, 0.6, "#Delta#hat#alpha_{-}", size=0.042)
        text(0.56, 0.6, "#Delta#hat#alpha_{+}", size=0.042)
        arrow1.Draw()
        arrow2.Draw()

        c.Update()
        Quiet(c.SaveAs)("output/figures/delta_log_likelihood.png".format(rebinN))
        print("      figure saved in: output/figures/delta_log_likelihood.png".format(rebinN))

    return bestalpha, low, up

def ExpectedSignificance_ToyMC(b, s, db, Ntoys):
    ro.gROOT.Clear()
    ro.gROOT.Delete()

    # Histograms for number of events
    h_Nbkg = ro.TH1D("h_Nbkg", "Number of background events", 500, -0.5, 499.5);
    h_Nsig_plus_bgr = ro.TH1D("h_Nsig_plus_bgr", "Number of signal + background events", 500, -0.5, 499.5);

    random = ro.TRandom3(0)

    for i in range(int(Ntoys)):
        mean_b = random.Gaus(b, db)
        mean_splusb = random.Gaus(b, db) + s

        if mean_b >= 0 and mean_splusb >=0:
            b_toy = random.Poisson(mean_b)
            splusb_toy = random.Poisson(mean_splusb)

            h_Nbkg.Fill(b_toy)
            h_Nsig_plus_bgr.Fill(splusb_toy)

        else:
            print("\n  ExpectedSignificance_ToyMC::ERROR\n")
            print("      The mean of b of s+b is smaller than 0 -> disregarding toy\n")
            print("      Note ; could also have used log_normal distribution rather than truncated Gaussian\n")

    nEvents_median = s+b
    pvalue = IntegrateFromRight(h_Nbkg, nEvents_median)
    sigma = ro.Math.gaussian_quantile_c(pvalue,1)

    print("      p-value = {:.3f}, sigificance = {:.3f}".format(pvalue, sigma))
    return pvalue

def count_events(nr):
    hist_sig = GetMassDistribution(nr)
    w = 7.15
    low = hist_sig.FindBin(125-0.5*w)
    up = hist_sig.FindBin(125+0.5*w)
    mean = hist_sig.Integral(low, up)
    print("nuber of events", mean)
    return mean

def Get_TestStatistic(h_mass_dataset, h_template_bkg, h_template_sig):
    #   Computes the likelihood ratio test-statistic for a dataset:
	#	(h_mass_dataset)
    #   from expected distributions for background and signal:
	#	(h_template_bkg, h_template_sig)
	#
    #   Returns the test statistic: X

    n = h_mass_dataset.GetNbinsX()

    loglikelihood_sb = 0
    loglikelihood_b = 0
    for i in range(n):
        data = h_mass_dataset.GetBinContent(i)
        signal = h_template_sig.GetBinContent(i)
        bkg = h_template_bkg.GetBinContent(i)
        loglikelihood_sb += np.log(ro.TMath.Poisson(data, signal+bkg))
        loglikelihood_b += np.log(ro.TMath.Poisson(data, bkg))

    X = -2*loglikelihood_sb - (-2*loglikelihood_b)
    #print(X)
    return X

def GenerateToyDataSet(h_mass_template):
    #   Generates toy datasets by drawing a random poisson number in each bin
    #   using the bin content of h_mass_template as central value.
    #   Returns the toy dataset: h_mass_toydataset

    n = h_mass_template.GetNbinsX()
    random = ro.TRandom3(0)

    h_mass_toydataset = h_mass_template.Clone("h_mass_toydataset")
    h_mass_toydataset.Reset()

    for i in range(n):
        central_value = h_mass_template.GetBinContent(i)
        h_mass_toydataset.SetBinContent(i, random.Poisson(central_value))

    return h_mass_toydataset

def Get_TestStatistic_Distribution(hist_template, ntoys, mode):
    # makes test-statistic distribution
    hist_signal = GetMassDistribution(0)
    hist_bkg = GetMassDistribution(2)

    h_teststatistic = ro.TH1F("test_statistic_{}".format(mode),"", 300,-30.,30.)
    for i in range(ntoys):
        if i%1000==0: print("calculating ... ({:}/{:})".format(i,ntoys))
        hist_toy = GenerateToyDataSet(hist_template)
        X = Get_TestStatistic(hist_toy, hist_bkg, hist_signal) #number
        h_teststatistic.Fill(X)

    return h_teststatistic

def Quantiles(hist, type=0):
    # From federico, must fix
    frac_1sigma = ro.Math.gaussian_cdf(-1.,1.,0.)
    frac_2sigma = ro.Math.gaussian_cdf(-2.,1.,0.)
    probs = np.array( [frac_2sigma, frac_1sigma, 0.5, 1-frac_1sigma, 1-frac_2sigma] )

    Xvalues = np.zeros(5)
    Xvalues_out = np.zeros(5)

    # Extract quantiles
    hist.GetQuantiles(5, Xvalues, probs)

    for i in range(0,5):
        Xvalues_out[i]=Xvalues[i]
    if type !=0:
        print("Quantiles for {}:".format(type))
        print("      +2 sigma: t = {:.2f}".format(Xvalues_out[0]))
        print("      +1 sigma: t = {:.2f}".format(Xvalues_out[1]))
        print("      +0 sigma: t = {:.2f}".format(Xvalues_out[2]))
        print("      -1 sigma: t = {:.2f}".format(Xvalues_out[3]))
        print("      -2 sigma: t = {:.2f}\n".format(Xvalues_out[4]))

    return Xvalues_out

def Generate_toys():
# Part 3 and 4
    hist_data = GetMassDistribution(3)
    hist_signal = GetMassDistribution(0)
    hist_bkg = GetMassDistribution(2)
    hist_sb = hist_bkg.Clone("hist_sb")
    hist_sb.Add(hist_signal, 1)

    # Monte carlo simulation for test-statistic distributions
    ntoys = 10000
    print("Generating test-statistic distribution for background-only hypothesis")
    hist_distribution_b = Get_TestStatistic_Distribution(hist_bkg, ntoys, "bkg")

    print("Generating test-statistic distribution for signal+background hypothesis")
    hist_distribution_sb = Get_TestStatistic_Distribution(hist_sb, ntoys, "sb")

    myfile = ro.TFile("output/histograms/test_statistic_distribution.root","RECREATE")
    hist_distribution_b.Write()
    hist_distribution_sb.Write()
    #hist_distribution_data.Write()
    myfile.Close()

def analyze_distributions(h_teststat_bkg, h_teststat_sb, t_data, plot=0):
	# Get the quantiles and print
	quant_b = Quantiles(h_teststat_bkg, "b-only")
	quant_sb = Quantiles(h_teststat_sb, "s+b")

	#histograms for the 1 and 2 sigma indications
	hist_1sigma = h_teststat_bkg.Clone("hist_1sigma"); hist_1sigma.Reset()
	hist_2sigma = h_teststat_bkg.Clone("hist_2sigma"); hist_2sigma.Reset()

	bins_quantiles = [h_teststat_bkg.FindBin(elm) for elm in quant_b]

	for i in range(h_teststat_bkg.GetNbinsX()):
		content = h_teststat_bkg.GetBinContent(i)
		if bins_quantiles[0] <= i <= bins_quantiles[4]:
			hist_2sigma.SetBinContent(i, content)
		if bins_quantiles[1] <= i <= bins_quantiles[3]:
			hist_1sigma.SetBinContent(i, content)

	# Part 5: Find 1-CLb for the distributions
	CLb_b = stat.IntegrateFromRight(h_teststat_bkg, quant_b[2])
	CLb_sb = stat.IntegrateFromRight(h_teststat_bkg, quant_sb[2])
	CLb_data = stat.IntegrateFromRight(h_teststat_bkg, t_data)

	zb_b = ro.Math.gaussian_quantile_c(1.-CLb_b, 1)
	zb_sb = ro.Math.gaussian_quantile_c(1.-CLb_sb, 1)
	zb_data = ro.Math.gaussian_quantile_c(1.-CLb_data, 1)

	# Part 6: Find CLs+b for the distributions
	CLsb_b = stat.IntegrateFromRight(h_teststat_sb, quant_b[2])
	CLsb_sb = stat.IntegrateFromRight(h_teststat_sb, quant_sb[2])
	CLsb_data = stat.IntegrateFromRight(h_teststat_sb, t_data)

	median_b = quant_b[2]
	median_sb = quant_sb[2]
	# Print a summary:
	print("        Median     Median     Median ")
	print("        b-only      s+b        data")
	print("t:     {:.2e}  {:.2e}  {:.2e}".format(median_b, median_sb, t_data))
	print("1-CLb: {:.2e}   {:.2e}   {:.2e}".format(1-CLb_b, 1-CLb_sb, 1-CLb_data))
	print("      ({:.6f}) ({:.6f}) ({:.6f})".format(zb_b, zb_sb, zb_data))
	print("CLs+b: {:.2e}   {:.2e}   {:.2e}".format(CLsb_b, CLsb_sb, CLsb_data))
	print("CLs:   {:.2e}   {:.2e}   {:.2e}".format(CLsb_b/CLb_b, CLsb_sb/CLb_sb, CLsb_data/CLb_data))

	if plot == "plot":
		print("\nPlotting the test-statistic distributions:")
		c = ro.TCanvas("c", "c", 1000, 600)

		h_teststat_sb.SetLineColor(ro.kRed)
		h_teststat_sb.SetTitle(";t=-2ln(Q); number of pseudo-experiments")
		h_teststat_bkg.SetLineColor(ro.kBlack)
		h_teststat_bkg.SetTitle(";t=-2ln(Q); number of pseudo-experiments")
		hist_1sigma.SetFillColor(ro.kGreen)
		hist_2sigma.SetFillColor(ro.kYellow)

		arrow = ro.TArrow(t_data, 2, t_data, 200, 0.01,"<|")
		legend = ro.TLegend(0.75, 0.79, 0.93, 0.93)
		legend.AddEntry(h_teststat_bkg, "b-only", "l")
		legend.AddEntry(h_teststat_sb, "s+b", "l")
		legend.AddEntry(hist_1sigma, "1#sigma", "f")
		legend.AddEntry(hist_2sigma, "2#sigma", "f")

		h_teststat_bkg.Draw("hist")
		hist_2sigma.Draw("same hist")
		hist_1sigma.Draw("same hist")
		h_teststat_sb.Draw("same hist")
		arrow.Draw()
		legend.Draw()
		text(0.39, 0.83, "data", 0.033)
		c.Update()
		Quiet(c.SaveAs)("output/figures/distribution_test_statistic.png")
		print("      figure saved in output/figures/distribution_test_statistic.png\n")
		c.Draw()
	#return c

if __name__ == "__main__":
	#s = count_events(0)
	#b = 5.134
	#db = 0.305
	#s = 5.96
	#b = 7.10
	#db = 0.41
	#Ntoys = 1e6
	#p = ExpectedSignificance_ToyMC(b, s, db ,Ntoys)
	#hist_data = GetMassDistribution(3)
	hist_signal = GetMassDistribution(0)
	MassPlot(20)
	#hist_bkg = GetMassDistribution(2)

	#xx = Get_TestStatistic(hist_data, hist_bkg, hist_signal)

	#hist_test = GenerateToyDataSet(hist_signal)
	# print("TestStatistics = ",higgs.GetTestStatistics(h_data, h_bgd, h_sig))
