import ROOT as ro
import numpy as np
from statistics import IntegratePoisson, IntegrateFromRight
import statistics as stat
import ctypes
from array import array

def Quiet(func,level = ro.kInfo + 1):
    # Stops ROOT from printing when preforming tasks such as saving images
    def qfunc(*args,**kwargs):
        oldlevel = ro.gErrorIgnoreLevel
        ro.gErrorIgnoreLevel = level
        try:
            return func(*args,**kwargs)
        finally:
            ro.gErrorIgnoreLevel = oldlevel
    return qfunc

def text(x, y, text_string, size=0.022):
    # Draw text in the plots
    t = ro.TLatex()
    t.SetNDC(1)
    t.SetTextAlign(13);
    t.SetTextSize(size);
    t.SetTextFont(42);
    t.DrawLatex(x, y, text_string)

def GetMassDistribution(type, scaleFactor = 1, rebinN=1):
    # Returns the scaled histogram from file "Histograms_fake.root"
    # type: 0 = Higgs125,
    #       1 = Higgs200,
    #       2 = ZZ SM background,
    #       3 = data
    # scaleFactor: 1 if no scaling is neccesary

    ro.TH1.AddDirectory(ro.kFALSE)
    histograms = ["h_m4l_Higgs125_fake", "h_m4l_Higgs200_fake", "h_m4l_ZZ_fake", "h_m4l_data_fake"]
    infile = ro.TFile("input/Histograms_fake.root", "READ")
    h_mass = infile.Get(histograms[type]).Clone("h_mass")
    #h_mass.Sumw2()
    h_mass.Rebin(rebinN)
    h_mass.SetTitle("; m_{4l} [GeV]"+";Events/ {:.1f} [GeV]".format(h_mass.GetBinWidth(1)))
    nbins = h_mass.GetNbinsX()
    for i in range(nbins):
        bin_cont = h_mass.GetBinContent(i)
        h_mass.SetBinContent(i,scaleFactor*bin_cont)
    infile.Close()
    return h_mass

def MassPlot(rebinN):
    # Make a plot of the mass distribution
    # rebiN: rebin histograms for nicer plotting

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
    #
    # Searches for the mass window width that optimises the expected significance
    #   LumiScale - If you want to scale the luminosity,
    #   plot="no" - If a figure is desired, do plot="plot"
    #
    # Returns the maximum expected significance and the corresponding window width
    #   max_exp, bestwidth_exp

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

        legend = ro.TLegend(0.75, 0.79, 0.9, 0.93)
        legend.AddEntry(hist_window_exp, "Expected", "l")

        hist_window_obs.SetLineColor(1)
        hist_window_exp.SetLineColor(4)

        if( abs(LumiScale-1.00) < 0.01 ):
            # Draw the observed luminosity
            hist_window_obs.Draw("L same")
            legend.AddEntry(hist_window_obs, "Observed", "l")

        hist_window_exp.Draw("L same")
        legend.Draw()
        c.Update()
        Quiet(c.SaveAs)("output/figures/expected_mass_window_{}.png".format(LumiScale))
        print("      figure saved in: output/figures/expected_mass_window_{}.png".format(LumiScale))

    # Find the mass window that maximises the expected significance
    max_exp = hist_window_exp.GetMaximum()
    bin_exp = hist_window_exp.GetMaximumBin()
    bestwidth_exp = hist_window_exp.GetXaxis().GetBinCenter(bin_exp)

    max_obs = hist_window_obs.GetMaximum()
    bin_obs = hist_window_obs.GetMaximumBin()
    bestwidth_obs = hist_window_obs.GetXaxis().GetBinCenter(bin_obs)

    if plot != "no" :print("      maximum expected significance: {:.2f} at mass window {:.2f}".format(max_exp, bestwidth_exp))
    if plot != "no" :print("      maximum observed significance: {:.2f} at mass window {:.2f}\n".format(max_obs, bestwidth_obs))

    return max_exp, bestwidth_exp

def find_fit_parameter(hist, name = "", plot=""):
    # Given a likelihood function -2logL calculates the minimum and
    # the 1 sigma uncertainties and returns these.
    # Input
    #   hist: the -2lnL histogram
    #   name: either "bkg" or "sig"
    #   plot: "plot" if you wish to plot the likelihood function
    #         "sideband" if you wish to plot the likelihood function for the sideband fit
    #                    note that this can only be done if rebinN=LumiScale=1 in sideband_fit()
    # Returns
    #   numpy array with the fit-parameter, lower unc. and upper unc.

    hist.GetYaxis().SetTitle("-2 ln L")
    min_bin = hist.GetMinimumBin()
    min_value = hist.GetBinContent(min_bin)
    scalefactor = hist.GetBinCenter(min_bin)

    h_deltaL = hist.Clone("hist_rescaled")
    h_deltaL.Reset()

    for i in range(1, hist.GetNbinsX()+1):
        cont = min_value - hist.GetBinContent(i)
        h_deltaL.SetBinContent(i, cont)

    low_bin = h_deltaL.FindFirstBinAbove(-1)-1
    up_bin = h_deltaL.FindLastBinAbove(-1)+1
    low_value = hist.GetBinCenter(low_bin)
    up_value = hist.GetBinCenter(up_bin)

    low = scalefactor - low_value
    up = up_value - scalefactor

    print("      Optimal {}. scale factor: {:.3f}, with unc.: -{:.3f}, +{:.3f} ".format(name, scalefactor, low, up))

    if plot == "plot":
        c = ro.TCanvas("c_{}".format(name), "c_{}".format(name), 1000, 600)
        c.cd()
        hist.GetYaxis().SetRangeUser(min_value, hist.GetMaximum())
        hist.Draw("L")
        c.Draw()

        if name == "bkg":
            text(0.768, 0.94, "#hat{#alpha} ="+"{:.3f}".format(scalefactor), 0.042)
            text(0.75, 0.89, "#Delta#hat{#alpha}_{-} ="+"{:.3f}".format(low), 0.042)
            text(0.75, 0.84, "#Delta#hat{#alpha}_{+}="+"{:.3f}".format(up), 0.042)
        if name == "sig":
            text(0.768, 0.94, "#hat{#mu} ="+"{:.3f}".format(scalefactor), 0.042)
            text(0.75, 0.89, "#Delta#hat{#mu}_{-} ="+"{:.3f}".format(low), 0.042)
            text(0.75, 0.84, "#Delta#hat{#mu}_{+}="+"{:.3f}".format(up), 0.042)
        c.Update()
        Quiet(c.SaveAs)("output/figures/proj_loglik_{}.png".format(name))
        print("      figure saved in: output/figures/proj_loglik_{}.png".format(name))

    if (plot=="sideband"):
        c = ro.TCanvas("c{}".format(rebinN), "c", 1000, 600)
        c.cd()
        hist.SetAxisRange(1364, 1368, "Y")
        hist.SetAxisRange(1, 1.25, "X")
        hist.Draw("L")

        line = ro.TLine(bestalpha, 1364, bestalpha, min_value+1)
        line.SetLineStyle(7); line.SetLineColor(ro.kGray+2)
        line.Draw()

        line1 = ro.TLine(1, min_value, 1.25, min_value)
        line1.SetLineStyle(7); line1.SetLineColor(ro.kGray+2)
        line1.Draw()

        line2 = ro.TLine(1, min_value+1, 1.25, min_value+1)
        line2.SetLineStyle(7); line2.SetLineColor(ro.kGray+2)
        line2.Draw()

        line3 = ro.TLine(low_value, 1364, low_value, min_value+1)
        line3.SetLineStyle(7); line3.SetLineColor(ro.kGray+2)
        line3.Draw()

        line4 = ro.TLine(up_value, 1364, up_value, min_value+1)
        line4.SetLineStyle(7); line4.SetLineColor(ro.kGray+2)
        line4.Draw()

        arrow1 = ro.TArrow(low_value, min_value+1, bestalpha, min_value+1, 0.01,"<|>")
        arrow2 = ro.TArrow(bestalpha, min_value+1, up_value, min_value+1, 0.01,"<|>")

        arrow1.Draw()
        arrow2.Draw()

        text(0.38, 0.4, "#Delta#hat#alpha_{-}", size=0.042)
        text(0.56, 0.4, "#Delta#hat#alpha_{+}", size=0.042)

        text(0.04, 0.4, "-2lnL_{max} +1", size=0.024)
        text(0.04, 0.23, "-2lnL_{max}", size=0.024)

        text(0.472, 0.1, "#hat#alpha", size=0.033)
        text(0.26, 0.1, "#hat{#alpha} - #Delta#hat#alpha_{-}", size=0.033)
        text(0.65, 0.1, "#hat{#alpha} + #Delta#hat#alpha_{+}", size=0.033)

        c.Update()
        #c.Draw()
        Quiet(c.SaveAs)("output/figures/neg_loglikelihood_{}.png".format(rebinN))
        print("      figure saved in: output/figures/neg_loglikelihood_{}.png".format(rebinN))

    sf_result = np.array([scalefactor, low, up])
    return sf_result

def sideband_fit(rebinN=1, LumiScale=1):
    #
    # Preforme a likelihood fit to the sidebands, scan the likelihood and find
    # the best parameter for the background and its uncertainties. Also makes a
    # plot of the negative log likelihood as a function of nuicence parameter.
    #
    # Input
    #   rebinN: rebin the histogram (rebinN)
    #   LumiScale: scale the luminosity
    # Returns
    #   bestalpha: The fit parameter that minimises -2lnL,
    #   low: the negative error,
    #   up: the positive error

    hist_bkg = GetMassDistribution(2, LumiScale, rebinN)
    hist_data = GetMassDistribution(3, LumiScale, rebinN)

    nbins = hist_data.GetNbinsX()
    scale_bkg = ro.TH1F("scale_bkg", ";bkg. scale, #alpha; ln L", nbins, 0.1, 3.1)
    hist_likelihood = ro.TH1F("like", ";bkg. scale, #alpha; ln L", nbins, 0.1, 3.1)
    hist_loglik = ro.TH1F("loglik", ";bkg. scale, #alpha; -2ln L", nbins, 0.1, 3.1)

    print("Preforming sideband fit:")
    for i in range(1, nbins+1):
    # Loop over scale factors
        alpha = scale_bkg.GetBinCenter(i)
        logL = 0.
        for j in range(1, nbins+1):
            # Loop over bins
            mass = hist_data.GetBinCenter(j)
            n = hist_data.GetBinContent(j)
            if (mass >= 150. and mass <= 400.):
                func = alpha*hist_bkg.GetBinContent(j)
                logL += np.log(ro.TMath.Poisson(n,func))

        hist_likelihood.SetBinContent(i, logL)
        hist_loglik.SetBinContent(i, -2*logL)

    if rebinN == LumiScale == 1:
        sf_result = find_fit_parameter(hist_loglik, "bkg", "sideband")
    else:
        sf_result = find_fit_parameter(hist_loglik, "bkg")

    bestalpha = sf_result[0]
    low = sf_result[1]
    up = sf_result[2]

    return bestalpha, low, up

def count_events(hist, signal_width, mass = 125.):
    lower_mass_bin = hist.FindBin(mass - 0.5*(signal_width))
    upper_mass_bin = hist.FindBin(mass + 0.5*(signal_width))
    nevents = hist.Integral(lower_mass_bin, upper_mass_bin)
    return nevents

def ExpectedSignificance_ToyMC(b, s, db, Ntoys):
    ro.gROOT.Clear()
    ro.gROOT.Delete()

    print("Calculating the expected significance using {:.2e} toys:".format(Ntoys))
    # Histograms for number of events
    h_Nbkg = ro.TH1D("h_bkg", "Number of background events", 500, -0.5, 499.5);
    h_Nsig_plus_bgr = ro.TH1D("h_sig_bgr", "Number of signal + background events", 500, -0.5, 499.5);

    random = ro.TRandom3(0)

    for i in range(int(Ntoys)):
        mean_b = random.Gaus(b, db)
        mean_splusb = random.Gaus(b, db) + s

        if mean_b >= 0 and mean_splusb >=0:
            b_toy = random.Poisson(mean_b)
            splusb_toy = random.Poisson(mean_splusb)

            h_Nbkg.Fill(b_toy)
            h_Nsig_plus_bgr.Fill(splusb_toy)

    nEvents_median = s+b
    pvalue = IntegrateFromRight(h_Nbkg, nEvents_median)
    #sigma = ro.Math.gaussian_quantile_c(pvalue,1)

    #print("      p-value = {:.3f}, sigificance = {:.3f}".format(pvalue, sigma))
    return pvalue

def Get_TestStatistic(h_mass_dataset, h_template_bkg, h_template_sig):
    #   Computes the likelihood ratio test-statistic for a dataset:
    #    (h_mass_dataset)
    #   from expected distributions for background and signal:
    #    (h_template_bkg, h_template_sig)
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

def Get_TestStatistic_Distribution(ntoys, mode, sf_bkg= 1, sf_sig=1, rebinN=1):
    # makes test-statistic distribution
    hist_signal = GetMassDistribution(0, sf_sig, rebinN)
    hist_bkg = GetMassDistribution(2, sf_bkg, rebinN)
    hist_sb = hist_bkg.Clone("hist_sb")
    hist_sb.Add(hist_signal, 1)

    h_teststatistic = ro.TH1F("test_statistic_{}".format(mode),"", 300,-30.,30.)

    if sf_bkg == 5: h_teststatistic = ro.TH1F("test_statistic_{}".format(mode),"", 300,-90.,90.)

    if mode == "bkg":
        hist_template = hist_bkg.Clone("hist_toy")
        print("Generating distribution for b-only hypothesis")
    if mode == "sb":
        hist_template = hist_sb.Clone("hist_toy")
        print("Generating distribution for s+b hypothesis")

    for i in range(ntoys):
        if i%1000==0: print("calculating ... ({:}/{:})".format(i,ntoys))
        #if i%500==0: print("calculating ... ({:}/{:})".format(i,ntoys))
        hist_toy = GenerateToyDataSet(hist_template)
        X = Get_TestStatistic(hist_toy, hist_bkg, hist_signal) #number
        h_teststatistic.Fill(X)

    return h_teststatistic

def Quantiles(hist, type=0):
    # Finds the median an the 1 and 2 sigma quantiles for the t-distribution.
    # Input
    #   hist: t-dist for the hypothesis
    #   type: option to print the values.
    # Returns
    #   numpy array with the quantiles (-2,-1,median,1,2)

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

def Generate_toys(ntoys, sf_bkg= 1, sf_sig=1, rebinN=1):
    # Calls on the Get_TestStatistic_Distribution function to generate the
    # test-statistic distributions and writes the histograms to file so
    # we do not have to do this again every time the program is edited.
    #
    # TODO: Clean up, this is a very messy way to do this

    hist_distribution_b = Get_TestStatistic_Distribution(ntoys, "bkg", sf_bkg, sf_sig, rebinN)
    hist_distribution_sb = Get_TestStatistic_Distribution(ntoys, "sb", sf_bkg, sf_sig, rebinN)

    #Write the historams to file
    myfile = ro.TFile("output/histograms/test_statistic_distribution_sfb{}_sfs{}_toys{}_bin{}.root".format(sf_bkg, sf_sig, ntoys, rebinN),"RECREATE")
    hist_distribution_b.Write()
    hist_distribution_sb.Write()
    myfile.Close()


def analyze_distributions(h_teststat_bkg, h_teststat_sb, t_data, plot=0):
    # Calculates the CLb, CLsb and CLs values using the test-statistic
    # distributions for bkg and singal and the test statistic from data,
    # and prints a summary.
    # Input:
    #   h_teststat_bkg: t-distribution for the null hypothesis
    #   h_teststat_sb:  t-distribution for the alterative hypothesis
    #   t_data: the value of the test statistic for the data
    #   plot="plot" if you want to plot the test-statistic distributions.
    #
    #  TODO: this function could calculate the test-statistic for the data

    quant_b = Quantiles(h_teststat_bkg, "b-only")
    quant_sb = Quantiles(h_teststat_sb, "s+b")

    median_b = quant_b[2]
    median_sb = quant_sb[2]

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

    # Find 1-CLb for the distributions
    CLb_b = stat.compute_CLb(h_teststat_bkg, median_b)
    CLb_sb = stat.compute_CLb(h_teststat_bkg, median_sb)    #
    CLb_data = stat.compute_CLb(h_teststat_bkg, t_data)

    zb_b = ro.Math.gaussian_quantile_c(1.-CLb_b, 1)
    zb_sb = ro.Math.gaussian_quantile_c(1.-CLb_sb, 1)
    zb_data = ro.Math.gaussian_quantile_c(1.-CLb_data, 1)

    # Find CLs+b for the distributions
    CLsb_b = stat.compute_CLsb(h_teststat_sb, median_b)     #
    CLsb_sb = stat.compute_CLsb(h_teststat_sb, median_sb)
    CLsb_data = stat.compute_CLsb(h_teststat_sb, t_data)

    zsb_b = ro.Math.gaussian_quantile_c(CLsb_b, 1)
    zsb_sb = ro.Math.gaussian_quantile_c(CLsb_sb, 1)
    zsb_data = ro.Math.gaussian_quantile_c(CLsb_data, 1)

    CLs_b = CLsb_b/(CLb_b)
    CLs_sb = CLsb_sb/(CLb_sb)
    CLs_data = CLsb_data/(CLb_data)

    # Print a summary:
    print("        Median     Median     Median ")
    print("        b-only      s+b        data")
    print("t:     {:.5f}  {:.5f}  {:.5f}".format(median_b, median_sb, t_data))
    print("1-CLb: {:.5f}   {:.5f}   {:.5f}".format(1-CLb_b, 1-CLb_sb, 1-CLb_data))
    print("      ({:.6f}) ({:.6f}) ({:.6f})".format(zb_b, zb_sb, zb_data))
    print("CLs+b: {:.5f}   {:.5f}   {:.5f}".format(CLsb_b, CLsb_sb, CLsb_data))
    print("      ({:.6f}) ({:.6f}) ({:.6f})".format(zsb_b, zsb_sb, zsb_data))
    print("CLs:   {:.5f}   {:.5f}   {:.5f}".format(CLs_b, CLs_sb, CLs_data))

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
        legend = ro.TLegend(0.75, 0.79, 0.90, 0.93)
        legend.AddEntry(h_teststat_bkg, "b-only", "l")
        legend.AddEntry(h_teststat_sb, "s+b", "l")
        legend.AddEntry(hist_1sigma, "1#sigma", "f")
        legend.AddEntry(hist_2sigma, "2#sigma", "f")

        #h_teststat_bkg.SetAxisRange(-30,30,"X")
        #h_teststat_sb.SetAxisRange(-30,30,"X")
        #hist_2sigma.SetAxisRange(-30,30,"X")
        #hist_1sigma.SetAxisRange(-30,30,"X")

        h_teststat_bkg.Draw("hist")
        hist_2sigma.Draw("same hist")
        hist_1sigma.Draw("same hist")
        h_teststat_sb.Draw("same hist")
        #arrow.Draw()
        legend.Draw()
        #text(0.39, 0.83, "data", 0.033)
        c.Update()
        Quiet(c.SaveAs)("output/figures/distribution_test_statistic_t{}.png".format(int(abs(t_data))))
        print("      figure saved in output/figures/distribution_test_statistic_t{}.png\n".format(int(abs(t_data))))
        c.Draw()
    return zb_sb, zsb_b

def plotMassSideband(rebinN, scale, dscale):
    # Plot the sideband region in the bkg mass spectrum with the scaled
    # bkg. distribution and the datapoints
    # Input
    #   rebinN: rebin factor for prettier plotting
    #   scale: the scale factor
    #   dscale: error in the scalefactor
    #
    # TODO: fix the error stuff in this function

    print("\nPlotting side band region:")
    hist_data = GetMassDistribution(3, 1, rebinN)
    hist_bkg = GetMassDistribution(2, 1, rebinN)

    hist_scaled_bkg = GetMassDistribution(2, scale, rebinN)
    hist_scaled_bkg.SetFillColor(2)
    hist_scaled_bkg.SetLineColor(2)

    hist_bkg.SetFillColor(0)
    hist_bkg.SetLineColor(1)
    #hist_bkg.SetFillStyle(3004)
    """
    # Set the errors on the histograms
    for i in range(hist_scaled_bkg.GetNbinsX()):
        db = hist_bkg.GetBinError(i)
        if b != 0:
            scaled_b, d_scaled_b = stat.multiplication_error(scale, b, dscale, db)
            hist_scaled_bkg.SetBinError(i, d_scaled_b)
    """
    hist_bkg.SetAxisRange(150, 400, "X")
    hist_data.SetAxisRange(150, 400, "X")
    hist_scaled_bkg.SetAxisRange(150, 400, "X")

    hist_error_band = hist_scaled_bkg.Clone("new")
    hist_error_band.SetFillColor(1)
    hist_error_band.SetLineColor(0)
    hist_error_band.SetFillStyle(3004)
    hist_error_band.SetMarkerColor(0)
    hist_error_band.SetMarkerStyle(0)

    legend = ro.TLegend(0.65, 0.75, 0.93, 0.93)
    legend.AddEntry(hist_bkg, "Unscaled Background", "l")
    legend.AddEntry(hist_scaled_bkg, "Scaled Background (scale={:.2f})".format(scale), "f")
    #legend.AddEntry(hist_error_band, "unc. on scaled bkg", "f")
    legend.AddEntry(hist_data, "Data", "p")

    c = ro.TCanvas("c", "c", 1000, 600)
    hist_data.Draw("e")
    hist_scaled_bkg.Draw("hist same")
    #hist_error_band.Draw("E3 same")
    hist_bkg.Draw("hist same")
    hist_data.Draw("e same")
    legend.Draw()
    c.Draw()
    c.Update()
    Quiet(c.SaveAs)("output/figures/sideband_scaled.png")
    print("      figure saved in: output/figures/sideband_scaled.png\n")

def signal_fit(alpha, rebinN=1, LumiScale=1):
    # Do a fit to find the signal scale factor, the bacground parameter is set to 1
    #
    # TODO: Not actually used in the analysis. fix.

    hist_sig = GetMassDistribution(0, LumiScale)
    hist_bkg = GetMassDistribution(2, LumiScale)
    hist_data = GetMassDistribution(3, LumiScale)

    hist_sig.Rebin(rebinN)
    hist_bkg.Rebin(rebinN)
    hist_data.Rebin(rebinN)

    nbins = hist_data.GetNbinsX()

    scale_sig = ro.TH1F("scale_sig", "", nbins, 0.1, 3.1)

    hist_lik_sig = ro.TH1F("Lsig", ";sig. scale, #mu; -2#Delta ln L", nbins, 0.1, 3.1)
    hist_likelihood = ro.TH1F("Lsig2", ";#mu; -2#Delta ln L", nbins, 0.1, 3.1)

    print("Preforming likelihood fit in signal region:")
    for i in range(1, nbins+1):
        # Loop over scale factors
        mu = scale_sig.GetBinCenter(i)

        logL = 0.
        for j in range(1, nbins+1):
            mass = hist_data.GetBinCenter(j)
            n = hist_data.GetBinContent(j)
            b = hist_bkg.GetBinContent(j)
            s = hist_sig.GetBinContent(j)

            func = mu*s + alpha*b
            logL += np.log(ro.TMath.Poisson(n, func))

            #if (mass >= 121 and mass <= 129):
            #    func = mu*s + alpha*b
            #    logL += np.log(ro.TMath.Poisson(n, func))

        hist_lik_sig.SetBinContent(i, -2*logL)
        hist_likelihood.SetBinContent(i, logL)

    # Find the minimum
    min_bin = hist_lik_sig.GetMinimumBin()
    min_val = hist_lik_sig.GetBinContent(min_bin)
    bestmu = hist_lik_sig.GetBinCenter(min_bin)

    #print(bestmu)

    # Find the maximum likelihood value
    maximum_bin = hist_likelihood.GetMaximumBin()
    maximum_value = hist_likelihood.GetBinContent(maximum_bin)
    bestmu2 = hist_likelihood.GetBinCenter(maximum_bin)

    #print(bestmu2)

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

    low = bestmu - low_value
    up = up_value - bestmu

    print("      Optimal sig. scale factor: {:.3f}, with unc.: -{:.3f}, +{:.3f}".format(bestmu, low, up))

    c = ro.TCanvas("s", "s", 1000, 600)
    hist_lik_sig.Draw("l")
    #h_likelihood_rescaled.Draw("l")
    c.Draw()
    Quiet(c.SaveAs)("output/figures/signal_scale_likelihood.png".format(rebinN))
    print("      figure saved in: output/figures/signal_scale_likelihood.png".format(rebinN))

def find_fit_parameter_2d(hist, plot=""):

    min_bin = hist.GetMinimumBin()
    xx, yy, zz = ctypes.c_int(0), ctypes.c_int(0), ctypes.c_int(0)
    hist.GetBinXYZ(min_bin, xx, yy, zz)

    min_val = hist.GetBinContent(xx.value, yy.value)
    bestalpha = hist.GetXaxis().GetBinCenter(xx.value)
    bestmu = hist.GetYaxis().GetBinCenter(yy.value)

    dhist = hist.Clone("new")
    for i in range(1, hist.GetNbinsX()+1):
        for j in range(1, hist.GetNbinsX()+1):
            dhist.SetBinContent(i, j, min_val - hist.GetBinContent(i, j))

    proj_sig = dhist.ProjectionY("dmu")
    proj_bkg = dhist.ProjectionX("dalpha")

    co = ro.TCanvas("c0", "c0", 1000, 600)
    co.Divide(2)
    co.cd(1)
    proj_bkg.Draw("L")
    co.cd(2)
    proj_sig.Draw("L")
    co.Draw()
    #c.Update()
    Quiet(co.SaveAs)("output/figures/proj_loglik.png")
    print("      figure saved in: output/figures/proj_loglik_.png")

    return np.array([bestmu, low1, up1]), np.array([bestalpha, low2, up2])

def muFit(rebinN, nscale, plot = "", sf_bkg=1, sf_sig=1):
    # Prefomes a likelihood fit to both the signal and bacground scale factor.
    # uses the projections of -2lnL to find the best parameters and their 1 sigma
    # uncertainty. Writes the histograms to file so these can be manipulated/displayed
    # without having to do the fit again. (takes some time)
    #
    # Input
    #   rebinN: rebin histogram
    #   nscale: number of scalefactors
    #   sf_bkg: scalefactor for the background
    #   sf_sig: scalefactor for the signal
    #   plot:   plot the countours and the projections

    hist_sig = GetMassDistribution(0, sf_sig, rebinN)
    hist_bkg = GetMassDistribution(2, sf_bkg, rebinN)
    hist_data = GetMassDistribution(3, 1, rebinN)

    nbins = hist_data.GetNbinsX()
    hist_scalefactors = ro.TH2D("scale",";bkg. scale, #alpha; sig.scale, #mu; -2lnL", nscale, 0.1, 3.1, nscale, 0.1, 3.1)

    hist_loglik = hist_scalefactors.Clone("loglik")
    hist_loglik.Reset()

    print("\nPreforming likelihood fit for both parameters:")
    for i in range(1, hist_scalefactors.GetNbinsX()+1):
        for j in range(1, hist_scalefactors.GetNbinsX()+1):
            # Loops over scalefactors
            alpha = hist_scalefactors.GetXaxis().GetBinCenter(i)
            mu = hist_scalefactors.GetYaxis().GetBinCenter(j)

            logL = 0.
            for k in range(1, nbins+1):
                # loop over bins in datasets
                n = hist_data.GetBinContent(k)
                b = hist_bkg.GetBinContent(k)
                s = hist_sig.GetBinContent(k)

                func = mu*s + alpha*b
                logL += np.log(ro.TMath.Poisson(n, func))

            hist_loglik.SetBinContent(i, j, -2*logL)

    min_bin = hist_loglik.GetMinimumBin()
    xx, yy, zz = ctypes.c_int(0), ctypes.c_int(0), ctypes.c_int(0)
    hist_loglik.GetBinXYZ(min_bin, xx, yy, zz)

    min_val = hist_loglik.GetBinContent(xx.value, yy.value)
    bestalpha = hist_loglik.GetXaxis().GetBinCenter(xx.value)
    bestmu = hist_loglik.GetYaxis().GetBinCenter(yy.value)

    print("      Optimal sig. scale factor: {:.3f}".format(bestmu))
    print("      Optimal bkg. scale factor: {:.3f}".format(bestalpha))

    dhist = hist_loglik.Clone("new")

    for i in range(1, hist_loglik.GetNbinsX()+1):
        for j in range(1, hist_loglik.GetNbinsX()+1):
            cont = hist_loglik.GetBinContent(i, j)
            if cont <=1: dhist.SetBinContent(i, j, 1)

    if plot=="plot":
        c_par = ro.TCanvas("c_par", "c_par", 1000, 600)
        c_par.cd()
        hist_loglik.SetAxisRange(0.9, 1.3, "X")
        hist_loglik.SetAxisRange(0.2, 2.8, "Y")
        m = ro.TMarker(bestalpha, bestmu, 20)
        #hist_loglik.Draw("colz")
        dhist.Draw("cont3")
        m.Draw()
        c_par.Update()
        Quiet(c_par.SaveAs)("output/figures/countour_fit_parameters_{}.png".format(rebinN))
        print("      figure saved in: output/figures/countour_fit_parameters_{}.png".format(rebinN))

    print("\nUsing the projections: ")
    # use projections to try to find some uncertainties
    h_proj_bkg = hist_loglik.ProjectionX("hpX")
    h_proj_sig = hist_loglik.ProjectionY("hpY")

    sf_sig = find_fit_parameter(h_proj_sig, "sig", "plot")
    sf_bkg = find_fit_parameter(h_proj_bkg, "bkg", "plot")
    return hist_loglik, sf_sig, sf_bkg

def lumi_scale(ntoys, rebinN = 1):
    # Scan over different luminosities and calculate the CLb value
    hist_lumi = ro.TH1D("lumi", "", 10, 0, 10)

    for i in range(1, hist_lumi.GetNbinsX()+1):
        LumiScale = hist_lumi.GetBinCenter(i)
        sf_bkg = LumiScale
        sf_sig = LumiScale
        print("Calculating for luminosity scale:", LumiScale)

        # Get test-stat for bkg
        hist_t_bkg = Get_TestStatistic_Distribution(ntoys, "bkg", sf_bkg, sf_sig, rebinN)
        median_bkg = Quantiles(hist_t_bkg)[2]

        #test-stat for sb
        hist_t_sb = Get_TestStatistic_Distribution(ntoys, "sb", sf_bkg, sf_sig,rebinN)
        median_sb = Quantiles(hist_t_sb)[2]

        # 1-CLb
        CLb_sb = 1 - stat.compute_CLb(hist_t_bkg, median_sb)
        if (CLb_sb) < 1e-9: CLb_sb = 1e-9

        significance = ro.Math.gaussian_quantile_c(CLb_sb, 1)
        print("\n      median bkg:{:.4f},   median sb: {:.4f}".format(median_bkg, median_sb))
        print("      expected CLb is {:.4f} ({:.4f})\n".format(1-CLb_sb, significance))

        hist_lumi.SetBinContent(i, significance)

    c = ro.TCanvas("c", "c", 1000, 600)
    c.cd()
    hist_lumi.Draw("P*")
    c.Draw()

    myfile = ro.TFile("output/histograms/clb_toys{}_bin{}.root".format(ntoys, rebinN),"RECREATE")
    hist_lumi.Write()
    myfile.Close()
    print ("histogram saved in output/histograms/clb_toys{}_bin{}.root".format(ntoys, rebinN))

def extrapolate(CL_list, type):
    # Takes three points and extrapolates
    # make graph
    print("Extrapolating using 3 points:")
    lumi, z_cl = array('d', [1.0, 1.5, 2.0]), array( 'd', CL_list)
    axis = "Expected {}".format(type)

    gr_lumi_sign = ro.TGraph(3, lumi, z_cl)

    fit = ro.TF1("fit", "pol1", 1, 6)
    gr_lumi_sign.Fit("fit", "Q0")
    func = gr_lumi_sign.GetFunction("fit")
    parameters = [func.GetParameter(i) for i in range(2)]

    hist_fit = ro.TH1F("h_sb", "", 100, 0, 6)
    hist_fit.Eval(func)

    lumi5 = 0
    z_val5 = 0
    for i in range(1, hist_fit.GetNbinsX()+1):
        z_val = hist_fit.GetBinContent(i)
        lumi = hist_fit.GetBinCenter(i)
        if z_val >= 5.0:
            lumi5 = lumi
            z_val5 = z_val
            break

    print("      {:.2f} sigma reached at {:.2f} x luminosity".format(z_val, lumi))

    hist_fit.SetLineStyle(7)
    hist_fit.SetLineColor(ro.kGray)
    hist_fit.GetXaxis().SetTitle("Luminosity scale")
    hist_fit.GetYaxis().SetTitle(axis)
    hist_fit.SetAxisRange(0,6,"Y")
    max = ro.TLine(0, 5, 6, 5)
    max.SetLineStyle(7)
    max.SetLineColor(ro.kRed)

    legend = ro.TLegend(0.40, 0.25, 0.85, 0.5)
    legend.AddEntry(hist_fit, "Extrapolation from points", "l")
    legend.AddEntry(gr_lumi_sign, "expected significance using 10^{4} toys")
    legend.AddEntry(max, "the famous 5#sigma limit")

    c = ro.TCanvas("c{}".format(type), "cc", 1000, 600)
    hist_fit.Draw("L")
    gr_lumi_sign.Draw("P")
    max.Draw()
    legend.Draw()
    c.Draw()
    c.Update()

    Quiet(c.SaveAs)("output/figures/extrapolation_{}.png".format(type))
    print("      figure saved in: output/figures/extrapolation_{}.png\n".format(type))
    #return c



if __name__ == "__main__":
    #lumi_scale(100, 10)
    #c11 = extrapolate("clb")
    #c22 = extrapolate("clsb")
    #sideband_fit(rebinN=1, LumiScale=1)
    #CL_list = [2.06, 2.45, 2.85]
    #extrapolate(CL_list, "1-CLb")
    muFit(10, 10)

    #f = ro.TFile("output/histograms/loglik_2d_rebinN10_nscale100.root", "READ")
    #hist_loglik = f.Get("loglik").Clone("h_loglik")
    #hist_loglik.SetDirectory(0)
    #f.Close()

    #find_fit_parameter_2d(hist_loglik, "plot")
