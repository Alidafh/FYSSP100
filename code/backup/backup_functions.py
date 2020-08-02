def sideband_fit(rebinN):
    #Sideband fit
    hist_bkg = GetMassDistribution(2)
    hist_data = GetMassDistribution(3)

    hist_bkg.Rebin(rebinN)
    hist_data.Rebin(rebinN)

    nbins_data = hist_data.GetNbinsX()
    scale_bkg = ro.TH1F("scale_bkg", ";#alpha; -2lnL", 100, 0.1, 3.1)

    print("Preforming sideband fit:")
    for i in range(1, scale_bkg.GetNbinsX()+1):
        alpha = scale_bkg.GetBinCenter(i)
        logL = 0.
        #logL = 1e-10
        for j in range(1, nbins_data+1):
            mass = hist_data.GetBinCenter(j)
            n = hist_data.GetBinContent(j)
            if (mass >= 150. and mass <= 400.):
                func = alpha*hist_bkg.GetBinContent(j)
                #logL += np.log(ro.TMath.Poisson(n, func))
                if func > 0.: logL += np.log(ro.TMath.Poisson(n, func))

        scale_bkg.SetBinContent(i, -2.*logL)

    # Find minimum
    min_bin = scale_bkg.GetMinimumBin()
    min_likelihood = scale_bkg.GetBinContent(min_bin)
    bestalpha = scale_bkg.GetBinCenter(min_bin)
    """
    for i in range(scale_bkg.GetNbinsX()):
        cont = scale_bkg.GetBinContent(i)
        if cont <= min_likelihood+1:
            print("yes")
            break
    """
    hist_limits= ro.TH1F("scale_bkg_rescaled", ";#alpha; -2lnL", 100, 0.1, 3.1)
    for i in range(1,scale_bkg.GetNbinsX()+1):
        hist_limits.SetBinContent(i, min_likelihood-scale_bkg.GetBinContent(i))

    left = scale_bkg.GetBinCenter(hist_limits.FindFirstBinAbove(-1)-1)
    right = scale_bkg.GetBinCenter(hist_limits.FindLastBinAbove(-1)+1)

    sigma1 = bestalpha - left
    sigma2 = right - bestalpha

    print("      Optimal bkg. scale factor: {:.4f}, with unc.: -{:.4f}, +{:.4f} ".format(bestalpha, sigma1, sigma2))

    c = ro.TCanvas("c0", "", 1000, 600)
    scale_bkg.Draw("L")
    min = 0.9
    max = 1.3
    scale_bkg.GetYaxis().SetRangeUser(421, 424)
    scale_bkg.SetAxisRange(min, max, "X")

    # Horizontal line on minimum
    linemin = ro.TLine(min , min_likelihood, max, min_likelihood)
    linemin.SetLineStyle(7); linemin.SetLineColor(ro.kGray+2)
    #text(0.099, 0.919, "ln L(#hat#alpha)", size=0.032)
    linemin.Draw()

    linemax = ro.TLine(min, min_likelihood+1, max, min_likelihood+1)
    linemax.SetLineStyle(7); linemax.SetLineColor(ro.kGray+2)
    #text(0.099, 0.919, "ln L(#hat#alpha)", size=0.032)
    linemax.Draw()

    c.Update()
    c.Draw()
    return c
