import ROOT as ro
import numpy as np

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
    t.Draw()

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

    c = ro.TCanvas("c", "c", 1000, 600)
    hist_sb.Draw("hist")
    hist_sb.GetYaxis().SetTitle("Events/ {:.1f} GeV".format(hist_sb.GetBinWidth(1)))
    hist_bkg.Draw("same")
    hist_data.Draw("e same")
    legend.Draw("same")
    c.Update()
    Quiet(c.SaveAs)("output/figures/mass_plot_125.png")


if __name__ == "__main__":
    ro.gROOT.SetBatch(ro.kTRUE)
    MassPlot(20)
