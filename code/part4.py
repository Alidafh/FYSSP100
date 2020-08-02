import ROOT as ro
import numpy as np
import analysis_functions as analysis
import ctypes
import statistics as stat


# method to analyse the distribution_test_statistic

# Get the test-statistic distributions from file
histograms = ["test_statistic_bkg", "test_statistic_sb", "test_statistic_data"]
infile = ro.TFile("output/histograms/test_statistic_distribution.root", "READ")
hist_distribution_b = infile.Get(histograms[0]).Clone("h_tdist_b")
hist_distribution_sb = infile.Get(histograms[1]).Clone("h_tdist_sb")
infile.Close()


def analyze_distributions(h_teststat_bkg, h_teststat_sb, t_data, plot=0):
	# Get the quantiles and print
	quant_b = analysis.Quantiles(h_teststat_bkg, "b-only")
	quant_sb = analysis.Quantiles(h_teststat_sb, "s+b")

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
		analysis.text(0.39, 0.83, "data", 0.033)
		c.Update()
		#analysis.Quiet(c.SaveAs)("output/figures/distribution_test_statistic.png")
		#print("      figure saved in output/figures/distribution_test_statistic.png\n")
		c.Draw()
	return c

###############################################################################
# New scale factors?
hist_data = analysis.GetMassDistribution(3)
hist_signal = analysis.GetMassDistribution(0)
hist_bkg = analysis.GetMassDistribution(2)

hist_sc_b = hist_bkg.Clone("newbkg"); hist_sc_b.Reset()
hist_sc_s = hist_signal.Clone("newsignal"); hist_sc_s.Reset()

mu = 2.5
alpha = 1.11

for i in range(hist_bkg.GetNbinsX()):
	b = hist_bkg.GetBinContent(i)
	s = hist_signal.GetBinContent(i)
	#hist_sc_b.SetBinContent(i, b*alpha)
	hist_sc_b.SetBinContent(i, b)
	hist_sc_s.SetBinContent(i, s*mu)


hist_sc_sb = hist_bkg.Clone("newsb")
hist_sc_sb.Add(hist_sc_s, 1)


t_data_scaled = analysis.Get_TestStatistic(hist_data, hist_sc_b, hist_sc_s)
print("      t = {:.3f}\n".format(t_data_scaled))

h_sc_sb = analysis.Get_TestStatistic_Distribution(hist_sc_sb, 100, "sb1")
h_sc_b = analysis.Get_TestStatistic_Distribution(hist_sc_b, 100, "bkg1")

c = analyze_distributions(h_sc_b, h_sc_sb, t_data_scaled, "plot")
