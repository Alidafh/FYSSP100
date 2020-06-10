import ROOT as ro
import numpy as np
from analysis_functions import GetMassDistribution


hist_higgs125 = GetMassDistribution(type=0, scaleFactor=1)

c = ro.TCanvas("c", "c", 1000, 800)
hist_higgs125.Draw()
c.Draw()
