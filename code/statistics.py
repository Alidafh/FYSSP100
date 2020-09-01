import numpy as np
import matplotlib.pyplot as plt
import ROOT as ro


def Poisson(lamb, k):
    #return np.log(lamb) * S - np.loggamma(S + 1) - lamb
    return (lamb**k * np.exp(-lamb))/np.math.factorial(k)

def CL_sb(lamb, nobs):
    CL_sb = 0
    for i in range(int(nobs)):
        CL_sb += Poisson(lamb, i)
    return CL_sb

def significance(s, b, n):
    Z_exp = np.nan_to_num(np.sqrt(2*(s+b)*np.log(1+(s/b)) - 2*s))
    Z_obs = np.nan_to_num(np.sqrt(2*n*np.log(n/b) - 2*(n-b)))
    return Z_exp, Z_obs


def sigFromP(pval, chisquare = "yes"):
    if chisquare != "yes": Z = ro.Math.gaussian_quantile_c(expected_p, 1)
    Z = np.sqrt(ro.Math.chisquared_quantile_c(pval*2,1))
    return Z

def pFromSig(sig):
    p = ro.Math.chisquared_cdf_c(pow(sig, 2), 1)/2
    return p

def test_sig(s,b):
    return s/np.sqrt(b)

#########################################################################

def IntegratePoissonFromRight(tobs, lamb):
	integral = 1.
	for i in range(0, int(tobs)):
		integral -= ro.TMath.Poisson(i, lamb)
	return integral


def IntegratePoisson(tobs, lamb):
	integral = 0.
	for i in range(0,int(tobs)):
		integral += ro.TMath.Poisson(i, lamb)
	return integral

def multiplication_error(a, b, da, db):
    z = a*b
    var = (z**2)*((da/a)**2 + (db/b)**2)
    std = np.sqrt(var)
    return z, std

def IntegrateFromRight(h_X_bgr, X_value):
    # calculates p-value
    Nbins = h_X_bgr.GetNbinsX()
    X_bin = h_X_bgr.FindBin(X_value)
    pvalue = h_X_bgr.Integral(X_bin, Nbins)/h_X_bgr.Integral()
    return pvalue

####################################
#   Confidence levels
####################################
def compute_CLb(hist, t_value):
    # Computes CLb from the test statistic for b-only and a test statistic value
    t_bin = hist.FindBin(t_value)
    CLb = hist.Integral(t_bin, hist.GetNbinsX())/hist.Integral()
    return CLb

def compute_CLsb(hist, t_value):
    t_bin = hist.FindBin(t_value)
    CLsb = hist.Integral(t_bin, hist.GetNbinsX())/hist.Integral()
    return CLsb

if __name__ == '__main__':
    s = 425.21142578125
    b = 19732.498291015625
    nobs = 20128.0
    poisson = Poisson(s+b, nobs)
    #CL = CL(s+b, nobs)
