import numpy as np
import pyUngewiss as pu

# from tikzplotlib import tikz_save


def HockettSherby(p, x):
    sigma = p[0] + p[1] - p[1] * np.exp(-p[2] * x ** p[3])
    return sigma


def FunHockettSherby(pUnc):
    Prob = pu.UncertainAnalysis(HockettSherby, pUnc)
    Prob.deltax = 1e-3
    Prob.epsStop = 1e-3
    nS = 251
    epsilonMax = 0.5
    epsilon = np.linspace(0, epsilonMax, nS)
    rFnUnc = [[]] * nS
    nEvaluation = 0
    for i, val in enumerate(epsilon):
        Prob.para = val
        Prob.calculate()
        rFnUnc[i] = Prob.rUnc
        nEvaluation += Prob.nEval
    Prob.rFnUnc = rFnUnc
    Prob.nEval = nEvaluation
    Prob.epsilon = epsilon
    return Prob


# Interval parameters
sigmaYInt = pu.UncertainNumber([240, 260])
sigmaPInt = pu.UncertainNumber([40, 60])
cHSInt = pu.UncertainNumber([8, 12])
nHSInt = pu.UncertainNumber([0.7, 0.8])
pInt = [sigmaYInt, sigmaPInt, cHSInt, nHSInt]

# Fuzzy parameters
nAlpha = 6
sigmaYFuzz = pu.UncertainNumber(
    [230, 240, 260, 270], Form='trapazoid', nalpha=nAlpha
)
sigmaPFuzz = pu.UncertainNumber(
    [35, 40, 60, 65], Form='trapazoid', nalpha=nAlpha
)
cHSFuzz = pu.UncertainNumber([7, 8, 12, 13], Form='trapazoid', nalpha=nAlpha)
nHSFuzz = pu.UncertainNumber(
    [0.6, 0.7, 0.8, 0.9], Form='trapazoid', nalpha=nAlpha
)
pFuzz = [sigmaYFuzz, sigmaPFuzz, cHSFuzz, nHSFuzz]

# Interval analysis
ProbInt = FunHockettSherby(pInt)
pu.plotUncertainFn(
    ProbInt.rFnUnc,
    ProbInt.epsilon,
    ylimits=[0, 340],
    xlimits=[0, np.max(ProbInt.epsilon)],
    xlabel='plastic strain $\\varepsilon_{pl}$ [-]',
    ylabel='uncertain stress $\\tilde{\\sigma}$ [MPa]',
    xAxisRot=False,
    color='tab:blue',
    pdpi=1000,
    fill=True,
    xsize=3.5,
    ysize=3,
)
