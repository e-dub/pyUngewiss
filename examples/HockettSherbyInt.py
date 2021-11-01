import numpy as np
import pyUngewiss as pu


def HockettSherby(p, x):
    sigma = p[0] + p[1] - p[1] * np.exp(-p[2] * x ** p[3])
    return [sigma]


sigmaY = pu.UncertainNumber([240, 260])
sigmaP = pu.UncertainNumber([40, 60])
cHS = pu.UncertainNumber([8, 12])
nHS = pu.UncertainNumber([0.7, 0.8])
pUnc = [sigmaY, sigmaP, cHS, nHS]
Prob = pu.UncertainAnalysis(HockettSherby, pUnc)
Prob.deltax = 1e-3
Prob.epsStop = 1e-3
nS = 100
epsilonMax = 0.5
epsilon = np.linspace(0, epsilonMax, nS)
rFnInt = [[]] * nS
nEvaluation = 0

for i, val in enumerate(epsilon):
    Prob.para = val
    Prob.calculate()
    rFnInt[i] = Prob.rUnc
    nEvaluation += Prob.nEval

pu.plotUncertainFn(
    rFnInt,
    epsilon,
    ylimits=[0, 340],
    xlimits=[0, epsilonMax],
    xlabel='plastic strain $\\varepsilon_{pl}$ [-]',
    ylabel='uncertain stress $\\tilde{\\sigma}$',
)
