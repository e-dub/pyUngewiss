import pyUngewiss as pu
import numpy as np


def SysEq(p, x):
    return (x + p) ** 2


def SensEq(p, r, g, x):
    return 2 * x + 2 * p


nAlpha = 11
pFuzz = pu.UncertainNumber([8, 9, 11, 12], Form='trapazoid', nalpha=nAlpha)
x = 1.0
ProbFD = pu.UncertainAnalysis(SysEq, pUnc=pFuzz)
ProbFD.Alg = 'NLPQLP'
ProbFD.nAlpha = nAlpha
ProbFD.deltax = (1e-6,)
ProbFD.paraNorm = 0
ProbFD.SBFA = False
ProbFD.Surr = (False,)
ProbFD.epsStop = 1e-6
ProbFD.para = x
ProbFD.SensCalc = 'FD'
ProbFD.calculate()
ProbAS = pu.UncertainAnalysis(SysEq, pUnc=pFuzz, SensEq=SensEq)
ProbAS.Alg = 'NLPQLP'
ProbAS.nAlpha = nAlpha
ProbAS.paraNorm = 0
ProbAS.SBFA = False
ProbAS.Surr = (False,)
ProbAS.epsStop = 1e-6
ProbAS.para = x
ProbAS.calculate()
print('-' * 60)
print(
    'Number of evaluations with finite differencing: '
    + str(ProbFD.OutputData['nEval'])
)
print(ProbFD.rUnc.Value)
print()
print(
    'Number of evaluations with analytical sensitivities: '
    + str(ProbAS.OutputData['nEval'])
)
print(ProbAS.rUnc.Value)
print()
print(
    'Reduction in computational effort (number of evaluations) '
    + str(
        round(
            np.double(ProbFD.OutputData['nEval'] - ProbAS.OutputData['nEval'])
            / ProbAS.OutputData['nEval']
            * 100,
            2,
        )
    )
    + ' %'
)
print(
    'Difference in solutions (first norm): '
    + str(
        np.format_float_scientific(
            np.sum(np.abs(ProbAS.rUnc.Value - ProbFD.rUnc.Value))
            / np.max(ProbAS.rUnc.Value),
            precision=4,
        )
    )
)
