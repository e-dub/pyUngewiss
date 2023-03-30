import numpy as np
import pyUngewiss as pu
from scipy.optimize import rosen, rosen_der


def SysEq(p, x):
    return rosen(p)


def SensEq(p, r, g, x):
    return rosen_der(p)


nLevel = 6
pFuzz = [[]] * 11
pFuzz[0] = pu.UncertainNumber([1, 2, 3, 4], Form='trapazoid', nLevel=nLevel)
pFuzz[1] = pu.UncertainNumber([1, 2, 3, 4], Form='trapazoid', nLevel=nLevel)
pFuzz[2] = pu.UncertainNumber([1, 2, 3, 4], Form='trapazoid', nLevel=nLevel)
pFuzz[3] = pu.UncertainNumber([1, 2, 3, 4], Form='trapazoid', nLevel=nLevel)
pFuzz[4] = pu.UncertainNumber([1, 2, 3, 4], Form='trapazoid', nLevel=nLevel)
pFuzz[5] = pu.UncertainNumber([1, 2, 3, 4], Form='trapazoid', nLevel=nLevel)
pFuzz[6] = pu.UncertainNumber([1, 2, 3, 4], Form='trapazoid', nLevel=nLevel)
pFuzz[7] = pu.UncertainNumber([1, 2, 3, 4], Form='trapazoid', nLevel=nLevel)
pFuzz[8] = pu.UncertainNumber([1, 2, 3, 4], Form='trapazoid', nLevel=nLevel)
pFuzz[9] = pu.UncertainNumber([1, 2, 3, 4], Form='trapazoid', nLevel=nLevel)
pFuzz[10] = pu.UncertainNumber([1, 2, 3, 4], Form='trapazoid', nLevel=nLevel)
x = 1.0
for i, p in enumerate(pFuzz):
    p.plotValue(xlabel='$p_{' + str(i) + '}$')
UncertainProbFD = pu.UncertainAnalysis(SysEq, pUnc=pFuzz)
UncertainProbFD.Alg = 'NLPQLP'
UncertainProbFD.nLevel = nLevel
UncertainProbFD.deltax = (1e-6,)
UncertainProbFD.paraNorm = 1
UncertainProbFD.SBFA = False
UncertainProbFD.Surr = (False,)
UncertainProbFD.epsStop = 1e-6
UncertainProbFD.para = x
UncertainProbFD.SensCalc = 'FD'
UncertainProbFD.calculate()
UncertainProbFD.rUnc.plotValue(xlabel='$r_{FD}$', color='r')

UncertainProbAS = pu.UncertainAnalysis(SysEq, pUnc=pFuzz, SensEq=SensEq)
UncertainProbAS.Alg = 'NLPQLP'
UncertainProbAS.nLevel = nLevel
UncertainProbAS.paraNorm = 0
UncertainProbAS.SBFA = False
UncertainProbAS.Surr = (False,)
UncertainProbAS.epsStop = 1e-6
UncertainProbAS.para = x
UncertainProbAS.calculate()
UncertainProbAS.rUnc.plotValue(xlabel='$r_{AS}$', color='r')
print('-' * 70)
print(
    'Number of evaluations with finite differencing: '
    + str(UncertainProbFD.OutputData['nEval'])
)
print(UncertainProbFD.rUnc.Value)
print()
print(
    'Number of evaluations with analytical sensitivities: '
    + str(UncertainProbAS.OutputData['nEval'])
)
print(UncertainProbAS.rUnc.Value)
print()
print(
    'Reduction in computational effort (number of evaluations) '
    + str(
        round(
            np.double(
                UncertainProbFD.OutputData['nEval']
                - UncertainProbAS.OutputData['nEval']
            )
            / UncertainProbAS.OutputData['nEval']
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
            np.sum(
                np.abs(UncertainProbAS.rUnc.Value - UncertainProbFD.rUnc.Value)
            )
            / np.max(UncertainProbAS.rUnc.Value),
            precision=4,
        )
    )
)
