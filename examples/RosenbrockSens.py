import numpy as np
import pyUngewiss as ung
from scipy.optimize import rosen, rosen_der


def SysEq(p, x):
    return(rosen(p))


def SensEq(p, r, g, x):
    return(rosen_der(p))


nAlpha = 6
pFuzz = [[]]*11
pFuzz[0] = ung.UncertainNumber([1, 2, 3, 4], Form="trapazoid", nalpha=nAlpha)
pFuzz[1] = ung.UncertainNumber([1, 2, 3, 4], Form="trapazoid", nalpha=nAlpha)
pFuzz[2] = ung.UncertainNumber([1, 2, 3, 4], Form="trapazoid", nalpha=nAlpha)
pFuzz[3] = ung.UncertainNumber([1, 2, 3, 4], Form="trapazoid", nalpha=nAlpha)
pFuzz[4] = ung.UncertainNumber([1, 2, 3, 4], Form="trapazoid", nalpha=nAlpha)
pFuzz[5] = ung.UncertainNumber([1, 2, 3, 4], Form="trapazoid", nalpha=nAlpha)
pFuzz[6] = ung.UncertainNumber([1, 2, 3, 4], Form="trapazoid", nalpha=nAlpha)
pFuzz[7] = ung.UncertainNumber([1, 2, 3, 4], Form="trapazoid", nalpha=nAlpha)
pFuzz[8] = ung.UncertainNumber([1, 2, 3, 4], Form="trapazoid", nalpha=nAlpha)
pFuzz[9] = ung.UncertainNumber([1, 2, 3, 4], Form="trapazoid", nalpha=nAlpha)
pFuzz[10] = ung.UncertainNumber([1, 2, 3, 4], Form="trapazoid", nalpha=nAlpha)
x = 1.
for i, p in enumerate(pFuzz):
    p.plotValue(xlabel="$p_{"+str(i)+"}$")
UncertainProbFD = ung.UncertainAnalysis(SysEq, pUnc=pFuzz)
UncertainProbFD.Alg = "NLPQLP"
UncertainProbFD.nAlpha = nAlpha
UncertainProbFD.deltax = 1e-6,
UncertainProbFD.paraNorm = 1
UncertainProbFD.SBFA = False
UncertainProbFD.Surr = False,
UncertainProbFD.epsStop = 1e-6
UncertainProbFD.para = x
UncertainProbFD.SensCalc = "FD"
UncertainProbFD.calculate()
UncertainProbFD.rUnc.plotValue(xlabel="$r_{FD}$", color="r")

UncertainProbAS = ung.UncertainAnalysis(SysEq, pUnc=pFuzz, SensEq=SensEq)
UncertainProbAS.Alg = "NLPQLP"
UncertainProbAS.nAlpha = nAlpha
UncertainProbAS.paraNorm = 0
UncertainProbAS.SBFA = False
UncertainProbAS.Surr = False,
UncertainProbAS.epsStop = 1e-6
UncertainProbAS.para = x
UncertainProbAS.calculate()
UncertainProbAS.rUnc.plotValue(xlabel="$r_{AS}$", color="r")
print("-"*70)
print("Number of evaluations with finite differencing: " +
      str(UncertainProbFD.OutputData["nEval"]))
print(UncertainProbFD.rUnc.Value)
print()
print("Number of evaluations with analytical sensitivities: " +
      str(UncertainProbAS.OutputData["nEval"]))
print(UncertainProbAS.rUnc.Value)
print()
print("Reduction in computational effort (number of evaluations) " +
      str(round(np.double(UncertainProbFD.OutputData["nEval"] -
                          UncertainProbAS.OutputData["nEval"]) /
                UncertainProbAS.OutputData["nEval"]*100, 2)) + " %")
print("Difference in solutions (first norm): " +
      str(np.format_float_scientific(np.sum(np.abs(UncertainProbAS.rUnc.Value-UncertainProbFD.rUnc.Value)) /
                                     np.max(UncertainProbAS.rUnc.Value), precision=4)))
