import numpy as np
import pyUngewiss as pu


def SysEq(p, x):
    y = x ** 2 + p
    return y


pInt = pu.UncertainNumber([-10, 10])
nS = 200
x = np.linspace(-10, 10, nS)
rFnInt = [[]] * nS
Prob = pu.UncertainAnalysis(SysEq, pInt)
Prob.deltax = 1e-6
Prob.epsStop = 1e-6
for ii in range(len(x)):
    Prob.para = x[ii]
    Prob.calculate()
    rFnInt[ii] = Prob.rUnc
plt, _ = pu.plotUncertainFn(rFnInt, x)
plt.show()
