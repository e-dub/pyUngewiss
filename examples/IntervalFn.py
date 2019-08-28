import numpy as np
import pyUngewiss as ung


def SysEq(p, x):
    y = x**2+p
    return(y)


pInt = ung.UncertainNumber([-10, 10])
nS = 200
x = np.linspace(-10, 10, nS)
rFnInt = [[]]*nS
Prob = ung.UncertainAnalysis(SysEq, pInt)
Prob.deltax = 1e-6
Prob.epsStop = 1e-6
for ii in range(len(x)):
    Prob.para = x[ii]
    Prob.calculate()
    rFnInt[ii] = Prob.rUnc
ung.plotUncertainFn(rFnInt, x)
