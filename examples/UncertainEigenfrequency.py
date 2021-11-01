# -*- coding: utf-8 -*-
import pyUngewiss as pu
import numpy as np


def Eigenfrequency1DoF(p, x):
    m = p[0]
    k = p[1]
    omega0 = np.sqrt(k / m)
    f0 = omega0 / 2 / np.pi
    return f0


m = pu.UncertainNumber([2.0, 2.5])
k = pu.UncertainNumber([40000, 60000])
pUnc = [m, k]

Prob = pu.UncertainAnalysis(Eigenfrequency1DoF, pUnc)
Prob.deltax = 1e-3
Prob.epsStop = 1e-3

Prob.calculate()

m.printValue()
k.printValue()
plt, _ = pu.plotIntervals(
    [m.Value, k.Value], labels=['mass $m$ [kg]', 'stiffness $k$ [N/mm]']
)
plt.show()

f0Unc = Prob.rUnc
f0Unc.printValue()
f0Unc.plotValue(color='r', xlabel='eigenfrequency $f_0$ [Hz]')

Prob.calcRobustness()
print(Prob.SystemRobustness)
print(Prob.SystemRobustnessNorm)
