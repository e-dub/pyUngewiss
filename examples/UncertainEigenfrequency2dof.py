import pyUngewiss as pu
import numpy as np
import scipy.linalg as linalg


def Eigenfrequence2DoF(p, x):
    m1 = p[0]
    m2 = p[1]
    k1 = p[2]
    k2 = p[3]
    M = np.array([[m1,  0],
                  [ 0, m2]])
    K =  np.array([[k1+k2, -k2],
                   [  -k2,  k2]])
    L, Phi = linalg.eigh(K, M)
    omega = np.sqrt(L)
    f0 = omega/2/np.pi
    return(f0)

def Jacobian():
    #...
    #coming soon

# Uncertain parameters
m1 = 1
m2 = 0.5
k1 = 1e7
k2 = 2e7
m1 = pu.UncertainNumber([0.75, 1.25])   # 1
m2 = pu.UncertainNumber([0.25, 0.75])  # 0.5
k1 = pu.UncertainNumber([0.75e7, 1.25e7])  # 0.5
k2 = pu.UncertainNumber([1.75e7, 2.25e7])  # 0.5
pUnc = [m1, m2, k1, k2]

# Set up uncertain problem
Prob = pu.UncertainAnalysis(Eigenfrequence2DoF, pUnc)
Prob.nr = 2
Prob.deltax = 1e-3
Prob.epsStop = 1e-3
Prob.Alg = "NLPQLP"

# Calculate uncertain response
Prob.calculate()
print(Prob.nEval)

# Uncertain parameters
m1.printValue()
m2.printValue()
k1.printValue()
k2.printValue()
pu.plotIntervals([m1.Value, k1.Value],
                 labels=["mass $m$ [kg]", "stiffness $k$ [N/mm]"])

# Uncertain reponse
f0Unc = Prob.rUnc
f0Unc[0].printValue()
f0Unc[1].printValue()
pu.plotIntervals(f0Unc, color="r", labels=["eigenfrequency $f_{0,1}$",
                                           "eigenfrequency $f_{0,2}$"],
                 xlabel="frequency [Hz]")

# Robustness
Prob.calcRobustness()
print(Prob.SystemRobustness)
print(Prob.SystemRobustnessNorm)

# Shadow uncertainty
print(Prob.ShadowUncertainty)