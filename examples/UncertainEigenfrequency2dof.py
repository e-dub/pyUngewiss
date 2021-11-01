import pyUngewiss as pu
import numpy as np
import scipy.linalg as linalg


def Eigenfrequence2DoF(p, x):
    m1 = p[0]
    m2 = p[1]
    k1 = p[2]
    k2 = p[3]
    M = np.array([[m1, 0], [0, m2]])
    K = np.array([[k1 + k2, -k2], [-k2, k2]])
    L, Phi = linalg.eigh(K, M)
    omega = np.sqrt(L)
    f0 = omega / 2 / np.pi
    return f0


def Jacobian(p, f, g, x):
    m1 = p[0]
    m2 = p[1]
    k1 = p[2]
    k2 = p[3]
    M = np.array([[m1, 0], [0, m2]])
    K = np.array([[k1 + k2, -k2], [-k2, k2]])
    L, Phi = linalg.eigh(K, M)
    omega = np.sqrt(L)
    kNabla = np.array(
        [
            [[0, 0], [0, 0]],
            [[0, 0], [0, 0]],
            [[1, 0], [0, 0]],
            [[1, -1], [-1, 1]],
        ]
    )
    mNabla = np.array(
        [
            [[1, 0], [0, 0]],
            [[0, 0], [0, 1]],
            [[0, 0], [0, 0]],
            [[0, 0], [0, 0]],
        ]
    )
    lambdaNabla = [[]] * len(L)
    omegaNabla = [[]] * len(L)
    f0Nabla = [[]] * len(L)
    eq = 'nHermitian'
    for i in range(len(L)):
        if eq == 'Hermitian':
            lambdaNabla[i] = Phi[:, i].T @ (kNabla - L[i] * mNabla) @ Phi[:, i]
            omegaNabla[i] = lambdaNabla[i] / 2 / omega[i]
            f0Nabla[i] = omegaNabla[i] / 2 / np.pi
        else:  # non-Hermitian
            omegaNabla[i] = (
                -L[i] * Phi[:, i].T @ mNabla @ Phi[:, i]
                + Phi[:, i].T @ kNabla @ Phi[:, i]
            ) / (2 * omega[i] * Phi[:, i].T @ M @ Phi[:, i])
            f0Nabla[i] = omegaNabla[i] / 2 / np.pi
    return f0Nabla


# Uncertain parameters
m1 = 1
m2 = 0.5
k1 = 1e7
k2 = 2e7
p = [m1, m2, k1, k2]

m1 = pu.UncertainNumber([0.75, 1.25])   # 1
m2 = pu.UncertainNumber([0.25, 0.75])  # 0.5
k1 = pu.UncertainNumber([0.75e7, 1.25e7])  # 0.5
k2 = pu.UncertainNumber([1.75e7, 2.25e7])  # 0.5
pUnc = [m1, m2, k1, k2]

# Set up uncertain problem
Prob = pu.UncertainAnalysis(Eigenfrequence2DoF, pUnc, SensEq=Jacobian)
Prob.nr = 2
Prob.deltax = 1e-6
Prob.epsStop = 1e-6
Prob.Alg = 'NLPQLP'

# Calculate uncertain response
Prob.calculate()
print(Prob.nEval)

# Uncertain parameters
m1.printValue()
m2.printValue()
k1.printValue()
k2.printValue()
pu.plotIntervals(
    [m1.Value, k1.Value], labels=['mass $m$ [kg]', 'stiffness $k$ [N/mm]']
)

# Uncertain reponse
f0Unc = Prob.rUnc
f0Unc[0].printValue()
f0Unc[1].printValue()
print('In Hertz')
print((f0Unc[0].Value ** 0.5) / 2 / np.pi)
print((f0Unc[1].Value ** 0.5) / 2 / np.pi)
plt, _ = pu.plotIntervals(
    f0Unc,
    color='r',
    xlabel='frequency [Hz]',
    labels=['eigenfrequency $f_{0,1}$', 'eigenfrequency $f_{0,2}$'],
)
plt.show()

# Robustness
Prob.calcRobustness()
print(Prob.SystemRobustness)
print(Prob.SystemRobustnessNorm)

# Shadow uncertainty
print(Prob.ShadowUncertainty)
Prob.ShadowUncertainty[0, 3, 0, :]
