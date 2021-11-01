import numpy as np
import pyUngewiss as pu

# Inputs
x1 = 1.8783
x2 = 1.7639
ux1 = 0.001
ux2 = 0.001


# Function of estimated damping ratio based on two measures
def DampEst(p, x):
    x1 = p[0]
    x2 = p[1]
    Lambda = np.log(x1 / x2)
    zeta = np.sqrt(Lambda ** 2 / (4 * np.pi ** 2 + Lambda ** 2))
    return zeta


# Interval uncertainty problem: Worst-case interval defined as +/-6sigma
x1Unc = pu.UncertainNumber([x1 - 6 * ux1, x1 + 6 * ux1])
x2Unc = pu.UncertainNumber([x2 - 6 * ux2, x2 + 6 * ux2])
Prob = pu.UncertainAnalysis(DampEst, pUnc=[x1Unc, x2Unc])
Prob.nr = 1
Prob.deltax = 1e-6
Prob.epsStop = 1e-6
Prob.Alg = 'NLPQLP'

# Calculate uncertain response
Prob.calculate()

# Calculation with standard uncertainty using first-order Taylor series
zeta = DampEst([x1, x2], [])
dzetadx1 = 1 / (2 * np.pi * x1)
dzetadx2 = 1 / (2 * np.pi * x2)
uzeta = (
    np.sqrt(ux1 ** 2 * dzetadx1 ** 2 + ux2 ** 2 * dzetadx2 ** 2) / 2 / np.pi
)
zeta1sigma = [zeta - uzeta, zeta + uzeta]
zeta2sigma = [zeta - 2 * uzeta, zeta + 2 * uzeta]
zeta3sigma = [zeta - 3 * uzeta, zeta + 3 * uzeta]
zeta4sigma = [zeta - 4 * uzeta, zeta + 4 * uzeta]
zeta5sigma = [zeta - 5 * uzeta, zeta + 5 * uzeta]
zeta6sigma = [zeta - 6 * uzeta, zeta + 6 * uzeta]

# Output
print('value of zeta = ' + str(zeta))
print('interval of estimated damping ratio, worst case')
Prob.rUnc.printValue()
print()
print('standard uncertainty in zeta = ' + str(uzeta))
print('convidence interval of 1sigma ' + str(zeta1sigma))
print('convidence interval of 2sigma ' + str(zeta2sigma))
print('convidence interval of 3sigma ' + str(zeta3sigma))
print('convidence interval of 4sigma ' + str(zeta4sigma))
print('convidence interval of 5sigma ' + str(zeta5sigma))
print('convidence interval of 6sigma ' + str(zeta6sigma))
