"""
Test function for interval analysis, in which the classic problem of interval
depedency is shown with IntPy and PyInteval, if avaialble, and how FuzzAnPy
avoids this problem.
"""
import pyUngewiss as pu
from termcolor import colored

try:
    import intpy

    IntPy = True
except:
    IntPy = False
    print("Package 'IntPy' not installed!")
    rInt1 = []
try:
    import interval

    PyInterval = True
except:
    PyInterval = False
    print("Pacakge 'PyInterval' not installed!")
    rInt2 = []
try:
    import uncertainties

    PyUncertainties = True
except:
    PyUncertainties = False
    print("Package 'uncertainties' not installed!")


def SysEq(p, x):
    return p - p


nAlpha = 1
pL = 1.0
pU = 5.0
pInt = pu.UncertainNumber([1, 5])
Prob = pu.UncertainAnalysis()
Prob.SysEq = SysEq
Prob.pUnc = pInt
Prob.calculate()

if IntPy:
    pInt1 = intpy.IReal(pL, pU)
    rInt1 = pInt1 - pInt1
if PyInterval:
    pInt2 = interval.interval(pL, pU)
    rInt2 = pInt2 - pInt2
if PyUncertainties:
    pMean = (pU - pL) / 2.0
    pDelta = pMean - pL
    pUnc = uncertainties.ufloat(pMean, pDelta)
    rUnc = pUnc - pUnc

print('-' * 50)
print(
    colored(
        'Interval analysis via alpha-level optimization with package '
        + "'pyUngewiss' optimization gives: "
        + str(Prob.rUnc.Value),
        'blue',
    )
)
print(
    colored('Number of evaluations: ' + str(Prob.OutputData['nEval']), 'blue')
)
print()
if IntPy:
    print(
        colored(
            "Standard interval analyis with package 'IntPy' gives: "
            + str(rInt1),
            'red',
        )
    )
    print()
if PyInterval:
    print(
        colored(
            "Standard interval analyis with package 'PyInterval' gives: "
            + str(rInt2),
            'red',
        )
    )
    print(colored('Uncertainty OVERESTIMATED!!!', 'red'))
    print()
if PyUncertainties:
    print(
        colored(
            "Uncertainty analyis with package 'uncertainties' gives: "
            + str(rUnc),
            'green',
        )
    )
print()
print('Check complete!')
