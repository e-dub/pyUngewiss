'''
Test function for interval analysis, in which the classic problem of interval
depedency is shown with IntPy and PyInteval, if avaialble, and how FuzzAnPy
avoids this problem.
'''
import pyUngewiss as ung
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
    return(p-p)


nAlpha = 1
pL = 1.
pU = 5.
pInt = ung.UncertainNumber([1, 5])
Prob = ung.UncertainAnalysis()
Prob.SysEq = SysEq
Prob.pUnc = pInt
Prob.calculate()

if IntPy:
    pInt1 = intpy.IReal(pL, pU)
    rInt1 = pInt1-pInt1
if PyInterval:
    pInt2 = interval.interval(pL, pU)
    rInt2 = pInt2-pInt2
if PyUncertainties:
    pMean = (pU-pL)/2.
    pDelta = pMean-pL
    pUnc = uncertainties.ufloat(pMean, pDelta)
    rUnc = pUnc-pUnc

print("-"*50)
print(colored("Interval analysis via alpha-level optimization with " +
              "'pyUngewiss' optimization gives: " +
              str(Prob.rUnc.Value), "blue"))
print(colored("Though it needs the following number of evaluations: " +
              str(Prob.OutputData["nEval"]), "blue"))
print()
if IntPy:
    print(colored("Standard interval analyis with 'IntPy' gives: " +
                  str(rInt1), "red"))
    print()
if PyInterval:
    print(colored("Standard interval analyis with 'PyInterval' gives: " +
                  str(rInt2), "red"))
    print(colored("OVERESTIMATED!!!", "red"))
    print()
if PyUncertainties:
    print(colored("Uncertainty analyis with 'uncertainties' gives: " +
                  str(rUnc), "green"))
print()
print("Check complete!")