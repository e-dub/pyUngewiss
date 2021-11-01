"""
pyUngewiss - Python library for UNcertainty analysis in liGhtwEight dsiGn with
IntervalS and fuzzy numberS

Previously called FuzzAnPy, but now renamed to current name to represent
concentration beyond fuzzy numbers, making intervals default.

This is based on the following work:
Wehrle (2015) Design optimization of lightweight space-frame structures
considering crashworthiness and parameter uncertainty


TODO Save convergence and further information
TODO Pretty plots
TODO Hook up with complicated example

TODO algoptions as object
"""
from pyUngewiss.UncertainNumber import UncertainNumber
from pyUngewiss.OptAlgOptions import AlgOptions

try:
    import pygmo as pg
except:
    pass
try:
    import pyOpt
except:
    pass
try:
    import cma
except:
    pass
import numpy as np
from numpy.linalg import lstsq


np.seterr(divide='ignore', invalid='ignore')

__title__ = (
    'Python library for UNcertainty analysis in liGhtwEight desiGn '
    + 'with IntervalS and fuzzy numberS'
)
__shorttitle__ = 'pyUngewiss'
__version__ = '1.0 - Initial public release'
__all__ = 'pyUngewiss'
__author__ = 'E. J. Wehrle'
__copyright__ = 'Copyright 2019, 2020, 2021: E. J. Wehrle'
__email__ = 'Erich.Wehrle(a)unibz.it'
__license__ = 'GNU General Public License v3.0'
__url__ = 'github.org/e-dub/pyUngewiss'


def printSplash():
    print('')
    print(__shorttitle__ + ' - ' + __title__)
    print('')
    print('Version:     ' + __version__)
    print('Internet:    ' + __url__)
    print('License:     ' + __license__)
    print('Copyright:   ' + __copyright__)
    print('')


class UncertainAnalysis(object):
    def __init__(self, SysEq=[], pUnc=[], SensEq=[]):
        self.SysEq = SysEq
        self.pUnc = pUnc
        self.SensEq = SensEq
        self.nr = 1
        self.Alg = 'NLPQLP'
        if not pUnc:
            self.nAlpha = 1
        else:
            if isinstance(pUnc, list):
                self.nAlpha = np.shape(pUnc[0].Value)[0]
            else:
                self.nAlpha = np.shape(pUnc.Value)[0]
        self.deltax = 1e-2
        self.paraNorm = True
        self.para = []
        self.SBFA = False
        self.Surr = 'Kriging'
        self.SensCalc = 'FD'
        self.epsStop = 1.0e-4
        self.UncertainModel = 'UncertainSystem'
        self.PrintOut = True
        self.nEval = 0

    def calculate(self):
        if hasattr(self.SensEq, '__call__'):
            self.SensCalc = 'OptSensEq'

        def MinSysEq(x):
            x = np.array(x)
            if not self.SBFA:
                r = self.SysEq(x, self.para)
            else:
                r = Surrogate(x, ir)
            if self.Alg == 'MMA':
                g = np.array([0.0])
            else:
                g = []
            self.nEval += 1
            if np.size(r) > 1 or type(r) is list:
                f = r[ir]
            else:
                f = r
            fail = 0
            return f, g, fail

        def MaxSysEq(x):
            f, g, fail = MinSysEq(x)
            f = -f
            return f, g, fail

        def MinSysEqNorm(xNorm):
            xNorm = np.array(xNorm)
            x = denormalize(xNorm, xL, xU)
            f, g, fail = MinSysEq(x)
            return f, g, fail

        def MaxSysEqNorm(xNorm):
            f, g, fail = MinSysEqNorm(xNorm)
            f = -f
            return f, g, fail

        def MinSensEq(x, f, g):
            drdx = self.SensEq(x, f, g, self.para)
            if self.Alg == 'MMA':
                dgdx = np.zeros([1, np.size(x)])
            else:
                dgdx = []
            if self.nr > 1:
                # dfdx = np.array(drdx)[:, ir].reshape([1, len(x)])
                dfdx = drdx[ir].reshape([1, len(x)])
            else:
                dfdx = np.array(drdx).reshape([1, len(x)])
            fail = 0
            return dfdx, dgdx, fail

        def MinSensEqNorm(xNorm, f, g):
            x = denormalize(xNorm, xL, xU)
            dfxdx, dgdx, fail = MinSensEq(x, f, g)
            dfdx = dfxdx * (xU - xL)
            return dfdx, dgdx, fail

        def MaxSensEq(x, f, g):
            dfdx, dgdx, fail = MinSensEq(x, f, g)
            dfdx *= -1
            return dfdx, dgdx, fail

        def MaxSensEqNorm(xNorm, f, g):
            dfdx, dgdx, fail = MinSensEqNorm(xNorm, f, g)
            dfdx *= -1
            return dfdx, dgdx, fail

        def normalize(x, xL, xU):
            xNorm = (x - xL) / (xU - xL)
            return xNorm

        def denormalize(xNorm, xL, xU):
            x = (
                xNorm[
                    0 : np.size(xL),
                ]
                * (xU - xL)
                + xL
            )
            return x

        def DefineProb(
            x0min,
            x0max,
            xL,
            xU,
            OptModel,
            DesVarNorm,
        ):
            if DesVarNorm == 1:
                x0minnorm = normalize(x0min, xL, xU)
                x0maxnorm = normalize(x0max, xL, xU)
                MinProb = pyOpt.Optimization(
                    OptModel, MinSysEqNorm, obj_set=None
                )
                MaxProb = pyOpt.Optimization(
                    OptModel, MaxSysEqNorm, obj_set=None
                )
                if np.size(x0min) == 1:
                    MinProb.addVar(
                        'x', 'c', value=x0minnorm, lower=0.0, upper=1.0
                    )
                    MaxProb.addVar(
                        'x', 'c', value=x0maxnorm, lower=0.0, upper=1.0
                    )
                    self.nx = 1
                elif np.size(x0min) > 1:
                    for ii in range(np.size(x0min)):
                        MinProb.addVar(
                            'x' + str(ii + 1),
                            'c',
                            value=x0minnorm[ii],
                            lower=0.0,
                            upper=1.0,
                        )
                        MaxProb.addVar(
                            'x' + str(ii + 1),
                            'c',
                            value=x0maxnorm[ii],
                            lower=0.0,
                            upper=1.0,
                        )
                    self.nx = ii + 1
            elif DesVarNorm == 0:
                MinProb = pyOpt.Optimization(OptModel, MinSysEq)
                MaxProb = pyOpt.Optimization(OptModel, MaxSysEq)
                if np.size(x0min) == 1:
                    MinProb.addVar('x', 'c', value=x0min, lower=xL, upper=xU)
                    MaxProb.addVar('x', 'c', value=x0max, lower=xL, upper=xU)
                    self.nx = 1
                elif np.size(x0min) > 1:
                    for ii in range(np.size(x0min)):
                        MinProb.addVar(
                            'x' + str(ii + 1),
                            'c',
                            value=x0min[ii],
                            lower=xL[ii],
                            upper=xU[ii],
                        )
                        MaxProb.addVar(
                            'x' + str(ii + 1),
                            'c',
                            value=x0max[ii],
                            lower=xL[ii],
                            upper=xU[ii],
                        )
                    self.nx = ii + 1
            MinProb.addObj('f')
            MaxProb.addObj('f')
            if self.Alg == 'MMA':
                MinProb.addCon('g', type='i')
                MaxProb.addCon('g', type='i')
            return (MinProb, MaxProb)

        def AlphaLevelOpt(OptProb, OptAlg, SensCalc, Name):
            if self.Alg in [
                'MMA',
                'GCMMA',
                'CONMIN',
                'KSOPT',
                'SLSQP',
                'PSQP',
                'KSOPT',
                'SOLVOPT',
                'ALGENCAN',
                'NLPQLP',
            ]:
                if SensCalc == 'OptSensEq':
                    if self.paraNorm:
                        if Name[-3:] == 'Min':
                            fOpt, xOpt, info = OptAlg(
                                OptProb,
                                sens_type=MinSensEqNorm,
                                store_hst=Name,
                            )
                        elif Name[-3:] == 'Max':
                            fOpt, xOpt, info = OptAlg(
                                OptProb,
                                sens_type=MaxSensEqNorm,
                                store_hst=Name,
                            )
                    else:
                        if Name[-3:] == 'Min':
                            fOpt, xOpt, info = OptAlg(
                                OptProb, sens_type=MinSensEq, store_hst=Name
                            )
                        elif Name[-3:] == 'Max':
                            fOpt, xOpt, info = OptAlg(
                                OptProb, sens_type=MaxSensEq, store_hst=Name
                            )
                else:
                    fOpt, xOpt, info = OptAlg(
                        OptProb,
                        sens_type=SensCalc,
                        sens_step=self.deltax,
                        store_hst=Name,
                    )
            else:
                fOpt, xOpt, info = OptAlg(OptProb, store_hst=Name)
            return fOpt, xOpt

        def calcShadowUncertainty(Name, xOpt, xL, xU, ir, ialpha):
            if hasattr(self, 'ShadowUncertainty') is False:
                self.ShadowUncertainty = np.zeros(
                    (self.nr, self.np, self.nAlpha, 2)
                )
            for i in ['Min', 'Max']:
                Hist = pyOpt.History(Name + i, 'r')
                fNablaAll = Hist.read([0, -1], ['grad_obj'])[0]['grad_obj']
                fNabla = fNablaAll[-1]
                epsActive = 1e-3
                gL = xL - xOpt
                gU = xOpt - xU
                gLU = np.hstack((gL, gU))
                gLNabla = -np.eye(self.np)
                gUNabla = np.eye(self.np)
                gLUNabla = np.vstack((gLNabla, gUNabla))
                gLUActiveIndex = -gLU <= epsActive
                gLUNablaActive = gLUNabla[gLUActiveIndex, :]
                lam = np.zeros((self.np * 2))
                # lam[gLUActiveIndex] = (fNabla@pinv(gLUNablaActive)).T
                lam[gLUActiveIndex] = lstsq(
                    gLUNablaActive.T, -fNabla, rcond=None
                )[0]
                if self.paraNorm:
                    if np.size(xL) == 1:
                        denorm = np.array([xU[0] - xL[0], xU[0] - xL[0]])
                    else:
                        denorm = np.concatenate((xU - xL, xU - xL), axis=0)
                    lam /= denorm
                else:
                    lam
                # print(lam)
                lam = lam[: self.np] + lam[self.np :]
                # print(lam)
                if i == 'Max':
                    self.ShadowUncertainty[ir, :, ialpha, 1] = lam
                elif i == 'Min':
                    self.ShadowUncertainty[ir, :, ialpha, 0] = lam

        # Start of Optimization loop
        rUnc = np.zeros([self.nr, self.nAlpha, 2])
        if type(self.pUnc) is list:
            pUnc = np.zeros((len(self.pUnc), self.nAlpha, 2))
            for i, val in enumerate(self.pUnc):
                if type(val) == np.ndarray:
                    pUnc[i] = val
                elif type(val) == UncertainNumber:
                    pUnc[i] = val.Value
                elif type(self.pUnc) == UncertainNumber:
                    pUnc = self.pUnc.Value
            self.np = i + 1
        else:
            pUnc = np.zeros((1, self.nAlpha, 2))
            pUnc[0, :, :] = self.pUnc.Value
            self.np = 1
        ptilde = np.zeros([self.nr, np.size(pUnc, 0), self.nAlpha, 2])
        # SU = np.zeros([self.nr, np.size(pUnc, 0)*2, self.nAlpha, 2])
        lambdaR = np.zeros([self.nr, np.size(pUnc, 0) * 2, self.nAlpha, 2])
        for ir in range(self.nr):
            for ialpha in reversed(range(self.nAlpha)):
                xL = pUnc[:, ialpha, 0]
                xU = pUnc[:, ialpha, 1]
                if abs(np.sum(abs(xU - xL) / (xL + np.spacing(1)))) < 0.001:
                    f, g, fail = MinSysEq(xL)
                    rUnc[ir, ialpha, 0] = f
                    rUnc[ir, ialpha, 1] = f
                    pMin = xU
                    pMax = xU
                else:
                    if ialpha == self.nAlpha - 1:
                        x0min = (xU + xL) / 2
                        x0max = x0min
                    else:
                        for j, val in enumerate(pMin):
                            if (
                                abs(val - pUnc[j, ialpha + 1, 0])
                                / (pMax[j] + np.spacing(1))
                                < 0.01
                            ):
                                x0min[j] = pUnc[j, ialpha, 0]
                            if (
                                abs(val - pUnc[j, ialpha + 1, 1])
                                / (pMax[j] + np.spacing(1))
                                < 0.01
                            ):
                                x0min[j] = pUnc[j, ialpha, 1]
                        for j, val in enumerate(pMax):
                            if (
                                abs(val - pUnc[j, ialpha + 1, 0])
                                / (pMax[j] + np.spacing(1))
                                < 0.01
                            ):
                                x0max[j] = pUnc[j, ialpha, 0]
                            if (
                                abs(val - pUnc[j, ialpha + 1, 1])
                                / (pMax[j] + np.spacing(1))
                                < 0.01
                            ):
                                x0max[j] = pUnc[j, ialpha, 1]
                    Name = 'alphaOpt_p' + str(ir + 1) + '_alpha' + str(ialpha)
                    if self.Alg == 'CMAES':
                        ResMin = cma.CMAEvolutionStrategy(
                            x0min, 0.5, {'bounds': [xL, xU]}
                        ).optimize(
                            MinSysEq,
                            min_iterations=25,
                            iterations=100,
                            verb_disp=0,
                        )
                        ResMax = cma.CMAEvolutionStrategy(
                            x0max, 0.5, {'bounds': [xL, xU]}
                        ).optimize(
                            MaxSysEq,
                            min_iterations=25,
                            iterations=100,
                            verb_disp=0,
                        )
                        rMin = ResMin[1]
                        pMin = ResMin[0]
                        rMax = ResMax[1]
                        pMax = ResMax[0]
                    elif self.Alg == '2Step':
                        self.Alg = 'ALHSO'
                        alphaLevelOptAlg = eval('pyOpt.' + self.Alg + '()')
                        alphaLevelOptAlg = AlgOptions(
                            alphaLevelOptAlg, self.Alg, self.epsStop
                        )
                        alphaLevelOptAlg.setOption('maxoutiter', 10)
                        MinProb, MaxProb = DefineProb(
                            x0min,
                            x0max,
                            xL,
                            xU,
                            self.UncertainModel,
                            self.paraNorm,
                        )
                        rMin, pMin = AlphaLevelOpt(
                            MinProb,
                            alphaLevelOptAlg,
                            self.SensCalc,
                            Name + 'Min',
                        )
                        rMax, pMax = AlphaLevelOpt(
                            MaxProb,
                            alphaLevelOptAlg,
                            self.SensCalc,
                            Name + 'Max',
                        )
                        self.Alg = 'NLPQLP'
                        alphaLevelOptAlg = eval('pyOpt.' + self.Alg + '()')
                        alphaLevelOptAlg = AlgOptions(
                            alphaLevelOptAlg, self.Alg, self.epsStop
                        )
                        MinProb, MaxProb = DefineProb(
                            pMin,
                            pMax,
                            xL,
                            xU,
                            self.UncertainModel,
                            self.paraNorm,
                        )
                        rMin, pMin = AlphaLevelOpt(
                            MinProb,
                            alphaLevelOptAlg,
                            self.SensCalc,
                            Name + 'Min',
                        )
                        rMax, pMax = AlphaLevelOpt(
                            MaxProb,
                            alphaLevelOptAlg,
                            self.SensCalc,
                            Name + 'Max',
                        )
                    elif self.Alg[:5] == 'PyGMO':
                        ngen = 50
                        nindiv = 10
                        if self.Alg == 'PyGMO_monte_carlo':
                            algo = pg.algorithm(
                                pg.monte_carlo(iters=AlgOptions.iter)
                            )
                        elif self.Alg == 'PyGMO_de':
                            algo = pg.algorithm(
                                pg.de(
                                    gen=ngen,
                                    F=1,
                                    CR=1,
                                    variant=2,
                                    ftol=1e-03,
                                    xtol=1e-03,
                                )
                            )
                        else:
                            algo = eval(
                                'pg.algorithm(pg.'
                                + self.Alg[6:]
                                + '('
                                + str(AlgOptions.gen)
                                + '))'
                            )

                        class minFn:
                            def fitness(self, x):
                                f, g, fail = MinSysEq(x)
                                return np.array([f])

                            def get_bounds(self):
                                return (xL, xU)

                        class maxFn:
                            def fitness(self, x):
                                f, g, fail = MaxSysEq(x)
                                return np.array([f])

                            def get_bounds(self):
                                return (xL, xU)

                        probMin = pg.problem(minFn())
                        popMin = pg.population(probMin, nindiv)
                        popMin = algo.evolve(popMin)
                        pMin = popMin.champion_x
                        rMin = popMin.champion_f[0]
                        probMax = pg.problem(maxFn())
                        popMax = pg.population(prob=probMax, size=nindiv)
                        popMax = algo.evolve(popMax)
                        pMax = popMax.champion_x
                        rMax = popMax.champion_f[0]
                    else:
                        alphaLevelOptAlg = eval('pyOpt.' + self.Alg + '()')
                        alphaLevelOptAlg = AlgOptions(
                            alphaLevelOptAlg, self.Alg, self.epsStop
                        )
                        MinProb, MaxProb = DefineProb(
                            x0min,
                            x0max,
                            xL,
                            xU,
                            self.UncertainModel,
                            self.paraNorm,
                        )
                        rMin, pMin = AlphaLevelOpt(
                            MinProb,
                            alphaLevelOptAlg,
                            self.SensCalc,
                            Name + 'Min',
                        )
                        rMax, pMax = AlphaLevelOpt(
                            MaxProb,
                            alphaLevelOptAlg,
                            self.SensCalc,
                            Name + 'Max',
                        )

                    pMin = pMin[0 : np.size(xL)]
                    pMax = pMax[0 : np.size(xL)]
                    pMin = np.resize(
                        pMin,
                        [
                            np.size(xL),
                        ],
                    )
                    pMax = np.resize(
                        pMax,
                        [
                            np.size(xL),
                        ],
                    )
                    if self.paraNorm:
                        pMinNorm = pMin
                        pMaxNorm = pMax
                        pMin = denormalize(pMinNorm, xL, xU)
                        pMax = denormalize(pMaxNorm, xL, xU)
                    rMax = -rMax
                    rUnc[ir, ialpha, 0] = rMin
                    rUnc[ir, ialpha, 1] = rMax
                    ptilde[ir, :, ialpha, 0] = pMin
                    ptilde[ir, :, ialpha, 1] = pMax
                    if self.Alg in ['NLPQLP', 'SLSQP', 'MMA']:
                        calcShadowUncertainty(Name, pMin, xL, xU, ir, ialpha)
        OutputData = {}
        OutputData['pUnc'] = ptilde
        OutputData['nEval'] = self.nEval
        OutputData['lambdaR'] = lambdaR
        # self.SU = SU
        if self.nr > 1:
            self.rUnc = [[]] * len(rUnc)
            for i, val in enumerate(rUnc):
                if self.nAlpha == 1:
                    self.rUnc[i] = UncertainNumber(
                        val, Form='interval', nalpha=self.nAlpha
                    )
                else:
                    self.rUnc[i] = UncertainNumber(
                        val, Form='empirical', nalpha=self.nAlpha
                    )
        else:
            if self.nAlpha == 1:
                self.rUnc = UncertainNumber(
                    rUnc[0][0], Form='interval', nalpha=self.nAlpha
                )
            else:
                self.rUnc = UncertainNumber(
                    rUnc[0], Form='empirical', nalpha=self.nAlpha
                )
        self.OutputData = OutputData

    def calcRobustness(self):
        self.pAreaSum = 0
        self.pAreaNormSum = 0
        self.rAreaSum = 0
        self.rAreaNormSum = 0
        if type(self.pUnc) is list:
            for val in self.pUnc:
                val.calcArea()
                self.pAreaSum += val.Area
                self.pAreaNormSum += val.AreaNorm
        elif type(self.pUnc) == UncertainNumber:
            self.pUnc.calcArea()
            self.pAreaSum = self.pUnc.Area
            self.pAreaNormSum = self.pUnc.AreaNorm
        if type(self.rUnc) is list:
            for val in self.rUnc:
                val.calcArea()
                self.rAreaSum += val.Area
                self.rAreaNormSum += val.AreaNorm
        elif type(self.rUnc) == UncertainNumber:
            self.rUnc.calcArea()
            self.rAreaSum = self.rUnc.Area
            self.rAreaNormSum = self.rUnc.AreaNorm
        self.SystemRobustness = np.divide(self.pAreaSum, self.rAreaSum)
        self.SystemRobustnessNorm = np.divide(
            self.pAreaNormSum, self.rAreaNormSum
        )


def ShadowUncertaintyPrices(SP, SU):
    SUP = np.abs(SP * SU)
    return SUP


if __name__ == '__main__':
    printSplash()
    print('Start FuzzAnPy from file containing system equations!')
    print('See documentation for further help.')
    print()
    print()
    print('Quick test 1: pUnc - pUnc = rUnc (should be 0)')
    print('-' * 70)

    def SysEq1(p, x):
        return p - p

    pInt = UncertainNumber([1, 5])
    UncertainProblem = UncertainAnalysis()
    UncertainProblem.pUnc = pInt
    UncertainProblem.SysEq = SysEq1
    UncertainProblem.calculate()
    print('rUnc = pUnc - pUnc')
    print('rUnc = ' + str(UncertainProblem.rUnc.Value))
    UncertainProblem.calcRobustness()
    print()
    print('System robustness:')
    print(UncertainProblem.SystemRobustness)
    print()
    print()
    print('Quick test 2: pUnc1 * pUnc2 = rUnc and check shadow uncertainty')
    print('-' * 70)

    def SysEq2(p, x):
        return p[0] * p[1]

    pInt1 = UncertainNumber([1, 5])
    pInt2 = UncertainNumber([1, 5])
    UncertainProblem = UncertainAnalysis()
    UncertainProblem.pUnc = [pInt1, pInt2]
    UncertainProblem.SysEq = SysEq2
    UncertainProblem.calculate()
    print('rUnc = pUnc - pUnc')
    print('rUnc = ' + str(UncertainProblem.rUnc.Value))
    UncertainProblem.calcRobustness()
    print()
    print('Shadow uncertainty:')
    print(UncertainProblem.ShadowUncertainty[0, :, :, :])
    print()
    print('Reset pUnc1 to change upper bound -1:')
    pInt1 = UncertainNumber([1, 5])
    UncertainProblem.pUnc = [pInt1, pInt2]
    UncertainProblem.calculate()
    print('New uncertain respose rUnc:')
    print('rUnc = ' + str(UncertainProblem.rUnc.Value))
    print()
    print()
    print('Quick test 3: Rosenbock banana function with different algorithms')
    print('-' * 70)
    from scipy.optimize import rosen, rosen_der

    def SysEq3(p, x):
        return rosen(p)

    def SensEq3(p, r, g, x):
        return rosen_der(p)

    nAlpha = 3
    pFuzz = [[]] * 2
    pFuzz[0] = UncertainNumber([1, 2, 3, 4], Form='trapazoid', nalpha=nAlpha)
    pFuzz[1] = UncertainNumber([1, 2, 3, 4], Form='trapazoid', nalpha=nAlpha)
    Prob = UncertainAnalysis(SysEq2, pUnc=pFuzz, SensEq=SensEq3)
    Prob.Alg = 'PyGMO_de'
    Prob.nAlpha = nAlpha
    Prob.paraNorm = 0
    Prob.epsStop = 1e-6
    Prob.para = 1
    Prob.calculate()
    print('Result with stochastic algorithm:')
    print('nEval = ' + str(Prob.nEval))
    print('rUnc = \n' + str(Prob.rUnc.Value))
    print()
    Prob.Alg = 'NLPQLP'
    Prob.nEval = 0
    Prob.calculate()
    print('Result with gradient-based algorithm:')
    print('nEval = ' + str(Prob.nEval))
    print('rUnc = \n' + str(Prob.rUnc.Value))
    print()
    print(
        'Result with two-step methods stochastic and gradient-based algorithms:'
    )
    Prob.Alg = '2Step'
    Prob.nEval = 0
    Prob.calculate()
    print('nEval = ' + str(Prob.nEval))
    print('rUnc = \n' + str(Prob.rUnc.Value))
