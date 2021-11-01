import pyUngewiss as pu
import numpy as np
import unittest


class TestMethods(unittest.TestCase):
    def test_FuzzyNumber(self):
        nAlpha = 3
        pFuzz = pu.UncertainNumber(
            [1, 2, 3, 4], Form='trapazoid', nalpha=nAlpha
        )
        pFuzzTarget = np.array([[2.0, 3.0], [1.5, 3.5], [1.0, 4.0]])
        self.assertAlmostEqual(pFuzz.Value.tolist(), pFuzzTarget.tolist())

    def test_IntervalAnalysis(self):
        pInt = pu.UncertainNumber([1, 5])

        def SysEq1(p, x):
            return p - p

        UncertainProblem = pu.UncertainAnalysis()
        UncertainProblem.pUnc = pInt
        UncertainProblem.SysEq = SysEq1
        UncertainProblem.calculate()
        rUncTarget = np.array([[0.0, 0.0]])
        self.assertAlmostEqual(
            UncertainProblem.rUnc.Value.tolist(), rUncTarget.tolist()
        )

    def test_FuzzyAnalysis(self):
        from scipy.optimize import rosen, rosen_der

        def SysEq2(p, x):
            return rosen(p)

        def SensEq2(p, r, g, x):
            return rosen_der(p)

        pFuzz1 = pu.UncertainNumber([1, 2, 3, 4], Form='trapazoid', nalpha=3)
        pFuzz2 = pu.UncertainNumber([1, 2, 3, 4], Form='trapazoid', nalpha=3)
        pFuzz = [pFuzz1, pFuzz2]
        Prob = pu.UncertainAnalysis(SysEq2, pUnc=pFuzz, SensEq=SensEq2)
        Prob.Alg = 'NLPQLP'
        Prob.nAlpha = 3
        Prob.paraNorm = 0
        Prob.epsStop = 1e-6
        Prob.para = 1
        Prob.calculate()
        rUncTarget = np.array(
            [
                [1.01000000e02, 4.90400000e03],
                [2.50000000e-01, 1.15625000e04],
                [0, 2.25090000e04],
            ]
        )
        # self.assertAlmostEqual(Prob.rUnc.Value.tolist(),
        #                       rUncTarget.tolist())
        self.assertAlmostEqual(
            Prob.rUnc.Value.tolist()[0], rUncTarget.tolist()[0]
        )
        self.assertAlmostEqual(
            Prob.rUnc.Value.tolist()[1], rUncTarget.tolist()[1]
        )
        self.assertAlmostEqual(
            Prob.rUnc.Value.tolist()[2][0], rUncTarget.tolist()[2][0]
        )
        self.assertAlmostEqual(
            Prob.rUnc.Value.tolist()[2][1], rUncTarget.tolist()[2][1]
        )

    def test_Robustness(self):
        UncertainProblem = pu.UncertainAnalysis()
        UncertainProblem.pUnc = pu.UncertainNumber([1, 2])
        UncertainProblem.rUnc = pu.UncertainNumber([3, 6])
        UncertainProblem.calcRobustness()
        self.assertAlmostEqual(UncertainProblem.SystemRobustness, 1 / 3)


if __name__ == '__main__':
    unittest.main()
