"""
Install pyre via
pip3 install git+git://github.com/hackl/pyre.git
"""
import numpy as np
from scipy.linalg import solve
import pyre as pyre

def K(x):
    """ Stiffness matrix """
    tmp = np.zeros((2,2),np.float64)
    tmp[0,0] = x[0]+x[1]+x[2]+x[3]
    tmp[0,1] = -x[2]-x[3]
    tmp[1,0] = -x[2]-x[3]
    tmp[1,1] = x[2]+x[3]
    return tmp

def f(x):
    """ Load vector """
    tmp = np.zeros(2,np.float64)
    tmp[0] = x[4]
    tmp[1] = x[5]
    return tmp

def example_limitstatefunction(X1, X2, X3, X4, X5, X6):
    x0 = np.zeros(6)
    x0[0] = X1
    x0[1] = X2
    x0[2] = X3
    x0[3] = X4
    x0[4] = X5
    x0[5] = X6
    # Displacements
    u = np.linalg.solve(K(x0),f(x0))
    return 2.17-u[1]


limit_state = pyre.LimitState(example_limitstatefunction)
options = pyre.AnalysisOptions()
stochastic_model = pyre.StochasticModel()
stochastic_model.addVariable(pyre.Lognormal('X1', 10., 2.))
stochastic_model.addVariable(pyre.Lognormal('X2', 10., 2.))
stochastic_model.addVariable(pyre.Lognormal('X3', 10., 2.))
stochastic_model.addVariable(pyre.Lognormal('X4', 10., 2.))
stochastic_model.addVariable(pyre.Lognormal('X5', 10., 2.))
stochastic_model.addVariable(pyre.Lognormal('X6', 10., 2.))
Analysis = pyre.Form(analysis_options=options,
                     stochastic_model=stochastic_model,
                     limit_state=limit_state)

