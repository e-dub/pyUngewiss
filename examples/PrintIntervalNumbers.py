import pyUngewiss as pu
import numpy as np

print('Test of interval plotting:')
print('')
data = np.array(
    [
        [300, 900],
        [350, 850],
        [400, 800],
        [450, 750],
        [500, 700],
        [550, 650],
        [599, 601],
        [300, 900],
        [350, 850],
        [400, 800],
        [450, 750],
        [500, 700],
        [550, 650],
        [599, 601],
        [300, 900],
        [350, 850],
        [400, 800],
        [450, 750],
        [500, 700],
        [550, 650],
        [350, 850],
        [400, 800],
        [450, 750],
        [500, 700],
        [550, 650],
        [599, 601],
        [300, 900],
        [350, 850],
        [400, 800],
        [450, 750],
        [500, 700],
        [550, 650],
        [599, 601],
        [300, 900],
        [350, 850],
        [400, 800],
        [450, 750],
        [500, 700],
        [550, 650],
    ]
)
AIntArray = [[]] * len(data)
for i, val in enumerate(data):
    AIntArray[i] = pu.UncertainNumber((val))
labels = []
for ii in range(len(data)):
    labels.append('Stiffness ' + str(ii + 1))
ax = pu.plotIntervals(AIntArray, labels=labels, Units='N/mm')
data = np.array(
    [
        [300, 900],
        [400, 800],
        [450, 750],
        [500, 700],
        [550, 650],
        [599, 601],
        [300, 900],
    ]
)
labels = []
for ii in range(len(data)):
    labels.append('Stiffness ' + str(ii + 1))
ax = pu.plotIntervals(data, labels=labels, Units='N/mm')
data = np.array([[300, 900], [400, 800]])
labels = []
for ii in range(len(data)):
    labels.append('Stiffness ' + str(ii + 1))
ax = pu.plotIntervals(data, labels=labels, Units='N/mm')
