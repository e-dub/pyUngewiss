# -*- coding: utf-8 -*-
"""
A class and operations for uncertain values using intervals and fuzzy numbers

Previously called FuzzyNumber.py, but now renamed to current name to represent
concentration not only with fuzzy numbers, making intervals default.

All equations refer to the following work:
Wehrle (2015) Design optimization of lightweight space-frame structures
considering crashworthiness and parameter uncertainty

Area of uncertain value: eq. (4.14)
Robustness: eq. (4.15)-(4.17)
Normalize of uncertain value: eq. (4.18)
Membership functions: eq. (4.21)-(4.24)
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import rcParams
import seaborn as sns

np.seterr(divide='ignore', invalid='ignore')
DefaultColor = (0, 123 / 255, 228 / 255)
DefaultEdgeColor = (0, 123 / 255, 228 / 255)
pgf_with_pdflatex = {
    'pgf.texsystem': 'pdflatex',
    'pgf.preamble': [
        r'\usepackage[utf8x]{inputenc}',
        r'\usepackage[T1]{fontenc}',
        r'\usepackage{cmbright}',
    ],
}


class UncertainNumber(object):
    def __init__(self, values, Form='interval', nalpha=1):
        self.Form = Form
        self.nalpha = nalpha
        self.type = 'UncertainNumber'
        if self.Form == 'interval':
            values = np.array(values).reshape(
                2,
            )
            alpha = np.linspace(0, 1.0, self.nalpha)
            A = values[0] - (0) * alpha
            B = values[1] - (0) * alpha
            self.Value = np.array([A, B]).T
        elif self.Form == 'triangle':
            alpha = np.linspace(0, 1, self.nalpha)
            A = values[1] - (values[1] - values[0]) * alpha
            B = values[1] + (values[2] - values[1]) * alpha
            self.Value = np.array([A, B]).T
        elif self.Form == 'trapazoid':
            alpha = np.linspace(0, 1.0, self.nalpha)
            if values[1] == values[0]:
                A = values[1] - (0) * alpha
            else:
                A = values[1] - (values[1] - values[0]) * alpha
            if values[2] == values[3]:
                B = values[3] + (0) * alpha
            else:
                B = values[2] + (values[3] - values[2]) * alpha
            self.Value = np.array([A, B]).T
        elif self.Form == 'gauss-cuttoff':
            mean = values[0]
            sigmaleft = values[1]
            sigmaright = values[2]
            sigmatrunc = values[3]
            alpha = np.linspace(0, 1, self.nalpha)
            A = np.zeros([self.nalpha, 2])
            A[:, 0] = np.flip(
                mean - np.sqrt(-2 * sigmaleft ** 2 * np.log(alpha))
            )
            A[:, 1] = np.flip(
                mean + np.sqrt(-2 * sigmaleft ** 2 * np.log(alpha))
            )
            A[A < mean - sigmatrunc * sigmaleft] = (
                mean - sigmatrunc * sigmaleft
            )
            A[A > mean + sigmatrunc * sigmaright] = (
                mean + sigmatrunc * sigmaright
            )
            self.Value = A
        elif self.Form == 'empirical':
            self.Value = values
        else:
            raise ValueError('Unexpected membership function form:')

    def normalizeValue(self):
        self.ValueNorm = np.zeros(np.shape(self.Value))
        for i, v in enumerate(self.Value):
            self.ValueNorm[i, :] = v / ((v[1] + v[0]) / 2)

    def calcArea(self):
        self.normalizeValue()
        self.Area = np.zeros([np.size(self.Value, 0), 1])
        self.AreaNorm = np.zeros([np.size(self.Value, 0), 1])
        mu1 = np.flip(np.linspace(1, 0, self.nalpha))
        mu2 = np.linspace(1, 0, self.nalpha)
        ymu = [mu1, mu2]
        ymu = np.resize(
            ymu,
            [
                self.nalpha * 2,
            ],
        )
        xVal = [self.Value[:, 0], self.Value[:, 1]] + np.min(self.Value[:, :])
        xVal = np.resize(xVal, [self.nalpha * 2])
        self.Area = np.abs(np.trapz(y=ymu, x=xVal))
        if self.Area == 0.0:
            self.AreaNorm = 0.0
        else:
            xValNorm = [self.ValueNorm[:, 0], self.ValueNorm[:, 1]]
            xValNorm = np.resize(
                xValNorm,
                [
                    self.nalpha * 2,
                ],
            )
            self.AreaNorm = np.abs(np.trapz(y=ymu, x=xValNorm))

    def printValue(self):
        for ai in self.Value:
            print('[{},  {}]'.format(ai[0], ai[1]))

    def printValueNorm(self):
        for ai in self.ValueNorm:
            print('[{},  {}]'.format(ai[0], ai[1]))

    def plotValue(
        self,
        xlabel=[],
        ylabel=[],
        pdpi=100,
        fill=True,
        fontsize=10,
        xsize=2.5,
        ysize=2,
        xlimits=[],
        lang='EN',
        TextRender=True,
        AxisNameLong=True,
        frame=False,
        nYTicks=2,
        nXTicks=5,
        xAxisRot=False,
        trans=1.0,
        font='tex gyre pagella',
        Units=[],
        plotBuffer=0.1,
        color=DefaultColor,
        delta=0.5,
        padLabel=[],
    ):
        plt.rcParams['font.family'] = font
        rcParams.update({'figure.autolayout': True})
        rcParams.update(pgf_with_pdflatex)
        if not TextRender:
            rcParams.update({'svg.fonttype': 'none'})
        if self.Form == 'interval' or self.nalpha == 1:   # for interval plot
            xsize = 5
            ysize = 0.25
            fig = plt.figure(figsize=(xsize, ysize), dpi=pdpi)
            ax = fig.add_subplot(1, 1, 1)
            thick = 1 ** 0.1 * 0.3
            yplaces = [0.5 + i for i in reversed(range(1))]
            ax.set_yticks(yplaces)
            if not xlabel:
                xlabel = r'Uncertain value'
            if not ylabel:
                ylabel = r'Parameter'
            ax.set_yticklabels([ylabel], horizontalalignment='left')
            ax.set_ylim((0, 1))
            start = self.Value[0, 0]
            end = self.Value[0, 1]
            pos = yplaces[0]
            ax.add_patch(
                patches.Rectangle(
                    (start, pos - delta / 2.0), end - start, thick, color=color
                )
            )
            pmin = np.min(self.Value)
            pmax = np.max(self.Value)
            ax.plot(
                (pmin * (1 - plotBuffer), pmax * (1 + plotBuffer)),
                (0, 0),
                'w',
                alpha=0.0,
            )
            sns.despine()
            nlmax = len(ylabel)
            yax = ax.get_yaxis()
            if not padLabel:
                padLabel = min(nlmax * 5.5, 175)
            yax.set_tick_params(direction='out', pad=padLabel)
            if not Units:
                ax.set_xlabel(xlabel)
            else:
                ax.set_xlabel(xlabel + ' [' + Units + ']')
            if xAxisRot:
                plt.setp(ax.get_xticklabels(), rotation='vertical')
            plt.xlim(
                (
                    np.min(self.Value)
                    - (np.max(self.Value) - np.min(self.Value)) / 5
                ),
                (
                    np.max(self.Value)
                    + (np.max(self.Value) - np.min(self.Value)) / 5
                ),
            )
        else:
            fig = plt.figure(figsize=(xsize, ysize), dpi=pdpi)
            ax = fig.add_subplot(1, 1, 1)
            mu = np.linspace(1, 0, np.size(self.Value, 0))
            xplot = np.concatenate(
                (np.flipud(self.Value[:, 0]), self.Value[:, 1]), axis=0
            )
            muplot = np.concatenate((np.flipud(mu), mu), axis=0)
            if fill:
                ax.fill_between(
                    self.Value[:, 0],
                    mu,
                    facecolor=color,
                    alpha=trans,
                    linewidth=0.0,
                )
                ax.fill_between(
                    [self.Value[0, 1], self.Value[0, 0]],
                    [1, 1],
                    facecolor=color,
                    alpha=trans,
                    linewidth=0.0,
                )
                ax.fill_between(
                    self.Value[:, 1],
                    mu,
                    facecolor=color,
                    alpha=trans,
                    linewidth=0.0,
                )
            ax.plot(xplot, muplot, 'k-', linewidth=1)
            plt.ylim(0, 1.10)
            if xlimits != []:
                plt.xlim(xlimits[0], xlimits[1])
            plt.xlabel(xlabel)
            plt.yticks(np.arange(0, 1 + (nYTicks - 1), 1 / (nYTicks - 1)))
            xloc = plt.MaxNLocator(nXTicks)
            # plt.MaxNLocator(nXTicks)
            xloc = plt.AutoLocator()
            # xloc = plt.LinearLocator(numticks=nXTicks)
            ax.xaxis.set_major_locator(xloc)
            if not frame:
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)
                ax.get_xaxis().tick_bottom()
                ax.get_yaxis().tick_left()
            if xAxisRot:
                plt.setp(ax.get_xticklabels(), rotation='vertical')
            if AxisNameLong:
                if lang == 'EN':
                    if TextRender:
                        plt.ylabel(u'Membership $\\mu$')
                    elif not TextRender:
                        plt.ylabel(r'Membership \$\\mu\$')
                elif lang == 'DE':
                    if TextRender:
                        plt.ylabel(u'Zugehörigkeit $\\mu$')
                    elif not TextRender:
                        plt.ylabel(r'Zugeh\"origkeit \$\\mu\$')
            else:
                if TextRender:
                    h = plt.ylabel(u'$\\mu$')
                elif not TextRender:
                    h = plt.ylabel(r'\$\\mu\$')
                h.set_rotation(0)
            if xlabel == []:
                if lang == 'EN':
                    plt.xlabel(r'Value')
                elif lang == 'DE':
                    plt.xlabel(r'Wert')
            rcParams.update({'font.size': fontsize})
            # plt.tight_layout()
            self.Plot = plt


def plotIntervals(
    pUncList,
    labels=[],
    Units=[],
    show=False,
    xPlot=5,
    color=DefaultColor,
    delta=0.5,
    plotBuffer=0.1,
    xAxisRot=False,
    xlabel='Uncertain value',
    font='tex gyre pagella',
    padLabel=[],
):
    if isinstance(pUncList, list):
        data = [[]] * len(pUncList)
        for i, val in enumerate(pUncList):
            if type(val) == np.ndarray:
                data[i] = val
            elif type(val) == UncertainNumber:
                data[i] = val.Value
        if isinstance(color, list) is False:
            color = [color] * len(pUncList)
    elif type(pUncList) == UncertainNumber:
        data = pUncList.Value
    else:
        data = [pUncList]
        color = [color]
    nP = len(data)
    if np.shape(np.shape(data))[0] == 1:
        data = data.reshape(nP, 1)
        interval = False
    elif np.shape(np.shape(data))[0] == 3:
        interval = True
    else:
        interval = True
    yPlot = nP * 0.25
    thick = nP ** 0.1 * 0.3
    plt.rcParams['font.family'] = font
    plt.rcParams['figure.figsize'] = (xPlot, yPlot)
    if not labels:
        if np.shape(data)[0] == 1:
            labels = ['Parameter']
        else:
            for ii in range(np.shape(data)[0]):
                labels.append('Parameter ' + str(ii + 1))
    yplaces = [0.5 + i for i in reversed(range(nP))]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yticks(yplaces)
    ax.set_yticklabels(labels, horizontalalignment='left')
    ax.set_ylim((0, nP))
    for ii, val in enumerate(data):
        if interval:
            start = val[0, 0]
            end = val[0, 1]
            pos = yplaces[ii]
            ax.add_patch(
                patches.Rectangle(
                    (start, pos - delta / 2.0),
                    end - start,
                    thick,
                    color=color[ii],
                )
            )
        elif not interval:
            ax.plot(data[ii, 0], ii + 0.5, 'x', color=color[ii])
    pmin = np.min(data)
    pmax = np.max(data)
    ax.plot(
        (pmin * (1 - plotBuffer), pmax * (1 + plotBuffer)),
        (0, 0),
        'w',
        alpha=0.0,
    )
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.yaxis.set_ticks_position('none')
    ax.xaxis.set_ticks_position('bottom')
    plt.grid(False)
    nlmax = 0
    for ilabel in labels:
        nl = len(ilabel)
        if nl > nlmax:
            nlmax = nl
    yax = ax.get_yaxis()
    if not padLabel:
        padLabel = min(nlmax * 5.5, 175)
    yax.set_tick_params(direction='out', pad=padLabel)
    if not Units:
        ax.set_xlabel(xlabel)
    else:
        ax.set_xlabel(xlabel + ' [' + Units + ']')
    if xAxisRot:
        plt.setp(ax.get_xticklabels(), rotation='vertical')
    plt.xlim(
        np.min(data) - (np.max(data) - np.min(data)) / 5,
        np.max(data) + (np.max(data) - np.min(data)) / 5,
    )
    return plt, ax


def plotUncertainFn(
    pUncList,
    x=[],
    xlabel='Dependent parameter',
    ylabel='Value of uncertain function',
    color=DefaultColor,
    fontsize=12,
    font='tex gyre pagella',
    pdpi=100,
    xlimits=[],
    ylimits=[],
    xsize=7,
    ysize=5,
    outline=False,
    fill=True,
    nYTicks=5,
    nXTicks=5,
    xAxisRot=False,
    frame=True,
):
    if type(pUncList) is list:
        if len(pUncList) == 0:
            nAlpha = pUncList[0].nalpha
        else:
            nAlpha = pUncList[0].nalpha
        pUnc = np.zeros((len(pUncList), nAlpha, 2))
        for i, val in enumerate(pUncList):
            if type(val) == np.ndarray:
                pUnc[i] = val
            elif type(val) == UncertainNumber:
                pUnc[i] = val.Value
            elif type(pUncList) == UncertainNumber:
                pUnc = pUncList.Value
            elif type(val) == list:
                pUnc[i] = val[0].Value
    else:
        pUnc = np.zeros((1, nAlpha, 2))
        pUnc[0, :, :] = pUncList.Value
    rFuzz = pUnc
    plt.rcParams['font.family'] = font
    nalpha = np.size(rFuzz, 1)
    if len(x) == 0:
        x = np.linspace(0, np.size(rFuzz))
    fig1 = plt.figure(figsize=(xsize, ysize), dpi=pdpi)
    ax1 = fig1.add_subplot(1, 1, 1)
    if fill:
        for ii in reversed(range(np.size(rFuzz, 1) - 1)):
            ax1.fill_between(
                x,
                rFuzz[:, ii + 1, 1],
                rFuzz[:, ii, 1],
                facecolor=color,
                alpha=1.0 / (nalpha) * (nalpha - ii - 1) * 0.9,
                linewidth=0,
                edgecolor=DefaultEdgeColor,
            )
            ax1.fill_between(
                x,
                rFuzz[:, ii + 1, 0],
                rFuzz[:, ii, 0],
                facecolor=color,
                alpha=1.0 / (nalpha) * (nalpha - ii - 1) * 0.9,
                linewidth=0,
                edgecolor=DefaultEdgeColor,
            )
        ax1.fill_between(
            x,
            rFuzz[:, 0, 1],
            rFuzz[:, 0, 0],
            facecolor=color,
            alpha=1.0,
            linewidth=0,
        )
    if outline:
        ax1.plot(
            x, rFuzz[:, :, 0], x, rFuzz[:, :, 1], color=color, linewidth=1.0
        )
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if not frame:
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax1.get_xaxis().tick_bottom()
        ax1.get_yaxis().tick_left()
    if xlimits != [] and xlimits != []:
        plt.yticks(np.linspace(ylimits[0], ylimits[1], nYTicks))
    xloc = plt.MaxNLocator(nXTicks)
    ax1.xaxis.set_major_locator(xloc)
    if xAxisRot:
        plt.setp(ax1.get_xticklabels(), rotation='vertical')

    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    # ax1.spines['left'].set_position(('outward', 17))
    # ax1.spines['bottom'].set_position(('outward', 17))
    if xlimits != []:
        plt.xlim(xlimits[0], xlimits[1])
    if ylimits != []:
        plt.ylim(ylimits[0], ylimits[1])
    return plt, ax1


if __name__ == '__main__':
    print('Test printing of interval number:')
    AInt = UncertainNumber([10, 20])
    AInt.normalizeValue()
    AInt.calcArea()
    AInt.printValue()
    AInt.printValueNorm()
    AInt.plotValue()
    plotIntervals(AInt.Value)
    plotIntervals([AInt.Value, AInt.Value / 10], color=['r', 'b'])
    print()

    print('Test printing of trapazoidal fuzzy number:')
    ATrap = UncertainNumber([10, 20, 30, 40], Form='trapazoid', nalpha=6)
    ATrap.normalizeValue()
    ATrap.printValue()
    ATrap.printValueNorm()
    print()

    print('Test plotting all uncertain numbers:')
    ATrapSing = UncertainNumber([25, 25, 25, 25], Form='trapazoid', nalpha=6)
    ATrapSing.plotValue()
    ATrapInt = UncertainNumber([10, 10, 40, 40], Form='trapazoid', nalpha=6)
    ATrapInt.plotValue()
    ATrap = UncertainNumber([10, 20, 30, 40], Form='trapazoid', nalpha=6)
    ATrap.plotValue()
    ATri = UncertainNumber([10, 25, 40], Form='triangle', nalpha=6)
    ATri.plotValue()
    AGauss = UncertainNumber([25, 5, 5, 1.5], Form='gauss-cuttoff', nalpha=61)
    AGauss.plotValue()
    print()

    print('Test plotting several intervals in same plot:')
    m = UncertainNumber([2.0, 2.5])
    k = UncertainNumber([40000, 60000])
    pUnc = [m, k]
    plotIntervals([m, k], labels=['mass $m$ [kg]', 'stiffness $k$ [N/mm]'])
