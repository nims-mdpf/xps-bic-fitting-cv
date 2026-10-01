import os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

class FixedOrderFormatter(ScalarFormatter):
    def __init__(self, order_of_mag=0, useOffset=True, useMathText=True):
        self._order_of_mag = order_of_mag
        ScalarFormatter.__init__(self, useOffset=useOffset, useMathText=useMathText)
    def _set_orderOfMagnitude(self, range):
        self.orderOfMagnitude = self._order_of_mag

def fittingPlot(path):

    dat = np.loadtxt(path, skiprows=1, delimiter=',')
    fname = path.replace('csv', 'jpg')

    plt.rcParams['font.family'] = 'Times New Roman'# 'sans-serif'
    plt.rcParams["mathtext.fontset"] = "stix" 
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.major.width'] = 1.0
    plt.rcParams['ytick.major.width'] = 1.0
    plt.rcParams['font.size'] = 20
    plt.rcParams['axes.linewidth'] = 1.0

    fig = plt.figure(figsize=(7.0, 6.0))
    ax = fig.add_subplot(111)

    x = dat[:, 0]
    y = dat[:, 1]
    f = dat[:, 2]
    b = dat[:, 3]

    ax.plot(x, b, color='gray')

    for i in range( 4, len(dat[0, :]) ):
        peaks = dat[:, i]
        ax.plot(x, peaks, color='navy')

    ax.scatter(x, y, edgecolors='black', color='white')
    ax.plot(x, f, color='black')
    ax.yaxis.set_major_formatter(FixedOrderFormatter(3 ,useMathText=True))

    ax.set_xlim( np.min(x), np.max(x) )
    ax.invert_xaxis()

    plt.xlabel('Binding Energy [eV]')
    plt.ylabel('Intensity [counts]')

    plt.tight_layout()

    plt.savefig(fname, bbox_inches="tight", pad_inches=0.05, dpi=350)
    plt.close('all')
    # plt.show()
    
if __name__ == "__main__":
    path = '/Users/natsu/work/develop/deconv_multi_spectra_multi_peaks/code/result/opt_result_0003.csv'
    fittingPlot(path)
