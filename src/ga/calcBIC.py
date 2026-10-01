import os, sys
from fitting_plot import fittingPlot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.special as sc
from glob import glob
from matplotlib.ticker import ScalarFormatter

class FixedOrderFormatter(ScalarFormatter):
    def __init__(self, order_of_mag=0, useOffset=True, useMathText=True):
        self._order_of_mag = order_of_mag
        ScalarFormatter.__init__(self, useOffset=useOffset, useMathText=useMathText)
    def _set_orderOfMagnitude(self, range):
        self.orderOfMagnitude = self._order_of_mag

#BASE_PATH = '/Users/natsu/work/develop/deconv_multi_single_multi_peaks/code'
#BASE_PATH = '/home/shinotsuka/Simple-Peak-Fitting/code'
BASE_PATH = './'

def calcLogLikelihood(y, f):
    log_likelihood = np.sum( - f + y * np.log( f ) - sc.gammaln( y + 1 ) )
    return log_likelihood

def calcPenaltyTerm(N, K, num_param=4):
    num_parameter = ( K * num_param + 2 )
    penalty = ( num_parameter * np.log(N) )
    return penalty

def calcBIC(output, K, path):
    
    df = pd.read_csv(path)
    fittingPlot(path)

    N = len( df.index )

    yy = np.array( df['spectrum'] )
    ff = np.array( df['fitting'] )

    log_likelihood = calcLogLikelihood(yy, ff)
    penalty_term = calcPenaltyTerm(N, K)

    BIC = -2 * log_likelihood + penalty_term

    print( 'k = {} : BIC = {}'.format(K, BIC) )

    output = output.append(
        {
            'K': K,
            'log_likelihood': ( -2 * log_likelihood ),
            'penalty' : penalty_term,
            'BIC': BIC
        },
        ignore_index=True
    )

    return output


def selectNumPeaks(K_min, K_max):
    output = pd.DataFrame(columns=['K', 'log_likelihood', 'penalty', 'BIC'])
    for k in range(K_min, K_max):
        k_path = os.path.join( BASE_PATH, 'num_peaks_{:02}/result/opt_result.csv'.format(k) )
        if not os.path.exists(k_path):
            continue
        output = calcBIC(output, k, k_path)
        plotLogEnergy(
            os.path.join( BASE_PATH, 'num_peaks_{:02}/result/energy_log.csv'.format(k) )
        )

    print(output)

    output.to_csv('output.csv')
    plotBIC(output)

def plotBIC(output):
    plt.rcParams['font.family'] = 'Times New Roman'# 'sans-serif'
    plt.rcParams["mathtext.fontset"] = "stix" 
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.major.width'] = 1.0
    plt.rcParams['ytick.major.width'] = 1.0
    plt.rcParams['font.size'] = 18
    plt.rcParams['axes.linewidth'] = 1.0

    fig = plt.figure(figsize=(8.0, 6.0))
    ax = fig.add_subplot(111)

    ax.scatter(output['K'], output['BIC'], edgecolors='k', color='w', s=150, lw=2)
    ax.yaxis.set_major_formatter(FixedOrderFormatter(3 ,useMathText=True))

    ax.set_xlabel('Number of peaks $K$')
    ax.set_ylabel('Bayesian Information Criterion; BIC')

    ax.set_xticks(output['K'])

    plt.tight_layout()
    plt.savefig('BIC_plot.jpg', bbox_inches="tight", pad_inches=0.05, dpi=350)
    # plt.show()

def plotLogEnergy(fpath):
    log_e = pd.read_csv( fpath )
    for c in log_e.columns[1:]:
        plt.plot( log_e['id'], np.log10(log_e[c]) )

    plt.xlabel('loop')
    plt.ylabel('log10(E)')
    plt.xscale('log')

    fname = fpath.replace('csv', 'jpg')
    plt.savefig(fname, bbox_inches="tight", pad_inches=0.05, dpi=350)
    plt.close('all')

if __name__ == "__main__":
    K_min = 1
    K_max = 5
    selectNumPeaks(K_min, K_max+1)
    
