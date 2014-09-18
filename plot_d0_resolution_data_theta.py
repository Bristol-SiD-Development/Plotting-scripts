#! /usr/bin/python2
# -*- coding: utf-8 -*-

import fileinput as fi
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
import sys
import argparse
from matplotlib.pyplot import cm
from mpl_toolkits.mplot3d import Axes3D

import matplotlib.colors as col

import matplotlib.mlab as mlab


from matplotlib import rc

from scipy.optimize import curve_fit

import pickle


def gauss_disc(x, *p):
    A, mu, sigma = p
    return np.floor(A*np.exp(-(x-mu)**2/(2.*sigma**2)))
def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))


parser = argparse.ArgumentParser(description='Plot Stuff');

parser.add_argument('-i',help='input file(s) (defaults to stdin)',type=argparse.FileType('rb'),default=[sys.stdin],nargs='*')

args = parser.parse_args()


matplotlib.rcParams['text.latex.unicode']=True #for greek letters
matplotlib.rcParams['text.usetex']=True #looks better
matplotlib.rcParams['axes.color_cycle'] = ['blue', 'green', 'red', 'orange']
rc('font',**{'family':'serif','serif':['Palatino']})

main_fig = plt.figure()
main_ax = main_fig.add_subplot(111)
plt.xlabel(r"$\theta$")
plt.ylabel(r"$\sigma_{d0}$")
plt.tick_params()
#plt.title(r"", size=60)
plt.grid(True)

num_histograms_to_plot = 5

def degrees(radians):
    return radians * 180. / math.pi

for inFile in args.i:	
    theta_histogram_bins_list = pickle.Unpickler(inFile).load()
    fittedThetas, fittedVariances = [], []
    figcount = 0
    for i, t in enumerate(theta_histogram_bins_list):
        
        theta, hist, bins = t
        width = (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2.

        mean = np.sum(hist * center)/float(np.sum(hist))
        sigma_est = np.sqrt(np.sum(hist * (center - mean)**2))
        
        a_est = 60#np.trapz(hist, center)

        p0 = [a_est, mean, 0.05]

        coeff, var_matrix = curve_fit(gauss, center, hist, p0=p0,maxfev=1000000,xtol=0.01)
        
        sigma_cont = math.fabs(coeff[2])    
    
        #if sigma < 100:
        if i % int(len(theta_histogram_bins_list)/num_histograms_to_plot) == 0: 
            f = plt.figure()
            plt.title(inFile.name.split('_')[0] + " " + r"$\theta = $ " + str(theta))
            plt.grid(True)
            #plt.xlim((-0.5, 0.5))
            plt.plot(center, hist, lw=2, label="Data")
            space = np.linspace(center[0], center[-1], 100)
            plt.plot(space, gauss_disc(space, *coeff), lw=2, label="Fit")
            plt.legend()
            f.savefig("report/d0_resolution_subhist_{0}_{1}.pdf".format(inFile.name.split('_')[0], figcount), bbox_inches="tight",orientation='portrait')
            figcount += 1
        fittedThetas.append(theta)
        fittedVariances.append(sigma_cont)
            #print( "sigma_est={0} sigma_cont={1} sigma_disc={2}".format(sigma_est, sigma_cont, sigma_disc) )
        


    main_ax.plot([degrees(t) for t in fittedThetas], fittedVariances, label=inFile.name.split('_')[0], lw=2)

main_ax.legend()

main_fig.savefig("report/d0_resolution_theta.pdf", bbox_inches="tight",orientation='landscape')

#plt.show()
