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

from matplotlib import rc

def normableHist(x,bins,normalise=False): #Because lxplus has a really old numpy so hist lacks the density option
	h, b = np.histogram(x,bins);
	if normalise:
		s = float(h.sum())
		#print s
		h = [num/s for num in h]
	return h, b

def cleanList(l):
	thetas, PDGs, bs, cs = [], [], [], []

	for item in l:
                strings = item.split()
                theta, PDG, b, c =  float(strings[0]), int(strings[1]), float(strings[2]), float(strings[3])
                
                thetas.append(theta)
                PDGs.append(PDG)
                bs.append(b)
                cs.append(c)
                        
	return  thetas, PDGs, bs, cs

def plot_2d_hist(xs, ys, bins, xlabel=None, ylabel=None, title=None):
        freqMap,xedges,yedges = np.histogram2d(xs, ys ,bins=bins, range=[[0, 1], [0, 1]])

        freqMapS = freqMap.swapaxes(0,1)

        elements = (len(xedges)-1) * (len(yedges)-1 )
        xpos, ypos = np.meshgrid(xedges[:-1], yedges[:-1])

        xpos = xpos.flatten()
        ypos = ypos.flatten()
        zpos = np.zeros(elements)

        fig=plt.figure(figsize=(5, 5), dpi=150)
        ax=fig.add_subplot(111, projection='3d')

        dx = 1/float(len(xedges)-1)
        dy = 1/float(len(yedges)-1)
        dz = freqMapS.flatten()

        ax.bar3d(xpos.flatten(), ypos.flatten(), dz*0, dx, dy, dz, color="blue")
        
        plt.tick_params(labelsize = 20)

        plt.xlabel(xlabel, size=20)
        plt.ylabel(ylabel,size=20)

        plt.title(title,size=25)


def plot_1d_hist(xs, bins, title=None, xlabel=None ):
        plt.figure()
        n, bins, patches = plt.hist(xs, bins, facecolor='g', alpha=0.75)
        plt.tick_params(size=40)
        plt.xlabel(xlabel, size=40)
        plt.ylabel("Frequency", size=40)
        #plt.bar(bins[:-1], n)
        plt.title(title, size=40)
        plt.grid(True)
        

def binomial_error(p, N):
        return np.sqrt( p*(1-p) / float(N) )

parser = argparse.ArgumentParser(description='Plot Stuff');


parser.add_argument('-i',help='input file(s) (defaults to stdin)',type=argparse.FileType('r'),default=[sys.stdin],nargs='*')

args = parser.parse_args()
#print args
matplotlib.rcParams['text.latex.unicode']=True #for greek letters
matplotlib.rcParams['text.usetex']=True #looks better
matplotlib.rcParams['axes.color_cycle'] = ['blue', 'green', 'red', 'orange']
rc('font',**{'family':'serif','serif':['Palatino']})

b_efficiency_v_purity_fig = plt.figure(tight_layout=True)
b_efficiency_v_purity_ax = b_efficiency_v_purity_fig.add_subplot(111)
plt.xlim((0.6,1))
plt.xlabel(r"Efficiency")
plt.ylabel(r"Purity")
plt.tick_params()
plt.title(r"B tagging")
plt.grid(True)

c_efficiency_v_purity_fig = plt.figure(tight_layout=True)
c_efficiency_v_purity_ax = c_efficiency_v_purity_fig.add_subplot(111)
plt.xlabel(r"Efficiency")
plt.ylabel(r"Purity")
plt.tick_params()
plt.title(r"C tagging")
plt.grid(True)

sidloi3_cut = 45#math.pi/3.
modified_cut = 45#math.pi/3.

def degrees(radians):
        return radians * (180. / math.pi)


for inFile in args.i:	
	thetas, PDGs, bs, cs = cleanList(inFile.readlines()) #return a list of ~nan floats
        thetas = [degrees(t) for t in thetas]
        labelStr = inFile.name.decode('utf-8').split("_")[0]

        thresholds = np.linspace(0., 1., 20)

        b_tags_barrel = np.zeros_like(thresholds)
        c_tags_barrel = np.zeros_like(thresholds)
        
        b_correct_tags_barrel = np.zeros_like(thresholds)
        c_correct_tags_barrel = np.zeros_like(thresholds)

        b_efficiencies_barrel = np.zeros_like(thresholds)
        c_efficiencies_barrel = np.zeros_like(thresholds)

        b_purities_barrel = np.zeros_like(thresholds)
        c_purities_barrel = np.zeros_like(thresholds)

        
        b_tags_endcap = np.zeros_like(thresholds)
        c_tags_endcap = np.zeros_like(thresholds)
        
        b_correct_tags_endcap = np.zeros_like(thresholds)
        c_correct_tags_endcap = np.zeros_like(thresholds)

        b_efficiencies_endcap = np.zeros_like(thresholds)
        c_efficiencies_endcap = np.zeros_like(thresholds)

        b_purities_endcap = np.zeros_like(thresholds)
        c_purities_endcap = np.zeros_like(thresholds)

        whichDetector = inFile.name.split('_')[0]
        cut = 0
        if whichDetector == "sidloi3":
                cut = sidloi3_cut
        elif whichDetector == "modified":
                cut = modified_cut
        
        barrelString = r" barrel ($|90 - \theta| < " + "{0}".format(cut) + "$)"
        endcapString = r" endcap ($|90 - \theta| \geq " + "{0}".format(cut) + "$)"

        for i, threshold in enumerate(thresholds):
                numBs_barrel = 0
                numCs_barrel = 0

                numBs_endcap = 0
                numCs_endcap = 0

                for theta, PDG, b, c in zip(thetas, PDGs, bs, cs):
                        if abs(90 - theta) < cut:
                                if PDG == 5:
                                        numBs_barrel += 1
                                elif PDG == 4:
                                        numCs_barrel += 1

                                if b > threshold:
                                        b_tags_barrel[i] += 1
                                        if PDG == 5:
                                                b_correct_tags_barrel[i] += 1
                                if c > threshold:
                                        c_tags_barrel[i] += 1
                                        if PDG == 4:
                                                c_correct_tags_barrel[i] += 1
                        else:
                                if PDG == 5:
                                        numBs_endcap += 1
                                elif PDG == 4:
                                        numCs_endcap += 1

                                if b > threshold:
                                        b_tags_endcap[i] += 1
                                        if PDG == 5:
                                                b_correct_tags_endcap[i] += 1
                                if c > threshold:
                                        c_tags_endcap[i] += 1
                                        if PDG == 4:
                                                c_correct_tags_endcap[i] += 1
                
                
                if numBs_barrel != 0:
                        b_efficiencies_barrel[i] = b_correct_tags_barrel[i] / float(numBs_barrel)
                else:
                        b_efficiencies_barrel[i] = 0

                if numCs_barrel != 0:
                        c_efficiencies_barrel[i] = c_correct_tags_barrel[i] / float(numCs_barrel)
                else:
                        c_efficiencies_barrel[i] = 0
                
                if b_tags_barrel[i] == 0:
                        b_purities_barrel[i] = 1
                else:
                        b_purities_barrel[i] = b_correct_tags_barrel[i] / float(b_tags_barrel[i])
                if c_tags_barrel[i] == 0:
                        c_purities_barrel[i] = 1
                else:
                        c_purities_barrel[i] = c_correct_tags_barrel[i] / float(c_tags_barrel[i])
        
        
                if numBs_endcap != 0:
                        b_efficiencies_endcap[i] = b_correct_tags_endcap[i] / float(numBs_endcap)
                else:
                        b_efficiencies_endcap[i] = 0

                if numCs_endcap != 0:
                        c_efficiencies_endcap[i] = c_correct_tags_endcap[i] / float(numCs_endcap)
                else:
                        c_efficiencies_endcap[i] = 0
                
                if b_tags_endcap[i] == 0:
                        b_purities_endcap[i] = 1
                else:
                        b_purities_endcap[i] = b_correct_tags_endcap[i] / float(b_tags_endcap[i])
                if c_tags_endcap[i] == 0:
                        c_purities_endcap[i] = 1
                else:
                        c_purities_endcap[i] = c_correct_tags_endcap[i] / float(c_tags_endcap[i])
        
        
        """
        #Plot efficiencies
        plt.figure()
        plt.tick_params(labelsize = 40)
        plt.grid()

        plt.xlabel("Threshold", size=40)
        plt.ylabel("Efficiency",size=40)

        plt.plot(thresholds, b_efficiencies_barrel, color = "blue", label='b tagging',lw=2)
        plt.plot(thresholds, c_efficiencies_barrel, color = "red", label='c tagging',lw=2)

        plt.legend(fontsize=25)
        plt.title(labelStr, size=40)

        #Plot purities
        plt.figure()
        plt.tick_params(labelsize = 40)
        plt.grid()

        plt.xlabel("Threshold", size=40)
        plt.ylabel("Purity",size=40)

        plt.plot(thresholds, b_purities_barrel, color = "blue", label='b tagging',lw=2)
        plt.plot(thresholds, c_purities_barrel, color = "red", label='c tagging',lw=2)

        plt.legend(fontsize=25)
        plt.title(labelStr, size=40)
        """

        
        #l=zip(*[[e, p] for e, p in zip(b_efficiencies_barrel, b_purities_barrel) if e > 0.1])
        
        #b_efficiencies_barrel, b_purities_barrel = np.array(l[0]), np.array(l[1])
        

        b_base_line_barrel, = b_efficiency_v_purity_ax.plot(b_efficiencies_barrel, b_purities_barrel, label=labelStr + barrelString, lw=2)
        b_eff_err_barrel = binomial_error(b_efficiencies_barrel, numBs_barrel)
        b_pure_err_barrel = binomial_error(b_purities_barrel, numBs_barrel)

        b_efficiency_v_purity_ax.plot(b_efficiencies_barrel - b_eff_err_barrel , b_purities_barrel - b_pure_err_barrel, '--', label=None, lw=1, color= b_base_line_barrel.get_color())
        b_efficiency_v_purity_ax.plot(b_efficiencies_barrel + b_eff_err_barrel , b_purities_barrel + b_pure_err_barrel, '--', label=None, lw=1, color= b_base_line_barrel.get_color())

        #l=zip(*[[e, p] for e, p in zip(b_efficiencies_endcap, b_purities_endcap) if e > 0.1])
        
        #b_efficiencies_endcap, b_purities_endcap = np.array(l[0]), np.array(l[1])
        

        b_base_line_endcap, = b_efficiency_v_purity_ax.plot(b_efficiencies_endcap, b_purities_endcap, label=labelStr + endcapString, lw=2)

        b_eff_err_endcap = binomial_error(b_efficiencies_endcap, numBs_endcap)
        b_pure_err_endcap = binomial_error(b_purities_endcap, numBs_endcap)

        b_efficiency_v_purity_ax.plot(b_efficiencies_endcap - b_eff_err_endcap , b_purities_endcap - b_pure_err_endcap, '--', label=None, lw=1, color= b_base_line_endcap.get_color())
        b_efficiency_v_purity_ax.plot(b_efficiencies_endcap + b_eff_err_endcap , b_purities_endcap + b_pure_err_endcap, '--', label=None, lw=1, color= b_base_line_endcap.get_color())

        c_base_line_barrel, = c_efficiency_v_purity_ax.plot(c_efficiencies_barrel, c_purities_barrel, label=labelStr+ barrelString, lw=2)
        c_eff_err_barrel = binomial_error(c_efficiencies_barrel, numCs_barrel)
        c_pure_err_barrel = binomial_error(c_purities_barrel, numCs_barrel)

        c_efficiency_v_purity_ax.plot(c_efficiencies_barrel - c_eff_err_barrel , c_purities_barrel - c_pure_err_barrel, '--', label=None, lw=1, color= b_base_line_barrel.get_color())
        c_efficiency_v_purity_ax.plot(c_efficiencies_barrel + c_eff_err_barrel , c_purities_barrel + c_pure_err_barrel, '--', label=None, lw=1, color= b_base_line_barrel.get_color())

        c_base_line_endcap, = c_efficiency_v_purity_ax.plot(c_efficiencies_endcap, c_purities_endcap, label=labelStr+ endcapString, lw=2)


        c_eff_err_endcap = binomial_error(c_efficiencies_endcap, numCs_endcap)
        c_pure_err_endcap = binomial_error(c_purities_endcap, numCs_endcap)

        c_efficiency_v_purity_ax.plot(c_efficiencies_endcap - c_eff_err_endcap , c_purities_endcap - c_pure_err_endcap, '--', label=None, lw=1, color= b_base_line_endcap.get_color())
        c_efficiency_v_purity_ax.plot(c_efficiencies_endcap + c_eff_err_endcap , c_purities_endcap + c_pure_err_endcap, '--', label=None, lw=1, color= b_base_line_endcap.get_color())

b_efficiency_v_purity_ax.legend(fontsize=12, loc=3)
c_efficiency_v_purity_ax.legend(fontsize=12, loc=3)
b_efficiency_v_purity_fig.savefig("report/b_pure_eff_threshold_with_theta_cuts.pdf", bbox_inches="tight",orientation='landscape')
c_efficiency_v_purity_fig.savefig("report/c_pure_eff_threshold_with_theta_cuts.pdf", bbox_inches="tight",orientation='landscape')
#plt.show()


