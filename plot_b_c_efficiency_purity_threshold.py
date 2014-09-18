#! /usr/bin/python2
# -*- coding: utf-8 -*-
import fileinput as fi
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math as maths #localisation
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
print args

matplotlib.rcParams['text.latex.unicode']=True #for greek letters
matplotlib.rcParams['text.usetex']=True #looks better
matplotlib.rcParams['axes.color_cycle'] = ['blue', 'orange', 'green', 'red']
rc('font',**{'family':'serif','serif':['Palatino']})


b_efficiency_v_purity_fig = plt.figure()
b_efficiency_v_purity_ax = b_efficiency_v_purity_fig.add_subplot(111)
plt.xlabel("Efficiency")
plt.ylabel("Purity")
plt.xlim((0.6,1))
#plt.tick_params(labelsize=60)
plt.title("B tagging")
plt.grid(True)

c_efficiency_v_purity_fig = plt.figure()
c_efficiency_v_purity_ax = c_efficiency_v_purity_fig.add_subplot(111)
plt.xlabel("Efficiency")
plt.ylabel("Purity")
#plt.tick_params(labelsize=40)
plt.title("C tagging")
plt.grid(True)


for inFile in args.i:	
	thetas, PDGs, bs, cs = cleanList(inFile.readlines()) #return a list of ~nan floats

        numEvents = len(PDGs)
        numBs = len(filter(lambda x: x == 5, PDGs))
        numCs = len(filter(lambda x: x == 4, PDGs))
        numUDS = len(filter(lambda x: (x == 1) or (x == 2) or (x == 3), PDGs))
        numTaus = len(filter(lambda x: x == 15, PDGs))


        print("From the MC there were the following particles (unfiltered):\nUDS: {}\nB: {}\nC: {}\nTau: {}\nOther: {}\nTotal: {}\n".format(numUDS,
                                                                                                                                            numBs, 
                                                                                                                                            numCs, 
                                                                                                                                            numTaus, 
                                                                                                                                            numEvents - numUDS - numBs-numCs,
                                                                                                                                            numEvents))

        
        

        #Remove taus!
        #PDGs, bs, cs = zip(*[(p, b, c) for (p, b, c) in zip(PDGs[:], bs[:], cs[:]) if (p != 15)])
        #Remove everything else
        PDGs, bs, cs = zip(*[(p, b, c) for (p, b, c) in zip(PDGs[:], bs[:], cs[:]) if ((p == 1) or(p == 2) or(p == 3) or(p == 4) or(p == 5) or(p == 6))])
        
        #PDGs, bs, cs = zip(*[(h, p, b, c) for (h, p, b, c) in zip(hit_locs[:], PDGs[:], bs[:], cs[:]) if h == loc])
        
        numEvents = len(PDGs)
        numBs = len(filter(lambda x: x == 5, PDGs))
        numCs = len(filter(lambda x: x == 4, PDGs))
        numUDS = len(filter(lambda x: (x == 1) or (x == 2) or (x == 3), PDGs))
        numTaus = len(filter(lambda x: x == 15, PDGs))

        
        print("From the MC there were the following particles (filtered):\nUDS: {}\nB: {}\nC: {}\nTau: {}\nOther: {}\nTotal: {}\n".format(numUDS,
                                                                                                                                          numBs, 
                                                                                                                                          numCs, 
                                                                                                                                          numTaus, 
                                                                                                                                          numEvents - numUDS - numBs-numCs,
                                                                                                                                          numEvents))

        
        labelStr = inFile.name.decode('utf-8').split('_')[0] 

        thresholds = np.linspace(0., 1., 100)

        b_tags = np.zeros_like(thresholds)
        c_tags = np.zeros_like(thresholds)
        
        b_correct_tags = np.zeros_like(thresholds)
        c_correct_tags = np.zeros_like(thresholds)

        b_efficiencies = np.zeros_like(thresholds)
        c_efficiencies = np.zeros_like(thresholds)

        b_purities = np.zeros_like(thresholds)
        c_purities = np.zeros_like(thresholds)

        for i, threshold in enumerate(thresholds):
                for PDG, b, c in zip(PDGs, bs, cs):
                        if b > threshold :
                                b_tags[i] += 1
                                if PDG == 5:
                                        b_correct_tags[i] += 1
                        if c > threshold:
                                c_tags[i] += 1
                                if PDG == 4:
                                        c_correct_tags[i] += 1
                

                if numBs != 0:
                        b_efficiencies[i] = b_correct_tags[i] / float(numBs)
                else:
                        b_efficiencies[i] = 0

                if numCs != 0:
                        c_efficiencies[i] = c_correct_tags[i] / float(numCs)
                else:
                        c_efficiencies[i] = 0
                
                if b_tags[i] == 0:
                        b_purities[i] = 1
                else:
                        b_purities[i] = b_correct_tags[i] / float(b_tags[i])
                if c_tags[i] == 0:
                        c_purities[i] = 1
                else:
                        c_purities[i] = c_correct_tags[i] / float(c_tags[i])
        
        """
        #Plot efficiencies
        plt.figure()
        plt.tick_params(labelsize = 40)
        plt.grid()

        plt.xlabel("Threshold", size=40)
        plt.ylabel("Efficiency",size=40)

        plt.plot(thresholds, b_efficiencies, color = "blue", label='b tagging',lw=2)
        plt.plot(thresholds, c_efficiencies, color = "red", label='c tagging',lw=2)

        plt.legend(fontsize=40)
        plt.title(labelStr, size=40)

        #Plot purities
        plt.figure()
        plt.tick_params(labelsize = 40)
        plt.grid()

        plt.xlabel("Threshold", size=40)
        plt.ylabel("Purity",size=40)

        plt.plot(thresholds, b_purities, color = "blue", label='b tagging',lw=2)
        plt.plot(thresholds, c_purities, color = "red", label='c tagging',lw=2)

        plt.legend(fontsize=40)
        plt.title(labelStr, size=40)
        """
        
        b_base_line, = b_efficiency_v_purity_ax.plot(b_efficiencies, b_purities, label=labelStr, lw=1)

        b_eff_err = binomial_error(b_efficiencies, numBs)
        b_pure_err = binomial_error(b_purities, numBs)

        b_efficiency_v_purity_ax.plot(b_efficiencies - b_eff_err , b_purities - b_pure_err, '--', label=None, lw=1, color= b_base_line.get_color())
        b_efficiency_v_purity_ax.plot(b_efficiencies + b_eff_err , b_purities + b_pure_err, '--', label=None, lw=1, color= b_base_line.get_color())


        c_base_line, = c_efficiency_v_purity_ax.plot(c_efficiencies, c_purities, label=labelStr, lw=1)

        c_eff_err = binomial_error(c_efficiencies, numCs)
        c_pure_err = binomial_error(c_purities, numCs)

        c_efficiency_v_purity_ax.plot(c_efficiencies - c_eff_err , c_purities - c_pure_err, '--', label=None, lw=1, color= b_base_line.get_color())
        c_efficiency_v_purity_ax.plot(c_efficiencies + c_eff_err , c_purities + c_pure_err, '--', label=None, lw=1, color= b_base_line.get_color())

        


        """
        #Plot Efficiency v purity
        plt.figure()
        plt.tick_params(labelsize = 40)
        plt.grid()

        plt.xlabel("Efficiency", size=40)
        plt.ylabel("Purity",size=40)

        plt.plot(b_efficiencies, b_purities, color = "blue", label='b tagging',lw=2)
        plt.plot(c_efficiencies, c_purities, color = "red", label='c tagging',lw=2)

        plt.legend(fontsize=40)
        plt.title(labelStr, size=40)
        """
        """
        #plt.rc('text', usetex=True)
        #plt.rc('font', family='serif')

        #Plot 3d b / c tag histogram
        plot_2d_hist(bs, cs, bins=10, xlabel="B likeness", ylabel="C likeness", title=labelStr)

        #B_likeness and C_likeness of things split on their true flavour
        b_b_likeneses = [b for (p, b) in zip(PDGs, bs) if p == 5]
        b_c_likeneses = [c for (p, c) in zip(PDGs, cs) if p == 5]
        plot_2d_hist(b_b_likeneses, b_c_likeneses, bins=7, xlabel="B likeness", ylabel="C likeness", title="B decays"+labelStr)
        #plot_1d_hist(b_b_likeneses, bins=49, xlabel = "B likeness", title=r"H\rightarrow BB decays")
        #plot_1d_hist(b_c_likeneses, bins=49, xlabel = "C likeness", title=r"H\rightarrow BB decays")

        c_b_likeneses = [b for (p, b) in zip(PDGs, bs) if p == 4]
        c_c_likeneses = [c for (p, c) in zip(PDGs, cs) if p == 4]
        plot_2d_hist(c_b_likeneses, c_c_likeneses, bins=7, xlabel="B likeness", ylabel="C likeness", title="C decays"+labelStr)
        #plot_1d_hist(c_b_likeneses, bins=49, xlabel = "B likeness", title=r"H\rightarrow CC decays")
        #plot_1d_hist(c_c_likeneses, bins=49, xlabel = "C likeness", title=r"H\rightarrow CC decays")

        uds_b_likeneses = [b for (p, b) in zip(PDGs[:], bs[:]) if ((p == 1) or (p == 2) or (p == 3))]
        uds_c_likeneses = [c for (p, c) in zip(PDGs[:], cs[:]) if ((p == 1) or (p == 2) or (p == 3))]
        plot_2d_hist(uds_b_likeneses, uds_c_likeneses, bins=7, xlabel="B likeness", ylabel="C likeness", title="UDS decays"+labelStr)
        #plot_1d_hist(uds_b_likeneses, bins=49, xlabel = "B likeness", title=r"H\rightarrow UDS decays")
        #plot_1d_hist(uds_c_likeneses, bins=49, xlabel = "C likeness", title=r"H\rightarrow UDS decays")
        

        tau_b_likeneses = [b for (p, b) in zip(PDGs[:], bs[:]) if p == 15]
        tau_c_likeneses = [c for (p, c) in zip(PDGs[:], cs[:]) if p == 15]
        plot_2d_hist(tau_b_likeneses, tau_c_likeneses, bins=7, xlabel="B likeness", ylabel="C likeness", title="Tau decays")
        plot_1d_hist(tau_b_likeneses, bins=49, xlabel = "B likeness", title=r'H\rightarrow\tau\tau decays')
        plot_1d_hist(tau_c_likeneses, bins=49, xlabel = "C likeness", title=r"H\rightarrow\tau\tau decays")
        """
b_efficiency_v_purity_ax.legend(fontsize=12, loc=3)
c_efficiency_v_purity_ax.legend(fontsize=12, loc=3)

#b_efficiency_v_purity_fig.savefig("report/b_pure_eff_threshold.pdf", bbox_inches="tight",orientation='landscape')
#c_efficiency_v_purity_fig.savefig("report/c_pure_eff_threshold.pdf", bbox_inches="tight",orientation='landscape')

plt.show()
