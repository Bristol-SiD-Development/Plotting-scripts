#! /bin/env python
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

def degrees(radians):
    return radians * (180./np.pi)

def cleanList(l):

    thetas, good_tracks, fake_tracks, total_tracks = [], [], [], []
    for item in l:
        strings = item.split()
        theta, good_track, fake_track, total_track = float(strings[0]),  float(strings[1]),  float(strings[2]),  float(strings[3])

        if total_track != 0.:
            thetas.append(theta)
            good_tracks.append(good_track)
            fake_tracks.append(fake_track)
            total_tracks.append(total_track)

    return thetas, good_tracks, fake_tracks, total_tracks


parser = argparse.ArgumentParser(description='Plot Stuff');

parser.add_argument('-i',help='input file(s) (defaults to stdin)',type=argparse.FileType('r'),default=[sys.stdin],nargs='*')

args = parser.parse_args()

matplotlib.rcParams['text.latex.unicode']=True #for greek letters
matplotlib.rcParams['text.usetex']=True #looks better
matplotlib.rcParams['axes.color_cycle'] = ['blue', 'green', 'red', 'orange']
rc('font',**{'family':'serif','serif':['Palatino']})


integrated_efficiency_fig = plt.figure(tight_layout=True)
integrated_efficiency_ax = integrated_efficiency_fig.add_subplot(111)
plt.xlabel(r"$\theta$", size=20)
plt.ylabel("Integrated Efficiency", size=20)
#plt.xlim(0, np.pi)
#plt.vlines(np.pi/2., 0, 1)
plt.tick_params(labelsize=20)
plt.title("Changed definitions", size=20)
plt.grid(True)

integrated_efficiency_dif_fig = plt.figure(tight_layout=True)
integrated_efficiency_dif_ax = integrated_efficiency_dif_fig.add_subplot(111)
plt.xlabel(r"$\theta$", size=20)
plt.ylabel("Integrated Efficiency difference", size=20)
plt.tick_params(labelsize=20)
plt.title("Changed definitions", size=20)
plt.grid(True)


efficiency_fig = plt.figure(tight_layout=True)
efficiency_ax = efficiency_fig.add_subplot(111)
plt.xlabel(r"$\theta$", size=20)
plt.ylabel("Efficiency", size=20)
#plt.xlim(0, np.pi)
#plt.vlines(np.pi/2., 0, 1)
plt.tick_params(labelsize=20)
plt.title("Changed definitions", size=20)
plt.grid(True)


fake_rate_fig = plt.figure(tight_layout=True)
fake_rate_ax = fake_rate_fig.add_subplot(111)
plt.xlabel(r"$\theta$", size=20)
plt.ylabel("Fake Rate", size=20)
#plt.xlim(0, np.pi)
#plt.vlines(np.pi/2., 0, 0.2)
plt.tick_params(labelsize=20)
plt.title("Changed definitions", size=20)
plt.grid(True)

#thetass, good_tracks_mcs, good_tracks_recos, total_trackss, charged_mcss = [], [], [], [], []

integrated_efficiencies = {}

for inFile in args.i:	

    thetas, good_tracks, fake_tracks, total_tracks = cleanList(inFile.readlines()) #return a list of ~nan floats

    thetas = [degrees(t) for t in thetas]
    
    name = inFile.name.split('_')[0]
    integrated_efficiencies[name] = [0]

    for theta, good, false, total in zip(thetas[1:], good_tracks[1:], fake_tracks[1:], total_tracks[1:]):
        integrated_efficiencies[name].append(integrated_efficiencies[name][-1] + good/(float(total) *100.)) #we've got 100 theta bins


    integrated_efficiency_ax.plot(thetas, integrated_efficiencies[name], "-", label=name, lw=2)
    efficiency_ax.plot(thetas, [g/float(t) for  g, t in zip(good_tracks, total_tracks)], "-", label=name, lw=2)
    fake_rate_ax.plot(thetas,  [f/float(t) for f, t in zip(fake_tracks, total_tracks)], "-", label=name, lw=2)
    plt.legend()

integrated_efficiency_dif_ax.plot(thetas, [m -s for m, s in zip(integrated_efficiencies['modified'], integrated_efficiencies['sidloi3'])], label="modified - sidloi3",lw=2 )

efficiency_ax.legend()
fake_rate_ax.legend()
integrated_efficiency_dif_ax.legend(loc=2)
integrated_efficiency_ax.legend()

#integrated_efficiency_fig.savefig("report/track_int_eff_theta.pdf", bbox_inches="tight",orientation='landscape')
#integrated_efficiency_dif_fig.savefig("report/track_int_eff_dif_theta.pdf", bbox_inches="tight",orientation='landscape')
#efficiency_fig.savefig("report/track_eff_theta.pdf", bbox_inches="tight",orientation='landscape')
#fake_rate_fig.savefig("report/track_fake_rate_theta.pdf", bbox_inches="tight",orientation='landscape')

plt.show()
