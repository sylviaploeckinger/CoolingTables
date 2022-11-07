#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 16:24:13 2018

@author: sylviaploeckinger
"""

#####################################
# Import standard python libraries
#####################################
import sys
import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

#####################################
# Import local python files 
#####################################
from global_parameters import runname
from get_strings_for_plots import getfilenames

#####################################
# Path to the cooling table
#####################################
filebase = sys.argv[1]

#####################################
# Some plot properties 
#####################################
SMALL_SIZE = 12
VERYSMALL_SIZE = SMALL_SIZE - 3
MEDIUM_SIZE = SMALL_SIZE + 2
BIG_SIZE    = SMALL_SIZE + 4
BIGGER_SIZE = SMALL_SIZE + 6
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=VERYSMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIG_SIZE)  # fontsize of the figure title

cmap = plt.cm.get_cmap('viridis')

#####################################
# Plotting routine 
#####################################
def gui_cooling_plots(irun,iz,iZZ,idens,iplottype, idisplay):
    if iplottype > 1:
        return 
    
    ymin = -35.
    ymax = -18.

    ########################################
    # Split the cooling contributions into:
    # panel 1: primoridal elements
    # panel 2: atomic metal cooling
    # panel 3: remaining cooling channels
    # This split is hardcoded here
    ########################################
    primcool = [0,1,12,14,15,17,18,20]
    atomcool = [2,3,4,5,6,7,8,9,10,21]
    restcool = [11,13,16,19,21]    
    
    myhdf5file = '%s/%s.hdf5'%(filebase, runname[irun])
    with h5py.File(myhdf5file, "r") as f:
        RedshiftBins       = f['TableBins/RedshiftBins'][:]
        MetallicityBins    = f['TableBins/MetallicityBins'][:]
        TemperatureBins    = f['TableBins/TemperatureBins'][:]
        DensityBins        = f['TableBins/DensityBins'][:]
        IdentifierCooling = f['IdentifierCooling'][:]
        if iplottype == 0:
            Cooling = f['Tdep/Cooling'][:]
            Cool1D  = Cooling[iz, :, iZZ, idens, :]         
            xx      = TemperatureBins
            xlab    = 'log T [K]'
        else:
            Cooling = f['ThermEq/Cooling'][:]
            Cool1D  = Cooling[iz, iZZ, :, :]
            xx      = DensityBins
            xlab    = 'log n [cm$^{-3}$]'

    xmin = xx[0]
    xmax = xx[-1]

    titlestring, outputfile = getfilenames(runname[irun], RedshiftBins[iz], DensityBins[idens],\
                                             MetallicityBins[iZZ], iz, iZZ, idens, iplottype, 'Cool')
    
    fig = plt.figure()
    fig.set_size_inches(10,6.2,forward=True)
    fig.suptitle(titlestring)
    fig.subplots_adjust(left = 0.1, right = 0.95, bottom = 0.15, top = 0.7)
    gs = gridspec.GridSpec(1,3,wspace=0, hspace=0)

    ax = plt.subplot(gs[0])
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    ax.set_ylabel('log $\Lambda_{\mathrm{cool}}$/$n_{\mathrm{H}}^2$ [erg cm$^{3}$ s$^{-1}$]')

    ax.set_xlabel(xlab)
    ax.xaxis.set_ticks(np.arange(xx[0], xx[-2], 1.))
    ax.plot(xx, np.log10(np.power(10., Cool1D[..., -2]) + np.power(10., Cool1D[..., -1])), color = 'grey', lw = 2, label = 'Total')

    icount = 0
    c = cmap(np.linspace(0,1,len(primcool)-1))
    for i in range (len(primcool)):
        icool = primcool[i]
        if icool == primcool[-1]:
            linec = 'black'
            lines = ':'
        else:
            linec = c[icount]
            if icount%3 == 0:
                lines = '-'
            elif icount%3 == 1:
                lines = '--'
            else: 
                lines = '-.' 
                
        ax.plot(xx, Cool1D[..., icool], color = linec, lw = 2, ls = lines, label = '%s'%(IdentifierCooling[icool].decode('utf-8')))
        icount = icount + 1


    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, bbox_to_anchor=(0., 1.02, 1., .102), loc=3, mode = 'expand', borderaxespad=0., ncol = 2, fontsize = VERYSMALL_SIZE, handlelength = 4)

    ax = plt.subplot(gs[1])
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    ax.tick_params(labelleft=False)

    ax.set_xlabel(xlab)
    ax.xaxis.set_ticks(np.arange(xx[0], xx[-2], 1.))
    ax.plot(xx, np.log10(np.power(10., Cool1D[..., -2]) + np.power(10., Cool1D[..., -1])), color = 'grey', lw = 2, label = 'Total')
    icount = 0
    c = cmap(np.linspace(0,1,len(atomcool)-1))
    for i in range (len(atomcool)):
        icool = atomcool[i]
        if icool == atomcool[-1]:
            linec = 'black'
            lines = ':'
        else:
            linec = c[icount]
            if icount%3 == 0:
                lines = '-'
            elif icount%3 == 1:
                lines = '--'
            else: 
                lines = '-.'
                
        ax.plot(xx, Cool1D[..., icool], color = linec, lw = 2, ls = lines, label = '%s'%(IdentifierCooling[icool].decode('utf-8')))
        icount = icount + 1

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, bbox_to_anchor=(0., 1.02, 1., .102), loc=3, mode = 'expand', borderaxespad=0., ncol = 2, fontsize = VERYSMALL_SIZE, handlelength = 4)


    ax = plt.subplot(gs[2])
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    ax.tick_params(labelleft=False)

    ax.set_xlabel(xlab)
    ax.xaxis.set_ticks(np.arange(xx[0], xx[-1], 1.))
    ax.plot(xx, np.log10(np.power(10., Cool1D[..., -2]) + np.power(10., Cool1D[..., -1])), color = 'grey', lw = 2, label = 'Total')
    icount = 0
    c = cmap(np.linspace(0,1,len(restcool)-1))
    for i in range (len(restcool)):
        icool = restcool[i]
        if icool == restcool[-1]:
            linec = 'black'
            lines = ':'
        else:
            linec = c[icount]
            if icount%3 == 0:
                lines = '-'
            elif icount%3 == 1:
                lines = '--'
            else: 
                lines = '-.' 
                
        ax.plot(xx, Cool1D[..., icool], color = linec, lw = 2, ls = lines, label = '%s'%(IdentifierCooling[icool].decode('utf-8')))
        icount = icount + 1


    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, bbox_to_anchor=(0., 1.02, 1., .102), loc=3, mode = 'expand', borderaxespad=0., ncol = 2, fontsize = VERYSMALL_SIZE, handlelength = 4)

    if idisplay == 1:
        outputfile = "tmp.png"
    
    fig.savefig(outputfile, dpi = 100)
    print('Plot saved as: %s'%(outputfile))
    plt.close('all')

    if idisplay == 1:
        os.system("display %s &"%(outputfile))  

