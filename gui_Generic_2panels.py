#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 16:54:58 2018

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
VERYSMALL_SIZE = SMALL_SIZE - 5
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
def gui_generic_2panelplots(irun, iz, iZZ, idens, PlotType, PlotDict, idisplay):
    myhdf5file = '%s/%s.hdf5'%(filebase, runname[irun])
    with h5py.File(myhdf5file, "r") as f:
        RedshiftBins       = f['TableBins/RedshiftBins'][:]
        MetallicityBins    = f['TableBins/MetallicityBins'][:]
        TemperatureBins    = f['TableBins/TemperatureBins'][:]
        DensityBins        = f['TableBins/DensityBins'][:]
   
        if PlotType == 1:   # Thermal Equilibrium
            Q1   = f['ThermEq/'+PlotDict['dset1']][:]
            Q2   = f['ThermEq/'+PlotDict['dset2']][:]
        else:
            Q1 = f['Tdep/'+PlotDict['dset1']][:]
            Q2 = f['Tdep/'+PlotDict['dset2']][:]
        if PlotType == 2:
            Teq = f['ThermEq/Temperature'][:]
        
        if PlotDict['shortname'] == 'CO': 
            abundances = f['TotalAbundances'][:]
            lognCnH = abundances[iZZ,2]
            Q1 = Q1 - lognCnH
            Q2 = Q2 - lognCnH
       
    titlestring, outputfile = getfilenames(runname[irun], RedshiftBins[iz], DensityBins[idens],\
                                             MetallicityBins[iZZ], iz, iZZ, idens, PlotType, PlotDict['shortname'])
    
    fig = plt.figure()
    fig.suptitle(titlestring, fontsize = 12)
    
    if PlotType == 2:  # 2D plot
        fig.subplots_adjust(left = 0.12, right = 0.97, bottom = 0.2, top = 0.8)
        gs = gridspec.GridSpec(1,2, wspace = 0.3)
        ax = plt.subplot(gs[0]) 
        ax.set_title(PlotDict['top2D1'])
        ax.set_xlabel('log n [cm$^{-3}$]')
        ax.set_ylabel('log T [K]')
        ax.xaxis.set_ticks(np.arange(DensityBins[0], DensityBins[-1] + 2., 2.))
        ax.yaxis.set_ticks(np.arange(2., 9., 2.))
        extent = (DensityBins[0], DensityBins[-1], TemperatureBins[0], TemperatureBins[-1])  
        if PlotDict['shortname'] == 'CO':
            im = ax.imshow(Q1[iz, :, iZZ, :], interpolation='none', origin='lower', \
                 extent = extent, aspect = 'auto',vmin = -20, vmax = 0.)  
        else:
            im = ax.imshow(Q1[iz, :, iZZ, :], interpolation='none', origin='lower', \
                 extent = extent, aspect = 'auto')  
        ax.autoscale(False)
        ax.plot(DensityBins, Teq[iz, iZZ, :], color = 'white', linestyle = '--')    
        cb = plt.colorbar(im,ax=ax, orientation='horizontal', pad = 0.25)
        cb.set_label(PlotDict['label1'])         
        cb.set_ticks(np.arange(-20., 5., 5.))        
       
        ax = plt.subplot(gs[1]) 
        ax.set_title(PlotDict['top2D2'])
        ax.set_xlabel('log n [cm$^{-3}$]')
        ax.set_ylabel('log T [K]')
        ax.xaxis.set_ticks(np.arange(DensityBins[0], DensityBins[-1] + 2., 2.))
        ax.yaxis.set_ticks(np.arange(2., 9., 2.))
        extent = (DensityBins[0], DensityBins[-1], TemperatureBins[0], TemperatureBins[-1]) 
        if PlotDict['shortname'] == 'CO':
            im = ax.imshow(Q2[iz, :, iZZ, :], interpolation='none', origin='lower', \
                           extent = extent, aspect = 'auto', vmin = -20., vmax = 0.) 
        else:
            im = ax.imshow(Q2[iz, :, iZZ, :], interpolation='none', origin='lower', \
                           extent = extent, aspect = 'auto')  
        ax.autoscale(False)
        ax.plot(DensityBins, Teq[iz, iZZ, :], color = 'white', linestyle = '--')    
        cb = plt.colorbar(im,ax=ax, orientation='horizontal', pad = 0.25)
        cb.set_label(PlotDict['label2'])        
        cb.set_ticks(np.arange(-20., 5., 5.))        
    else: 
        fig.subplots_adjust(left = 0.2, right = 0.95, bottom = 0.2, top = 0.8)
        gs = gridspec.GridSpec(1,1)   
        ax = plt.subplot(gs[0]) 
        ax.set_title(PlotDict['name'])
        if PlotDict['shortname'] == 'CO':
            ax.set_ylim(-20.,0.5)
        if PlotType == 0: # constant density
            ax.set_xlabel('log T [K]')
            ax.set_ylabel(PlotDict['label'])
            ax.plot(TemperatureBins, Q1[iz,:,iZZ,idens], label = PlotDict['leglabel1'])
            ax.plot(TemperatureBins, Q2[iz,:,iZZ,idens], label = PlotDict['leglabel2'])
        elif PlotType == 1: # equilibrium temperature
            ax.set_xlabel('log n [cm$^{-3}$]')
            ax.xaxis.set_ticks(np.arange(DensityBins[0], DensityBins[-1] + 2., 2.))
            ax.set_ylabel(PlotDict['label'])
            ax.plot(DensityBins, Q1[iz,iZZ,:], label = PlotDict['leglabel1'])
            ax.plot(DensityBins, Q2[iz,iZZ,:], label = PlotDict['leglabel2'])
   
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc=3, fontsize = 12)
    
    if idisplay == 1:
        outputfile = "tmp.png"
    
    fig.savefig(outputfile, dpi = 100)
    print('Plot saved as: %s'%(outputfile))
    plt.close('all')

    if idisplay == 1:
        os.system("display %s &"%(outputfile))  
    
    
    
    
    
    
