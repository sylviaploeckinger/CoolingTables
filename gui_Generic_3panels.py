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
def gui_generic_3panelplots(irun, iz, iZZ, idens, PlotType, PlotDict, idisplay):
    myhdf5file = '%s/%s.hdf5'%(filebase, runname[irun])
    with h5py.File(myhdf5file, "r") as f:
        RedshiftBins       = f['TableBins/RedshiftBins'].value
        MetallicityBins    = f['TableBins/MetallicityBins'].value
        TemperatureBins    = f['TableBins/TemperatureBins'].value
        DensityBins        = f['TableBins/DensityBins'].value
   
        if PlotType == 1:   # Thermal Equilibrium
            Q   = f['ThermEq/'+PlotDict['dset']].value
        else:
            Q = f['Tdep/'+PlotDict['dset']].value
        if PlotType == 2:
            Teq = f['ThermEq/Temperature'].value
        
    Q1 = Q[...,0]
    Q2 = Q[...,1]
    Q3 = Q[...,2]
    
    titlestring, outputfile = getfilenames(runname[irun], RedshiftBins[iz], DensityBins[idens],\
                                             MetallicityBins[iZZ], iz, iZZ, idens, PlotType, PlotDict['shortname'])
    
    fig = plt.figure()
    fig.suptitle(titlestring, fontsize = 12)
    
    if PlotType == 2:  # 2D plot
        fig.set_size_inches(7.5,4,forward=True)
        fig.subplots_adjust(left = 0.1, right = 0.97, bottom = 0.15, top = 0.8)
        gs = gridspec.GridSpec(1,3, wspace = 0.3)
        ax = plt.subplot(gs[0]) 
        ax.set_title(PlotDict['top2D1'])
        ax.set_xlabel('log n [cm$^{-3}$]')
        ax.set_ylabel('log T [K]')
        ax.xaxis.set_ticks(np.arange(DensityBins[0], DensityBins[-1]+2., 2.))
        ax.yaxis.set_ticks(np.arange(2., 10., 2.))
        extent = (DensityBins[0], DensityBins[-1], TemperatureBins[0], TemperatureBins[-1])  
        im = ax.imshow(Q1[iz, :, iZZ, :], interpolation='none', origin='lower', \
               extent = extent, aspect = 'auto', vmin = np.float(PlotDict['cmin']), vmax = np.float(PlotDict['cmax']))  
        ax.autoscale(False)
        ax.plot(DensityBins, Teq[iz, iZZ, :], color = 'white', linestyle = '--')    
        cb = plt.colorbar(im,ax=ax, orientation='horizontal', pad = 0.25)
        cb.set_label(PlotDict['label1'])         
        
        ax = plt.subplot(gs[1]) 
        ax.set_title(PlotDict['top2D2'])
        ax.set_xlabel('log n [cm$^{-3}$]')
        ax.set_ylabel('')
        ax.xaxis.set_ticks(np.arange(DensityBins[0], DensityBins[-1]+2., 2.))
        ax.yaxis.set_ticks(np.arange(2., 10., 2.))
        extent = (DensityBins[0], DensityBins[-1], TemperatureBins[0], TemperatureBins[-1])  
        im = ax.imshow(Q2[iz, :, iZZ, :], interpolation='none', origin='lower', \
               extent = extent, aspect = 'auto', vmin = np.float(PlotDict['cmin']), vmax = np.float(PlotDict['cmax']))  
        ax.autoscale(False)
        ax.plot(DensityBins, Teq[iz, iZZ, :], color = 'white', linestyle = '--')    
        cb = plt.colorbar(im,ax=ax, orientation='horizontal', pad = 0.25)
        cb.set_label(PlotDict['label2'])        

        ax = plt.subplot(gs[2]) 
        ax.set_title(PlotDict['top2D3'])
        ax.set_xlabel('log n [cm$^{-3}$]')
        ax.set_ylabel('')
        ax.xaxis.set_ticks(np.arange(DensityBins[0], DensityBins[-1]+2., 2.))
        ax.yaxis.set_ticks(np.arange(2., 10., 2.))
        extent = (DensityBins[0], DensityBins[-1], TemperatureBins[0], TemperatureBins[-1])  
        im = ax.imshow(Q3[iz, :, iZZ, :], interpolation='none', origin='lower', \
               extent = extent, aspect = 'auto', vmin = np.float(PlotDict['cmin']), vmax = np.float(PlotDict['cmax']))  
        ax.autoscale(False)
        ax.plot(DensityBins, Teq[iz, iZZ, :], color = 'white', linestyle = '--')    
        cb = plt.colorbar(im,ax=ax, orientation='horizontal', pad = 0.25)
        cb.set_label(PlotDict['label3'])   
        
    else: 
        fig.subplots_adjust(left = 0.2, right = 0.95, bottom = 0.2, top = 0.8)
        gs = gridspec.GridSpec(1,1)   
        ax = plt.subplot(gs[0]) 
        ax.set_title(PlotDict['name'])
    
        if PlotType == 0: # constant density
            ax.set_xlabel('log T [K]')
            ax.set_ylabel(PlotDict['label'])
            ax.set_ylim(np.float(PlotDict['cmin']), np.float(PlotDict['cmax']))
            ax.plot(TemperatureBins, Q1[iz,:,iZZ,idens], label = PlotDict['top2D1'])
            ax.plot(TemperatureBins, Q2[iz,:,iZZ,idens], label = PlotDict['top2D2'])
            ax.plot(TemperatureBins, Q3[iz,:,iZZ,idens], label = PlotDict['top2D3'])
        elif PlotType == 1: # equilibrium temperature
            ax.set_xlabel('log n [cm$^{-3}$]')
            ax.set_ylabel(PlotDict['label'])
            ax.set_ylim(np.float(PlotDict['cmin']), np.float(PlotDict['cmax']))

            ax.plot(DensityBins, Q1[iz,iZZ,:], label = PlotDict['top2D1'])
            ax.plot(DensityBins, Q2[iz,iZZ,:], label = PlotDict['top2D2'])
            ax.plot(DensityBins, Q3[iz,iZZ,:], label = PlotDict['top2D3'])
   
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc=3, fontsize = 12)
    
    if idisplay == 1:
        outputfile = "tmp.png"
    
    fig.savefig(outputfile, dpi = 100)
    print('Plot saved as: %s'%(outputfile))
    plt.close('all')

    if idisplay == 1:
        os.system("display %s &"%(outputfile))  
    
    
    
    
    
