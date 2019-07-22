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
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cycler import cycler

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
lines = ['-', '-.', '--', ':']
colorcycle = ['b','r','g', 'y']
colorcycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

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
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIG_SIZE)  # fontsize of the figure title

cmap = plt.cm.get_cmap('viridis')
props = dict(boxstyle='round', facecolor='white')

def gui_hydext_plots(irun, iz, iZZ, idens, PlotType, PlotDict, idisplay):
    myhdf5file = '%s/%s.hdf5'%(filebase, runname[irun])
    with h5py.File(myhdf5file, "r") as f:
        RedshiftBins       = f['TableBins/RedshiftBins'].value
        MetallicityBins    = f['TableBins/MetallicityBins'].value
        TemperatureBins    = f['TableBins/TemperatureBins'].value
        DensityBins        = f['TableBins/DensityBins'].value            

        if PlotType == 0:         # constant density 
            xlab = 'log T [K]'
            ylab = PlotDict['label']
            xx   = TemperatureBins
            Q = f['Tdep/'+PlotDict['dset']].value
            Q1D  = Q[iz,:,iZZ,idens,:]
            
        if PlotType == 1:          # Thermal equilibrium   
            xlab = 'log n [cm$^{-3}$]'
            ylab = PlotDict['label']
            xx   = DensityBins
            Q   = f['ThermEq/'+PlotDict['dset']].value
            Q1D  = Q[iz,iZZ,:,:]
            
        if PlotType == 2:           # 2D
            Q = f['Tdep/'+PlotDict['dset']].value
            Q2D = Q[iz,:,iZZ,:,:]
            Teq = f['ThermEq/Temperature'].value   
            
            
    ymin = -10.8
    ymax =  0.5     
                                                   
                                                   
    top2D = []
    top2D.append(PlotDict['top2D1'])                                               
    top2D.append(PlotDict['top2D2'])                                               
    top2D.append(PlotDict['top2D3'])                                               
    top2D.append(PlotDict['top2D4'])                                               
    top2D.append(PlotDict['top2D5'])                                               
    top2D.append(PlotDict['top2D6'])                                               

    lab = []
    lab.append(PlotDict['label1'])
    lab.append(PlotDict['label2'])
    lab.append(PlotDict['label3'])
    lab.append(PlotDict['label4'])
    lab.append(PlotDict['label5'])
    lab.append(PlotDict['label6'])

            
    titlestring, outputfile = getfilenames(runname[irun], RedshiftBins[iz], DensityBins[idens],\
                                             MetallicityBins[iZZ], iz, iZZ, idens, PlotType, PlotDict['shortname'])
    
    fig = plt.figure()
    if PlotType < 2:  # 1D plots
        fig.set_size_inches(7,4,forward=True)

        fig.suptitle(titlestring, fontsize = 12)
        fig.subplots_adjust(left = 0.15, right = 0.75, bottom = 0.2, top = 0.8)
        gs = gridspec.GridSpec(1,1)
        ax = plt.subplot(gs[0]) 
        ax.set_title(PlotDict['name'])
    
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_ylim(ymin, ymax)
        ax.xaxis.set_ticks(np.arange(xx[0], xx[-1]+1., 1.))
        ax.set_prop_cycle(cycler('linestyle', lines)*#cycler(color='bgrcmyk'))
        cycler('color', colorcycle))


        for iion in range (len(Q1D[0,:])):            
            stringlabel = '%s'%(lab[iion])
            ax.plot(xx, Q1D[:,iion], label = stringlabel)       
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, bbox_to_anchor=(1.02, 0., 1.2, 1), loc=3, borderaxespad=0., ncol = 1, fontsize = 12)
        
    else:  # 2D
        extent = (DensityBins[0], DensityBins[-1], TemperatureBins[0], TemperatureBins[-1])  

        iix = np.int8(np.float(PlotDict['nrxpanels']))
        iiy = np.int8(np.float(PlotDict['nrypanels']))
        sx  = 1.8   # inch
        sy  = 1.3   # inch
        lm  = sx/2. # inch
        rm  = sx    # inch
        tm  = sy/2. # inch
        bm  = sy/2. # inch
        cm  = sx/4. # inch 
        cmin = ymin
        cmax = ymax
        hsp  = 0. # inch
        
        xsize_inch = np.float(iix) * sx + lm + rm
        ysize_inch = np.float(iiy) * sy + tm + bm + hsp
        fig.set_size_inches(xsize_inch,ysize_inch,forward=True)

        fig.suptitle(titlestring, fontsize = 12)
        fig.subplots_adjust(left = lm/xsize_inch, right = 1., bottom = bm/ysize_inch, top = 1.-tm/ysize_inch)

        gs  = gridspec.GridSpec(iiy, iix, wspace = 0 , hspace = hsp/ysize_inch, right = 1.-rm/xsize_inch)
        gs1 = gridspec.GridSpec(1,1, wspace = 0, hspace = 0, left = 1.-rm/xsize_inch, right = 1.-(rm-cm)/xsize_inch)

        iionmax = len(Q2D[0,0,:])

        for iion in range (iionmax):
            ax = plt.subplot(gs[iion])
            ax.tick_params(direction = 'in')
            ax.xaxis.set_ticks(np.arange(DensityBins[0]+2, DensityBins[-1]+2., 2.))
            ax.yaxis.set_ticks(np.arange(TemperatureBins[0]+1., TemperatureBins[-1]+2., 2.))
            
            stringlabel = '%s'%(lab[iion])
                
            ax.text(0.05, 0.95, stringlabel, transform=ax.transAxes, fontsize=14,
                        verticalalignment='top', horizontalalignment = 'left', bbox=props)
            
            if iion == 0:
                ax.text(0.05, 1.05, '%s'%PlotDict['name'], transform=ax.transAxes, fontsize=14,
                        verticalalignment='bottom', horizontalalignment = 'left')

            if iion%iix == 0:
                ax.set_ylabel('log T [K]')
            else:
                ax.tick_params(labelleft=False)
                
            if iion < iionmax - iix:
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel('log n$_{\mathrm{H}}$ [cm$^{-3}$]')
                
            im = ax.imshow(Q2D[:,:,iion], interpolation='none', origin='lower', \
                           extent = extent, \
                           aspect = 'auto', vmin = cmin, vmax = cmax, cmap = cmap)
            ax.autoscale(False)
            ax.plot(DensityBins, Teq[iz, iZZ, :], color = 'white', linestyle = '--')
            
        ax = plt.subplot(gs1[0])
        ax.set_axis_off()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="200%", pad=0.05)
        cb  = plt.colorbar(im, cax=cax, orientation='vertical', ticks = np.arange(int(cmin), cmax, 2.))
        cb.set_label('%s'%(PlotDict['label']))        
    
    fig.savefig(outputfile, dpi = 100)
    plt.close('all')

    if idisplay == 1:
        os.system("display %s &"%(outputfile))    
    
    
    
    
    
    
    
    
    
