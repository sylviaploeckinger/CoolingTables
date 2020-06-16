#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 16:43:31 2018

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
def gui_coolheat2D(irun, iz, iZZ, idisplay):
    cmin = -35.
    cmax = -18.
    cmin2= -3.
    cmax2= +3.

    myhdf5file = '%s/%s.hdf5'%(filebase, runname[irun])
    with h5py.File(myhdf5file, "r") as f:
        RedshiftBins       = f['TableBins/RedshiftBins'].value
        MetallicityBins    = f['TableBins/MetallicityBins'].value
        TemperatureBins    = f['TableBins/TemperatureBins'].value
        DensityBins        = f['TableBins/DensityBins'].value
        ThermEqT           = f['ThermEq/Temperature'].value
        Heating = f['Tdep/Heating'].value
        Cooling = f['Tdep/Cooling'].value
    extent = (DensityBins[0], DensityBins[-1], TemperatureBins[0], TemperatureBins[-1])

    Heat2D = np.zeros_like(Heating[iz, :, iZZ, :, 0])
    Heat2D = np.log10(np.power(10., Heating[iz, :, iZZ, :, -2]) + np.power(10., Heating[iz, :, iZZ, :, -1]))
    Cool2D = np.zeros_like(Cooling[iz, :, iZZ, :, 0])
    Cool2D = np.log10(np.power(10., Cooling[iz, :, iZZ, :, -2]) + np.power(10., Cooling[iz, :, iZZ, :, -1]))
    CooldivHeat      = np.zeros_like(Cooling[iz, :, iZZ, :, 0])
    CooldivHeat[...] = Cool2D[...] - Heat2D[...]

    for iden in range (len(DensityBins)):
        Heat2D[:,iden] = Heat2D[:,iden] + 2. * DensityBins[iden]
        Cool2D[:,iden] = Cool2D[:,iden] + 2. * DensityBins[iden]

    titlestring, outputfile = getfilenames(runname[irun], RedshiftBins[iz], 0.0,\
                                             MetallicityBins[iZZ], iz, iZZ, 0, 2, 'CoolHeat2D')
    
    neq = DensityBins
    Teq = ThermEqT[iz,iZZ,:]

    fig = plt.figure()
    fig.suptitle(titlestring, fontsize = 12)
    fig.subplots_adjust(left = 0.1, right = 0.95, bottom = 0.15, top = 0.90)
    gs = gridspec.GridSpec(1,3, wspace = 0.15, hspace = 0)

    ax = plt.subplot(gs[0])
    ax.set_xlabel('log n$_{\mathrm{H}}$ [cm$^{-3}$]')
    ax.set_ylabel('log T [K]')
    ax.tick_params(direction = 'inout')
    ax.xaxis.set_ticks(np.arange(DensityBins[0], DensityBins[-1]+2., 2.))

    im = ax.imshow(Cool2D, interpolation='none', origin='lower', \
               extent = extent, \
               aspect = 'auto', vmin = cmin, vmax = cmax, cmap = cmap)
    
    ax.autoscale(False)
    ax.plot(neq, Teq, color = 'white', linestyle = '--')   
    
    cb = plt.colorbar(im,ax=ax, orientation='horizontal', pad = 0.2, ticks = np.arange(cmin, cmax, 5.))
    cb.set_label('log $\Lambda$ \n[erg cm$^{-3}$ s$^{-1}$]')

    ax = plt.subplot(gs[1])
    ax.set_xlabel('log n$_{\mathrm{H}}$ [cm$^{-3}$]')
    ax.tick_params(labelleft=False, direction = 'inout')
    ax.xaxis.set_ticks(np.arange(DensityBins[0], DensityBins[-1]+2., 2.))

    im = ax.imshow(Heat2D, interpolation='none', origin='lower', \
               extent = extent, \
               aspect = 'auto', vmin = cmin, vmax = cmax, cmap = cmap)

    ax.autoscale(False)
    ax.plot(neq, Teq, color = 'white', linestyle = '--')

    cb = plt.colorbar(im,ax=ax, orientation='horizontal', pad = 0.2, ticks = np.arange(cmin, cmax, 5.))
    cb.set_label('log $\Gamma$ \n [erg cm$^{-3}$ s$^{-1}$]')

    ax = plt.subplot(gs[2])
    ax.set_xlabel('log n$_{\mathrm{H}}$ [cm$^{-3}$]')
    ax.tick_params(labelleft=False, direction = 'inout')
    ax.xaxis.set_ticks(np.arange(DensityBins[0], DensityBins[-1]+2., 2.))

    im = ax.imshow(CooldivHeat, interpolation='none', origin='lower', \
               extent = extent, \
               aspect = 'auto', vmin = cmin2, vmax = cmax2, cmap = 'RdBu')

    ax.autoscale(False)
    ax.plot(neq, Teq, color = 'white', linestyle = '--')

    cb = plt.colorbar(im,ax=ax, orientation='horizontal', pad = 0.2)
    cb.set_label('log $\Lambda / \Gamma$')

    '''
    handles, labels = ax.get_legend_handles_labels()
    l = ax.legend(handles, labels, loc=1, fontsize = 12)
    for text in l.get_texts():
        text.set_color("white")
    l.get_frame().set_facecolor('none')
    '''

    if idisplay == 1:
        outputfile = "tmp.png"
    
    fig.savefig(outputfile, dpi = 100)
    print('Plot saved as: %s'%(outputfile))
    plt.close('all')

    if idisplay == 1:
        os.system("display %s &"%(outputfile))  
