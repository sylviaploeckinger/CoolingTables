#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 16:35:13 2018

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
def gui_heating_plots(irun,iz,iZZ,idens,iplottype, idisplay):
    if iplottype > 1:
        return
    
    ymin = -35.
    ymax = -18.
    
    ########################################
    # Split the heating contributions into:
    # panel 1: primoridal elements
    # panel 2: atomic metal heating
    # panel 3: remaining heating channels
    # This split is hardcoded here
    ########################################
    primheat = [0,1,12,14,17,19,20,22]
    atomheat = [2,3,4,5,6,7,8,9,10,23]
    restheat = [11,13,15,16,18,21,23]
    
    myhdf5file = '%s/%s.hdf5'%(filebase, runname[irun])
    with h5py.File(myhdf5file, "r") as f:
        RedshiftBins       = f['TableBins/RedshiftBins'].value
        MetallicityBins    = f['TableBins/MetallicityBins'].value
        TemperatureBins    = f['TableBins/TemperatureBins'].value
        DensityBins        = f['TableBins/DensityBins'].value
        IdentifierHeating = f['IdentifierHeating'].value
            
        if iplottype == 0:
            Heating = f['Tdep/Heating'].value
            Heat1D  = Heating[iz, :, iZZ, idens, :]          
            xx      = TemperatureBins
            xlab    = 'log T [K]'
        else:
            Heating = f['ThermEq/Heating'].value
            Heat1D  = Heating[iz, iZZ, :, :]
            xx      = DensityBins
            xlab    = 'log n [cm$^{-3}$]'

    xmin = xx[0]
    xmax = xx[-1]            

    titlestring, outputfile = getfilenames(runname[irun], RedshiftBins[iz], DensityBins[idens],\
                                             MetallicityBins[iZZ], iz, iZZ, idens, iplottype, 'Heat')
    
    print ('------------------------------------------------------------------')
    print ('Outputfile: ')
    print (outputfile)
    print ('------------------------------------------------------------------')


    fig = plt.figure()
    fig.set_size_inches(10,6.2,forward=True)
    fig.suptitle(titlestring)
    fig.subplots_adjust(left = 0.1, right = 0.95, bottom = 0.15, top = 0.7)
    gs = gridspec.GridSpec(1,3,wspace=0, hspace=0)

    ax = plt.subplot(gs[0])
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)

    ax.set_xlabel(xlab)
    ax.set_ylabel('log $\Gamma$/$n_{\mathrm{H}}^2$ [erg cm$^{3}$ s$^{-1}$]')

    ax.xaxis.set_ticks(np.arange(xx[0], xx[-2], 1.))
    ax.plot(xx, np.log10(np.power(10., Heat1D[..., -2]) + np.power(10., Heat1D[..., -1])), color = 'grey', lw = 2, label = 'Total')
    icount = 0
    c = cmap(np.linspace(0,1,len(primheat)-1))
    
    for i in range (len(primheat)):
        iheat = primheat[i]
        if iheat == primheat[-1]:
            linec = 'black'
            lines = ':'
        else:
            linec = c[icount]
            if icount%3 == 0:
                linec = c[icount]
                lines = '-'
            elif icount%3 == 1:
                lines = '--'
            else: 
                lines = '-.'

        ax.plot(xx, Heat1D[..., iheat], color = linec, lw = 2, ls = lines, label = '%s'%(IdentifierHeating[iheat].decode('utf-8')))
        icount = icount + 1


    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, bbox_to_anchor=(0., 1.02, 1., .102), loc=3, mode = 'expand', borderaxespad=0., ncol = 2, fontsize = VERYSMALL_SIZE, handlelength = 4)


    ax = plt.subplot(gs[1])
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    ax.tick_params(labelleft=False)

    ax.set_xlabel(xlab)
    ax.xaxis.set_ticks(np.arange(xx[0], xx[-2], 1.))
    ax.plot(xx, np.log10(np.power(10., Heat1D[..., -2]) + np.power(10., Heat1D[..., -1])), color = 'grey', lw = 2, label = 'Total')
    icount = 0
    c = cmap(np.linspace(0,1,len(atomheat)-1))

    for i in range (len(atomheat)):
        iheat = atomheat[i]
        if iheat == atomheat[-1]:
            linec = 'black'
            lines = ':'
        else:
            linec = c[icount]
            if icount%3 == 0:
                linec = c[icount]
                lines = '-'
            elif icount%3 == 1:
                lines = '--'
            else: 
                lines = '-.'
                
        ax.plot(xx, Heat1D[..., iheat], color = linec, lw = 2, ls = lines, label = '%s'%(IdentifierHeating[iheat].decode('utf-8')))
        icount = icount + 1


    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, bbox_to_anchor=(0., 1.02, 1., .102), loc=3, mode = 'expand', borderaxespad=0., ncol = 2, fontsize = VERYSMALL_SIZE, handlelength = 4)


    ax = plt.subplot(gs[2])
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    ax.tick_params(labelleft=False)

    ax.set_xlabel(xlab)
    ax.xaxis.set_ticks(np.arange(xx[0], xx[-2], 1.))

    ax.plot(xx, np.log10(np.power(10., Heat1D[..., -2]) + np.power(10., Heat1D[..., -1])), color = 'grey', lw = 2, label = 'Total')
    icount = 0
    c = cmap(np.linspace(0,1,len(restheat)-1))
    for i in range (len(restheat)):
        iheat = restheat[i]
        if iheat == restheat[-1]:
            linec = 'black'
            lines = ':'
        else:
            linec = c[icount]
            if icount%3 == 0:
                linec = c[icount]
                lines = '-'
            elif icount%3 == 1:
                lines = '--'
            else: 
                lines = '-.'

        ax.plot(xx, Heat1D[..., iheat], color = linec, lw = 2, ls = lines, label = '%s'%(IdentifierHeating[iheat].decode('utf-8')))
        icount = icount + 1


    handles, labels = ax.get_legend_handles_labels()
    labels = [w.replace('Oatoms', 'OtherA') for w in labels]
    ax.legend(handles, labels, bbox_to_anchor=(0., 1.02, 1., .102), loc=3, mode = 'expand', borderaxespad=0., ncol = 2, fontsize = VERYSMALL_SIZE, handlelength = 4)

    fig.savefig(outputfile, dpi = 100)
    plt.close('all')

    if idisplay == 1:
        os.system("display %s &"%(outputfile))
