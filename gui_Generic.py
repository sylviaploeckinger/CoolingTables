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
BIG_SIZE	= SMALL_SIZE + 4
BIGGER_SIZE = SMALL_SIZE + 6
plt.rc('font', size=SMALL_SIZE)			 # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)	 # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)	 # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)	 # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)	 # fontsize of the tick labels
plt.rc('legend', fontsize=VERYSMALL_SIZE)	 # legend fontsize
plt.rc('figure', titlesize=BIG_SIZE)  # fontsize of the figure title

cmap = plt.cm.get_cmap('viridis')
props = dict(boxstyle='round', facecolor='white')

#####################################
# Plotting routine 
#####################################
def gui_generic_plots(irun, iz, iZZ, idens, PlotType, PlotDict, idisplay):
	myhdf5file = '%s/%s.hdf5'%(filebase, runname[irun])
	with h5py.File(myhdf5file, "r") as f:
		RedshiftBins	   = f['TableBins/RedshiftBins'][:]
		MetallicityBins    = f['TableBins/MetallicityBins'][:]
		TemperatureBins    = f['TableBins/TemperatureBins'][:]
		DensityBins		   = f['TableBins/DensityBins'][:]
   
		if PlotType == 1:	# Thermal Equilibrium
			Q	= f['ThermEq/'+PlotDict['dset']][:]
		else:
			Q = f['Tdep/'+PlotDict['dset']][:]
			Nref = f['Tdep/ShieldingColumnRef'][:]
		
	titlestring, outputfile = getfilenames(runname[irun], RedshiftBins[iz], DensityBins[idens],\
											 MetallicityBins[iZZ], iz, iZZ, idens, PlotType, PlotDict['shortname'])
	
	fig = plt.figure()
	fig.suptitle(titlestring, fontsize = 12)
	fig.subplots_adjust(left = 0.2, right = 0.95, bottom = 0.2, top = 0.8)
	gs = gridspec.GridSpec(1,1)
	ax = plt.subplot(gs[0]) 
	ax.set_title(PlotDict['name'])
	
	if PlotDict['name'] == 'Grid Fails':
		print ('GridFails ----- ')
		print ('Total number of grid fails	= %i'%( np.int(np.sum(Q))))
	
	if PlotType == 0: # constant density
		ax.set_xlabel('log T [K]')
		ax.set_ylabel(PlotDict['label'])
		ax.plot(TemperatureBins, Q[iz,:,iZZ,idens])
	elif PlotType == 1: # equilibrium temperature
		ax.set_xlabel(r'log n$_{\mathrm{H}}$ [cm$^{-3}$]')
		ax.xaxis.set_ticks(np.arange(DensityBins[0], DensityBins[-1] + 2., 2.))
		ax.set_ylabel(PlotDict['label'])
		ax.plot(DensityBins, Q[iz,iZZ,:])
	else:  # 2D
		ax.set_xlabel(r'log n$_{\mathrm{H}}$ [cm$^{-3}$]')
		ax.xaxis.set_ticks(np.arange(DensityBins[0], DensityBins[-1] + 2., 2.))
		ax.set_ylabel('log T [K]')
		extent = (DensityBins[0], DensityBins[-1], TemperatureBins[0], TemperatureBins[-1]) 
		if PlotDict['label'] == 'Fails':
			im = ax.imshow(Q[iz, :, iZZ, :], interpolation='none', origin='lower', \
			   extent = extent, aspect = 'auto', vmin = 0, vmax = 1.)  
		elif PlotDict['shortname'] == 'DG':
				im = ax.imshow(Q[iz, :, iZZ, :], interpolation='none', origin='lower', \
							 extent = extent, aspect = 'auto', vmin = -10., vmax = -2.)
				if (MetallicityBins[iZZ] == int(MetallicityBins[iZZ])):
					stringlabel = "log Z/Z$_{\odot}$ = %.0f"%(MetallicityBins[iZZ])
				else:
					stringlabel = "log Z/Z$_{\odot}$ = %.1f"%(MetallicityBins[iZZ])
				ax.text(0.05, 0.95, stringlabel, transform=ax.transAxes, fontsize=14, \
						 verticalalignment='top', horizontalalignment = 'left', bbox=props)
				ax.contour(DensityBins, TemperatureBins, Nref[iz, :, iZZ, :], levels=[20.56], \
						   colors = 'white', linestyles = 'dashed', origin = 'lower', \
						   linewidths = 2)
				ax.text(2., 4.2, 'N$_{\mathrm{ref}}$ = N$_{\mathrm{H,0}}$', color = 'white')
		elif PlotDict['shortname'] == 'Rad':
				im = ax.imshow(Q[iz, :, iZZ, :], interpolation='none', origin='lower', \
							   extent = extent, aspect = 'auto', vmin = -3., vmax = 4.) 
				if (RedshiftBins[iz] == int(RedshiftBins[iz])):
					stringlabel = "z = %.0f"%(RedshiftBins[iz])
				else:
					stringlabel = "z = %.1f"%(RedshiftBins[iz])
				ax.text(0.05, 0.95, stringlabel, transform=ax.transAxes, fontsize=14, \
						 verticalalignment='top', horizontalalignment = 'left', bbox=props) 
				ax.contour(DensityBins, TemperatureBins, Nref[iz, :, iZZ, :], levels=[20.56], \
						   colors = 'white', linestyles = 'dashed', origin = 'lower', \
						   linewidths = 2)
				ax.text(2., 4.2, 'N$_{\mathrm{ref}}$ = N$_{\mathrm{H,0}}$', color = 'white')
		else:
				im = ax.imshow(Q[iz, :, iZZ, :], interpolation='none', origin='lower', \
							   extent = extent, aspect = 'auto')				   
		ax.autoscale(False)
		cb = plt.colorbar(im,ax=ax, orientation='vertical', pad = 0.1)
		cb.set_label(PlotDict['label'])    

	if idisplay == 1:
		outputfile = "tmp.png"
    
	fig.savefig(outputfile, dpi = 100)
	print('Plot saved as: %s'%(outputfile))
	plt.close('all')
    
	if idisplay == 1:
		os.system("open %s &"%(outputfile))  
