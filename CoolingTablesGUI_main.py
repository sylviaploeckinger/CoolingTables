#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 17:27:25 2018

@author: sylviaploeckinger
"""

#####################################
# Import standard python libraries
#####################################
import warnings
warnings.filterwarnings('ignore') 
import matplotlib
matplotlib.use('TkAgg')
import sys
import h5py
import os
import numpy as np
from tkinter import *
from tkinter import messagebox

#####################################
# Import local python files 
#####################################
from global_parameters import runname, cwdpath
from alldictnames import *
from gui_Cool import gui_cooling_plots
from gui_Heat import gui_heating_plots
from gui_CoolHeat2D import gui_coolheat2D
from gui_Generic import gui_generic_plots
from gui_Generic_2panels import gui_generic_2panelplots
from gui_Generic_3panels import gui_generic_3panelplots
from gui_Generic_IonFractions import  gui_ionfraction_plots
from gui_Generic_ElecFractions import gui_elecfraction_plots
from gui_Generic_Depletion import gui_depletion_plots
from gui_Generic_HydExt import gui_hydext_plots
#from gui_CoolingTime2D import gui_coolingtime2D

#####################################
# Path to the cooling table
#####################################
filebase = sys.argv[1]

#####################################
# Create folder for plots
#####################################
# plots are stored in cwd/plots by default:
print ("All plots will be stored in: ")
print ("%s/plots"%(cwdpath))
if not os.path.exists(cwdpath+'/plots/'):
    os.makedirs(cwdpath+'/plots/')

#####################################
# Get bins of the cooling table 
#####################################
def loadparams():
    global RedshiftBins
    global MetallicityBins
    global TemperatureBins
    global DensityBins     
    print (filebase)
    print ("Reading: ", runname[irunselect])
    myhdf5file = '%s/%s.hdf5'%(filebase, runname[irunselect])
   
    with h5py.File(myhdf5file, "r") as f:
         RedshiftBins       = f['TableBins/RedshiftBins'][()]
         MetallicityBins    = f['TableBins/MetallicityBins'][()]
         TemperatureBins    = f['TableBins/TemperatureBins'][()]
         DensityBins        = f['TableBins/DensityBins'][()]
         
irunselect = 0
loadparams()


#####################################
# Set gui item properties 
# (sizes and positions)
#####################################
wb = 8    # width of plot buttons
wx = 800  # width of window
wy = 600  # height of window

lm = 50   # left margin
rm = 50   # right margin 
cm = 20   # central margin 
tm = 50   # top margin

syframe1 =  80     # height of select table frame
syframe2 = 120     # height of select redshift frame
syframe3 = 120     # height of select metallicity frame 
syframe4 = 120     # height of select density frame
syframe5 = syframe4
syframe6 = syframe2+cm+syframe3+cm+syframe1

sl = 200.    # length of scale bar

fx = (wx - (lm + cm + rm))/2

#####################################
# Initialize gui 
#####################################
root = Tk()
root.title("Cooling Tables")
root.geometry("%ix%i"%(wx, wy))

############################################################################################
# Create and place the frames 
############################################################################################

TableFrame       = Frame(root, bd=5, width=fx, height=syframe1, relief=RIDGE)
RedshiftFrame    = Frame(root, bd=5, width=fx, height=syframe2, relief=RIDGE)
MetallicityFrame = Frame(root, bd=5, width=fx, height=syframe3, relief=RIDGE)
DensityFrame     = Frame(root, bd=5, width=fx, height=syframe4, relief=RIDGE)
InfoFrame        = Frame(root, bd=5, width=fx, height=syframe5, relief=RIDGE)
PlotFrame        = Frame(root, bd=5, width=fx, height=syframe6, relief=RIDGE)
MainPlotFrame    = Frame(PlotFrame, width=fx-30, height=syframe6-60)

TableFrame.place      (x=rm, y=tm, anchor="nw")
RedshiftFrame.place   (x=rm, y=tm+cm+syframe1, anchor="nw")
MetallicityFrame.place(x=rm, y=tm+cm+syframe1+cm+syframe2, anchor="nw")
DensityFrame.place    (x=rm, y=tm+cm+syframe1+cm+syframe2+cm+syframe3, anchor="nw")
InfoFrame.place       (x=rm+fx+cm, y=tm+cm+syframe6, anchor="nw")
PlotFrame.place       (x=rm+fx+cm, y=tm, anchor="nw")
MainPlotFrame.place   (x=10, y=30, anchor="nw")

TableFrame.grid_propagate(0)
RedshiftFrame.grid_propagate(0)
MetallicityFrame.grid_propagate(0)
DensityFrame.grid_propagate(0)
InfoFrame.grid_propagate(0)
PlotFrame.grid_propagate(0)
MainPlotFrame.grid_propagate(0)

Slabel = Label(root, text="S. Ploeckinger (2020)")
Slabel.place(rely=0.98, relx=0.98, x=0, y=0, anchor=SE)

############################################################################################
# 1. Select table set
############################################################################################
def change_dropdown(*args):
    #print( runvar.get() )
    irunselect = runname.index(runvar.get())
    #print(irunselect,  runvar.get(), runname[irunselect] )
    loadparams()

runvar = StringVar(root)

rundic = {}
for irun in range(len(runname)):
    rundic[runname[irun]] = irun
runvar.set(runname[0])

l = Label(TableFrame, text="1. Select set of tables")
popupMenu = OptionMenu(TableFrame, runvar, *rundic)

l.grid(row = 1, column = 0, sticky=W)
popupMenu.grid(row = 2, column = 0, padx = 40, pady = 10)
runvar.trace('w', change_dropdown)

############################################################################################
# 2. Select redshift 
############################################################################################
indxred    = np.where(RedshiftBins < 49.)
#reionizvar = IntVar()

lred     = Label(RedshiftFrame, text="2. Select redshift: z")
if (len(indxred[0]) > 1):
    wred     = Scale(RedshiftFrame, from_=RedshiftBins[indxred[0][0]], to=RedshiftBins[indxred[0][-1]], length = sl, orient=HORIZONTAL, resolution = RedshiftBins[indxred[0][1]] - RedshiftBins[indxred[0][0]])
#cred     = Checkbutton(RedshiftFrame, text = "Before Reionization", variable = reionizvar)

lred.grid(row = 1, column = 0, sticky = W)
if (len(indxred[0]) > 1):
    wred.grid(row = 2, column = 0, sticky = E, padx = 50)
#cred.grid(row = 3, column = 0, sticky = W, pady = 20)

############################################################################################
# 3. Select metallicity 
############################################################################################
indxmet = np.where(MetallicityBins > -49.)
primvar = IntVar()

lmet     = Label(MetallicityFrame, text="3. Select metallicity: log Z/Zsol")
if (len(indxmet[0]) > 1):
    wmet     = Scale(MetallicityFrame, from_=MetallicityBins[indxmet[0][0]], to=MetallicityBins[indxmet[0][-1]], length = sl, orient=HORIZONTAL, resolution = MetallicityBins[indxmet[0][1]] - MetallicityBins[indxmet[0][0]])
cmet     = Checkbutton(MetallicityFrame, text = "Primoridal abundances", variable = primvar)

lmet.grid(row = 1, column = 0, sticky = W)
if (len(indxmet[0]) > 1):
    wmet.grid(row = 2, column = 0, sticky = E, padx = 50)
cmet.grid(row = 3, column = 0, sticky = W, pady = 20)

#####################################################
# 4. Select density
#####################################################
indxden = np.where(DensityBins > -49.)

lden     = Label(DensityFrame, text="4. Select density: log nH [cm-3]")
wden     = Scale(DensityFrame, from_=DensityBins[indxden[0][0]], to=DensityBins[indxden[0][-1]], length = sl, orient=HORIZONTAL, resolution = DensityBins[indxden[0][1]] - DensityBins[indxden[0][0]])

lden.grid(row = 1, column = 0, sticky = W)
wden.grid(row = 2, column = 0, sticky = E, padx = 50)

#####################################################
# 5. Select plot
#####################################################
def getsettings():
    global irunselect
    global ired, imet, iden
  
    irunselect = runname.index(runvar.get())
    iden = (np.abs(DensityBins - wden.get())).argmin()
    
    if primvar.get():
        imet  = (np.abs(MetallicityBins - (-50.))).argmin()
        tmet.set  ("Metallicity : primoridal")
    else:
        if (len(indxmet[0]) > 1):
            imet = (np.abs(MetallicityBins - wmet.get())).argmin()
            tmet.set  ("Metallicity : %.2f [log Z/Zsol]"%(MetallicityBins[imet]))
        elif (len(indxmet[0]) == 1):
            imet = indxmet[0][0]
            tmet.set  ("Metallicity : %.2f [log Z/Zsol]"%(MetallicityBins[imet]))
        else:
            imet  = (np.abs(MetallicityBins - (-50.))).argmin()
            tmet.set  ("Metallicity : primoridal")
            
    if (len(indxred[0]) > 1):
        ired = (np.abs(RedshiftBins - wred.get())).argmin()
        tred.set  ("Redshift    : %.2f [z]"%(RedshiftBins[ired]))
    elif (len(indxred[0]) == 1):
        ired = indxred[0][0]
        tred.set  ("Redshift    : %.2f [z]"%(RedshiftBins[ired]))
    else:
        ired = (np.abs(RedshiftBins - 50.)).argmin()
        tred.set  ("Redshift    : before reioniz")            
        
    ttable.set("Table set   : %s"%(runname[irunselect]))
    tden.set  ("Density     : %.2f [log nH/cm-3]"%(DensityBins[iden]))


def genericplots_display():
    genericplots(idisplay = 1)

def genericplots_save():
    genericplots(idisplay = 0)

def genericplots(idisplay):
    getsettings()
    SelectPlotTypeVar.get()
    iplottype = SelectPlotTypeVar.get()
    idsettype = plotvar.get()
    plotvar.get()
    if idsettype in dict_names_simple:
########### gui_generic_plots                from: gui_Generic.py
        gui_generic_plots(irunselect, ired, imet, iden, iplottype, dict_names[idsettype], idisplay)
    elif idsettype in dict_names_2panels:
########### gui_generic_2panelplots          from: gui_Generic_2panels.py
        gui_generic_2panelplots(irunselect, ired, imet, iden, iplottype, dict_names[idsettype], idisplay)
    elif idsettype in dict_names_3panels:
########### gui_generic_3panelplots          from: gui_Generic_3panels.py  
        gui_generic_3panelplots(irunselect, ired, imet, iden, iplottype, dict_names[idsettype], idisplay)
#    elif idsettype in dict_names_tcool:
#        if iplottype == 2:
########### gui_coolingtime2D                   from: gui_CoolintTime2D.py
#            gui_coolingtime2D(irunselect, ired, imet, idisplay)
    elif idsettype in dict_names_cool:
        if iplottype == 2:
########### gui_coolheat2D                   from: gui_CoolHeat2D.py
            gui_coolheat2D(irunselect, ired, imet, idisplay)
        else:
########### gui_cooling_plots                from: gui_Cool.py
            gui_cooling_plots(irunselect, ired, imet, iden, iplottype, idisplay)        
    elif idsettype in dict_names_heat:
        if iplottype == 2:
########### gui_coolheat2D                   from: gui_CoolHeat2D.py
            gui_coolheat2D(irunselect, ired, imet, idisplay)
        else:
########### gui_ionfraction_plots            from: gui_Heat.py
            gui_heating_plots(irunselect, ired, imet, iden, iplottype, idisplay)
    elif idsettype in dict_names_ionfrac:     
########### gui_ionfraction_plots            from: gui_Generic_IonFractions.py
        gui_ionfraction_plots(irunselect, ired, imet, iden, iplottype, dict_names[idsettype], idisplay)
    elif idsettype in dict_names_elec:
########### gui_elecfraction_plots           from: gui_Generic_ElecFractions.py
        gui_elecfraction_plots(irunselect, ired, imet, iden, iplottype, dict_names[idsettype], idisplay)
    elif idsettype in dict_names_depl:
########### gui_depletion_plots              from: gui_Generic_Depletion
        gui_depletion_plots(irunselect, ired, imet, iden, iplottype, dict_names[idsettype], idisplay)
    elif idsettype in dict_names_hydext:
########### gui_hydext_plots                 from: gui_Generic_HydExt
        gui_hydext_plots(irunselect, ired, imet, iden, iplottype, dict_names[idsettype], idisplay )
    else:
        print('Choice not implemented')
        
        
lplot = Label(PlotFrame, text="5. Select plot type")
lplot.grid(row = 0, column = 0, sticky = W)


def change_dropdown2(*args):
    subdict = dict_names[plotvar.get()]

plotvar = StringVar(root)
plotvar.set(dict_names['Cosmic ray rate']['name'])

lpop2    = Label     (MainPlotFrame, text="Test")
pop2Menu = OptionMenu(MainPlotFrame, plotvar, *dict_names)

pop2Menu.grid(row = 5, column = 0, sticky = W, padx = 0, pady = 10)
plotvar.trace('w', change_dropdown2)

SelectPlotTypeVar = IntVar()
RB2D       = Radiobutton(MainPlotFrame, text="2D (for set metallicity, redshift)", variable=SelectPlotTypeVar, value=2)
RBalongn   = Radiobutton(MainPlotFrame, text="1D (constant density)"               , variable=SelectPlotTypeVar, value=0)
RBalongTeq = Radiobutton(MainPlotFrame, text="1D (equilibrium temperature)"        , variable=SelectPlotTypeVar, value=1)
RBalongn.grid  (row = 7, column = 0, sticky = W, pady = 20)
RBalongTeq.grid(row = 8, column = 0, sticky = W, pady = 20)
RB2D.grid      (row = 6, column = 0, sticky = W, pady = 20)

display_button  = Button(MainPlotFrame, text="Display", command=genericplots_display, width = wb)
saveplot_button  = Button(MainPlotFrame, text="Save", command=genericplots_save, width = wb)

display_button.grid(row = 9, column = 0, sticky = W)
saveplot_button.grid(row = 9, column = 1, sticky = W)


#####################################################
# 6. Info box
#####################################################

linfoheader = Label(InfoFrame, text = 'Selection overview:')

ttable = StringVar()
tred = StringVar()
tmet = StringVar()
tden = StringVar()


ltable = Label(InfoFrame, textvariable = ttable)
lred   = Label(InfoFrame, textvariable = tred)
lmet   = Label(InfoFrame, textvariable = tmet)
lden   = Label(InfoFrame, textvariable = tden)

ttable.set("Table set   : ")
tred.set  ("Redshift    : ")
tmet.set  ("Metallicity : ")
tden.set  ("Density     : ")

linfoheader.grid(row = 0, column = 0, sticky = W)
ltable.grid(row = 1, column = 0, sticky = W)
lred.grid  (row = 2, column = 0, sticky = W)
lmet.grid  (row = 3, column = 0, sticky = W)
lden.grid  (row = 4, column = 0, sticky = W)

root.mainloop()
