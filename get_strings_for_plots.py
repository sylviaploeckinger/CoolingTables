#####################################
# Import standard python libraries
#####################################
import sys

#####################################
# Import local python files 
#####################################
from global_parameters import cwdpath

#####################################
# Path to the cooling table
#####################################
filebase = sys.argv[1]


def getfilenames(runnameloc, z, n, ZZ, iz, iZZ, idens, PlotType, shortstring):
    
    if z > 49.:
        redstring = 'befRe, '
    else:
        redstring = 'z=%4.1f, '%(z)
    
    densstring = 'log n$_{\mathrm{H}}$ [cm$^{-3}$]=%4.1f, '%(n)
    
    if ZZ < -49.:
        metstring = 'log Z/Z$_{\odot}$=prim'
    else:
        metstring = 'log Z/Z$_{\odot}$=%4.1f'%(ZZ)
        
    if PlotType == 0: # constant density
        titlestring = runnameloc + ', ' + redstring + densstring + metstring
        outputfile = '%s/plots/%s_%s_iz%2.2i_iZZ%2.2i_idens%2.2i.png'%(cwdpath, shortstring, runnameloc, iz, iZZ, idens)
    elif PlotType == 1: # Thermal Equilibrium
        titlestring = runnameloc + ', ' + redstring + metstring
        outputfile = '%s/plots/ThermEq_%s_%s_iz%2.2i_iZZ%2.2i.png'%(cwdpath, shortstring, runnameloc, iz, iZZ)
    else:
        titlestring = runnameloc + ', ' + redstring + metstring
        outputfile = '%s/plots/%s_%s_iz%2.2i_iZZ%2.2i.png'%(cwdpath, shortstring, runnameloc, iz, iZZ)
    return titlestring, outputfile
