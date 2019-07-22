#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 14:16:40 2018

@author: sylviaploeckinger
"""

def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z

def merge_three_dicts(x,y,z):
    xy = merge_two_dicts(x,y)
    xyz = merge_two_dicts(xy,z)
    return xyz
 
def merge_all_dicts(a,b,c,d,e,f,g,h,i):
    abcdic    = merge_three_dicts(a,b,c)
    defdic    = merge_three_dicts(d,e,f)
    abcdefg   = merge_three_dicts(abcdic, defdic, g)
    abcdefghi = merge_three_dicts(abcdefg, h, i)
    return abcdefghi


# build the dictionary
dict_names_simple = {'Cosmic ray rate'          : {'shortname':'CR', 
                                                   'name'     :'Cosmic ray rate', 
                                                   'dset'     :'CosmicRayRate', 
                                                   'label'    :'log CR rate [s$^{-1}$]'},
    
                     'D/G ratio'                : {'shortname':'DG', 
                                                   'name'     :'D/G ratio', 
                                                   'dset'     :'DGratio', 
                                                   'label'    :'log D/G'}, 
                                                   
                     'Gamma'                    : {'shortname':'Gamma', 
                                                   'name'     :'Gamma', 
                                                   'dset'     :'GammaHeat', 
                                                   'label'    :'$\gamma$'},
                                                   
                     'Fails'                    : {'shortname':'Fails', 
                                                   'name'     :'Grid Fails', 
                                                   'dset'     :'GridFails', 
                                                   'label'    :'Fails'},
                                                   
                     'Mean particle mass'       : {'shortname':'Mu', 
                                                   'name'     :'Mean particle mass', 
                                                   'dset'     :'MeanParticleMass', 
                                                   'label'    :'$\mu$'},
                                                   
                     'Radiation field strength' : {'shortname':'Rad', 
                                                   'name'     :'Radiation field strength', 
                                                   'dset'     :'RadField', 
                                                   'label'    :'log J/J$_{\mathrm{0}}$'},
                                                   
                     'Total H shielding column' : {'shortname':'NHtot', 
                                                   'name'     :'Total H shielding column', 
                                                   'dset'     :'ShieldingColumn', 
                                                   'label'    :'log N$_{\mathrm{H}}$ [cm$^{-2}$]'},
                                                   
                     'Normalization column'     : {'shortname':'NHnorm',
                                                   'name'     :'Reference column density', 
                                                   'dset'     :'ShieldingColumnNorm', 
                                                   'label'    :'log N$_{\mathrm{H}}$ [cm$^{-2}$]'},

                     'Reference column'         : {'shortname':'NHref',
                                                   'name'     :'Reference column density', 
                                                   'dset'     :'ShieldingColumnRef', 
                                                   'label'    :'log N$_{\mathrm{H}}$ [cm$^{-2}$]'},
                                                   
                     'Internal energy'          : {'shortname':'UfromT', 
                                                   'name'     :'Internal energy', 
                                                   'dset'     :'U_from_T', 
                                                   'label'    :'log U [erg g$^{-1}$]'}
                    }

dict_names_2panels = {'Visual extinction'       : {'shortname':'AV', 
                                                   'name'     :'Visual extinction', 
                                                   'dset1'    :'AVextend', 
                                                   'dset2'    :'AVpoint', 
                                                   'top2D1'   :'AV (extend)',
                                                   'top2D2'   :'AV (point)',                                                   
                                                   'leglabel1':'AV (extend)',
                                                   'leglabel2':'AV (point)',
                                                   'label'    :'log A$_{\mathrm{V}}$ [mag]',
                                                   'label1'   :'log A$_{\mathrm{V}}$ (ext) [mag]', 
                                                   'label2'   :'log A$_{\mathrm{V}}$ (point) [mag]' },
                     
                      'CO fractions'            : {'shortname':'CO', 
                                                   'name'     :'CO fractions',            
                                                   'dset1'    :'COFractionCol', 
                                                   'dset2'    :'COFractionVol', 
                                                   'top2D1'   :'Column density',
                                                   'top2D2'   :'Volume density',
                                                   'leglabel1':'log N$_{\mathrm{CO}}$/N$_{\mathrm{C}}$',
                                                   'leglabel2':'log n$_{\mathrm{CO}}$/n$_{\mathrm{C}}$',                                                   
                                                   'label'    :'log f', 
                                                   'label1'   :'log N$_{\mathrm{CO}}$/N$_{\mathrm{C}}$', 
                                                   'label2'   :'log n$_{\mathrm{CO}}$/n$_{\mathrm{C}}$'}
                     }
 
dict_names_3panels = {'Carbon column densities' : {'shortname':'CCol',
                                                   'name'     :'Carbon column densities',
                                                   'dset'     :'ColumnDensitiesC',
                                                   'label'    :'log N [cm$^{-2}$]',
                                                   'top2D1'   :'C I',
                                                   'top2D2'   :'C II',
                                                   'top2D3'   :'CO',
                                                   'cmin'     :'10',
                                                   'cmax'     :'20',
                                                   'label1'   :'log N$_{\mathrm{CI}}$ [cm$^{-2}$]',
                                                   'label2'   :'log N$_{\mathrm{CII}}$ [cm$^{-2}$]',
                                                   'label3'   :'log N$_{\mathrm{CO}}$ [cm$^{-2}$]'},
                                                   
                      'Hydrogen column densities':{'shortname':'HCol',
                                                   'name'     :'Hydrogen column densities',
                                                   'dset'     :'ColumnDensitiesH',
                                                   'label'    :'log N [cm$^{-2}$]',
                                                   'top2D1'   :'H I',
                                                   'top2D2'   :'H II',
                                                   'top2D3'   :'H2',  
                                                   'cmin'     :'16',
                                                   'cmax'     :'24',                                                   
                                                   'label1'   :'log N$_{\mathrm{HI}}$ [cm$^{-2}$]',
                                                   'label2'   :'log N$_{\mathrm{HII}}$ [cm$^{-2}$]',
                                                   'label3'   :'log N$_{\mathrm{H2}}$ [cm$^{-2}$]'}, 

                      'Hydrogen fractions (col)' :{'shortname':'HfracCol',
                                                   'name'     :'Hydrogen fractions (col)',
                                                   'dset'     :'HydrogenFractionsCol',
                                                   'label'    :'log N$_{\mathrm{x}}$/N$_{\mathrm{H}}$',
                                                   'top2D1'   :'H I',
                                                   'top2D2'   :'H II',
                                                   'top2D3'   :'H2',
                                                   'cmin'     :'-7.8',
                                                   'cmax'     :'0.2',                                                 
                                                   'label1'   :'log N$_{\mathrm{HI}}$/N$_{\mathrm{H}}$',
                                                   'label2'   :'log N$_{\mathrm{HII}}$/N$_{\mathrm{H}}$',
                                                   'label3'   :'log 2 N$_{\mathrm{H2}}$/N$_{\mathrm{H}}$'},

                      'Hydrogen fractions (vol)' :{'shortname':'HfracVol',
                                                   'name'     :'Hydrogen fractions (vol)',
                                                   'dset'     :'HydrogenFractionsVol',
                                                   'label'    :'log n$_{\mathrm{x}}$/n$_{\mathrm{H}}$',
                                                   'top2D1'   :'H I',
                                                   'top2D2'   :'H II',
                                                   'top2D3'   :'H2', 
                                                   'cmin'     :'-7.8',
                                                   'cmax'     :'0.2', 
                                                   'label1'   :'log n$_{\mathrm{HI}}$/n$_{\mathrm{H}}$',
                                                   'label2'   :'log n$_{\mathrm{HII}}$/n$_{\mathrm{H}}$',
                                                   'label3'   :'log 2 n$_{\mathrm{H2}}$/n$_{\mathrm{H}}$'}
                     }
                      
dict_names_cool =    {'Cooling'                  :{'shortname':'Cool',
                                                   'name'     :'Cooling',
                                                   'dset'     :'Cooling'}
                     }                           
                      
dict_names_heat =    {'Heating'                  :{'shortname':'Heat',
                                                   'name'     :'Heating',
                                                   'dset'     :'Heating'}
                     }   
                      
dict_names_ionfrac = {'Ion fractions: H'         :{'shortname':'IonFracH',
                                                   'name'     :'Ion fractions: H',
                                                   'label'    :'log n$_{\mathrm{x}}$/n$_{\mathrm{H}}$',
                                                   'nrxpanels':'2',
                                                   'nrypanels':'1',
                                                   'dset'     :'IonFractionsVol/00hydrogen'},
                      
                      'Ion fractions: He'        :{'shortname':'IonFracHe',
                                                   'name'     :'Ion fractions: He',
                                                   'label'    :'log n$_{\mathrm{x}}$/n$_{\mathrm{He}}$',
                                                   'nrxpanels':'3',
                                                   'nrypanels':'1',                                                   
                                                   'dset'     :'IonFractionsVol/01helium'},
                                                   
                      'Ion fractions: C'         :{'shortname':'IonFracC',
                                                   'name'     :'Ion fractions: C',
                                                   'label'    :'log n$_{\mathrm{x}}$/n$_{\mathrm{C}}$',
                                                   'nrxpanels':'4',
                                                   'nrypanels':'2',                                                  
                                                   'dset'     :'IonFractionsVol/02carbon'},
                                                   
                      'Ion fractions: N'         :{'shortname':'IonFracN',
                                                   'name'     :'Ion fractions: N',
                                                   'label'    :'log n$_{\mathrm{x}}$/n$_{\mathrm{N}}$',
                                                   'nrxpanels':'4',
                                                   'nrypanels':'2',                                                   
                                                   'dset'     :'IonFractionsVol/03nitrogen'},
                                                   
                      'Ion fractions: O'         :{'shortname':'IonFracO',
                                                   'name'     :'Ion fractions: O',
                                                   'label'    :'log n$_{\mathrm{x}}$/n$_{\mathrm{O}}$',
                                                   'nrxpanels':'5',
                                                   'nrypanels':'2',                                                   
                                                   'dset'     :'IonFractionsVol/04oxygen'},
                                                   
                      'Ion fractions: Ne'        :{'shortname':'IonFracNe',
                                                   'name'     :'Ion fractions: Ne',
                                                   'label'    :'log n$_{\mathrm{x}}$/n$_{\mathrm{Ne}}$',
                                                   'nrxpanels':'4',
                                                   'nrypanels':'3',                                                   
                                                   'dset'     :'IonFractionsVol/05neon'},
                                                   
                      'Ion fractions: Mg'        :{'shortname':'IonFracMg',
                                                   'name'     :'Ion fractions: Mg',
                                                   'label'    :'log n$_{\mathrm{x}}$/n$_{\mathrm{Mg}}$',
                                                   'nrxpanels':'5',
                                                   'nrypanels':'3',                                                   
                                                   'dset'     :'IonFractionsVol/06magnesium'},
                                                   
                      'Ion fractions: Si'        :{'shortname':'IonFracSi',
                                                   'name'     :'Ion fractions: Si',
                                                   'label'    :'log n$_{\mathrm{x}}$/n$_{\mathrm{Si}}$',
                                                   'nrxpanels':'5',
                                                   'nrypanels':'3',                                                   
                                                   'dset'     :'IonFractionsVol/07silicon'},
                                                   
                      'Ion fractions: S'         :{'shortname':'IonFracS',
                                                   'name'     :'Ion fractions: S',
                                                   'label'    :'log n$_{\mathrm{x}}$/n$_{\mathrm{S}}$',
                                                   'nrxpanels':'6',
                                                   'nrypanels':'3',                                                   
                                                   'dset'     :'IonFractionsVol/08sulphur'},
                                                   
                      'Ion fractions: Ca'        :{'shortname':'IonFracCa',
                                                   'name'     :'Ion fractions: Ca',
                                                   'label'    :'log n$_{\mathrm{x}}$/n$_{\mathrm{Ca}}$',
                                                   'nrxpanels':'6',
                                                   'nrypanels':'4',                                                   
                                                   'dset'     :'IonFractionsVol/09calcium'},
                                                   
                      'Ion fractions: Fe'        :{'shortname':'IonFracFe',
                                                   'name'     :'Ion fractions: Fe',
                                                   'label'    :'log n$_{\mathrm{x}}$/n$_{\mathrm{Fe}}$',
                                                   'nrxpanels':'6',
                                                   'nrypanels':'5',                                                   
                                                   'dset'     :'IonFractionsVol/10iron'}
                     }                      

dict_names_elec    = {'Electron fractions'       :{'shortname':'ElecFrac',
                                                   'name'     :'Electron fractions',
                                                   'label'    :'log n$_{\mathrm{e}}$/n$_{\mathrm{H}}$',
                                                   'nrxpanels':'5',
                                                   'nrypanels':'3',
                                                   'dset'     :'ElectronFractionsVol'}
                     }

dict_names_depl    = {'Depletion'                :{'shortname':'Depl',
                                                   'name'     :'Depletion',
                                                   'label'    :'log f$_{\mathrm{dust}}$',
                                                   'nrxpanels':'4',
                                                   'nrypanels':'3',
                                                   'dset'     :'Depletion'}
                     }

dict_names_hydext  = {'Hydrogen fractions (ext)' :{'shortname':'HfracVolExt',
                                                   'name'     :'Hydrogen fractions (ext)',
                                                   'label'    :'log n$_{\mathrm{x}}$/n$_{\mathrm{H}}$',
                                                   'nrxpanels':'3',
                                                   'nrypanels':'2',
                                                   'top2D1'   :'H I',
                                                   'top2D2'   :'H II',
                                                   'top2D3'   :'H2',
                                                   'top2D4'   :'H2+',
                                                   'top2D5'   :'H3+',
                                                   'top2D6'   :'H-',     
                                                   'label'    :'log n$_{\mathrm{x}}$/n$_{\mathrm{H}}$',
                                                   'label1'   :'n$_{\mathrm{HI}}$/n$_{\mathrm{H}}$',
                                                   'label2'   :'n$_{\mathrm{HII}}$/n$_{\mathrm{H}}$',
                                                   'label3'   :'2 n$_{\mathrm{H2}}$/n$_{\mathrm{H}}$',
                                                   'label4'   :'2 n$_{\mathrm{H2+}}$/n$_{\mathrm{H}}$',
                                                   'label5'   :'3 n$_{\mathrm{H3+}}$/n$_{\mathrm{H}}$',
                                                   'label6'   :'n$_{\mathrm{H-}}$/n$_{\mathrm{H}}$',                                                   
                                                   'dset'     :'HydrogenFractionsVolExtended'}
                     }

dict_names = merge_all_dicts(dict_names_simple, 
                             dict_names_2panels, 
                             dict_names_3panels, 
                             dict_names_cool,
                             dict_names_heat,
                             dict_names_ionfrac,
                             dict_names_elec,
                             dict_names_depl,
                             dict_names_hydext)
