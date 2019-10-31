import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm
from astropy.io import fits
import astropy.wcs as wcs
from scipy.interpolate import interp1d
import pyfits
import os
import subprocess
from astropy.table import Table
import emcee
import corner
import random


cluster_list=['PCL1002', 'S2CLS_random']
#cluster_list=['NGP1', 'NGP2', 'NGP3', 'NGP4', 'NGP5', 'NGP6', 'NGP7', 'NGP8', 'NGP9', 'PCL1002', 'S2CLS_random']
#cluster_list=['NGP2', 'NGP3', 'NGP4', 'NGP5', 'NGP6', 'NGP7', 'NGP8', 'NGP9']
#cluster_list=['Bootes1', 'EGS', 'Lockman', 'G12', 'NGP1', 'NGP2', 'NGP3', 'NGP4', 'NGP5', 'NGP6', 'NGP7', 'NGP8', 'NGP9', 'PCL1002', 'S2CLS_random']
cluster_list=['totalSFR']

#cluster_list=['PCL1002', 'S2CLS_random']



#=====================import template(s) redshift it, interpolate the model values at the correspoinding bands, and create parameter space (z and L --> delx and dely)==================

#template_name='aless' #=========================================================================HERE: change template name
template_name=['aless', 'arp220', 'hfls3', 'eyelash'] #=========================================================================HERE: change template name
#template_name=['aless'] #=========================================================================HERE: change template name

for w in range(0, len(template_name)):

    print template_name[w]
    if template_name[w]=='aless':
        template_file = open("./aless_average_seds.dat", "r")    
    if template_name[w]=='arp220':
        template_file = open("./arp220_template_donley2007_apj_660_167.txt", "r")    
    if template_name[w]=='hfls3':
        template_file = open("./HFLS3.txt", "r")   
    if template_name[w]=='eyelash':
        template_file = open("./cosmic_eyelash_sed_josh.txt", "r")  


    lines = template_file.readlines()
    template_file.close()

    response_curve=[]

    for i in range(0, len(lines)):
        separated_lines=lines[i].split() 
        response_curve.append(separated_lines)


    response_curve = np.array(response_curve)
    wavelength_micron=response_curve[:,0]
    flux_mjy=response_curve[:,1]


    wavelength_micron = np.array(wavelength_micron)
    flux_mjy = np.array(flux_mjy)
    wavelength_micron=np.array([float(i) for i in wavelength_micron])  #(wavelength in angstrom)
    flux_mjy=np.array([float(i) for i in flux_mjy])  #(wavelength in angstrom)


    #+++++++++++++++++++++For HFLS3+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=====================================: for HFLS3 template
    if template_name[w]=='hfls3':
        wavelength_micron=wavelength_micron/7.34  #redshift back to restframe (z=6.34)
        wavelength_micron=wavelength_micron*1.0e6  #convert from m to micron
        flux_mjy=flux_mjy*1000.0  #convert from Jy to mJy


    #+++++++++++++++++++++For Cosmic Eyelash+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=====================================: for Cosmic Eyelash template (no need if using Josh's template)
    #if template_name=='eyelash':

    #    flux_mjy_loop=[]
    #    for j in range(0, len(wavelength_micron)):

    #        flux_mjy_loop.append(1.0e-14*1000.0*(1.0e26)*flux_mjy[j]*(wavelength_micron[j]*1.0e-6)/2.99792458e8)


    #    flux_mjy=flux_mjy_loop

    #    flux_mjy = np.array(flux_mjy)
    #    flux_mjy=np.array([float(i) for i in flux_mjy])  #(wavelength in angstrom)

    #plt.plot(wavelength_micron, flux_mjy)
    #plt.xscale('log',nonposy='clip')
    #plt.yscale('log',nonposy='clip')
    #plt.show()



    for l in range(0, len(cluster_list)):

        print cluster_list[l]
        #execfile("./photoz_mcmc_181114.py")  #==================================================================================HERE
        #execfile("./photoz_mcmc_190410.py")  #==================================================================================HERE
        #execfile("./photoz_mcmc_190712.py")  #==================================================================================HERE
        #execfile("./photoz_mcmc_181116_field_casey.py")  #==================================================================================HERE
        execfile("./photoz_mcmc_totalSFR.py")  #==================================================================================HERE


os.system('spd-say "your program has finished"')

