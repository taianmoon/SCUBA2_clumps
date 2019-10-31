#import numpy as np
#import corner
#import pyfits
#import math
#import matplotlib.pyplot as plt
##from astropy.io import fits
#import scipy.constants
##from astropy.cosmology import FlatLambdaCDM
##import astropy.units as u
##cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)
#import cosmolopy.distance as cd
#cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 'omega_k_0':0.0, 'h':0.72}
#




#z_want_arp220 = 4.0  #====================================================Redshift for the SED

##--------------import template first -----------------------------
#
#template_file_arp220 = open("./arp220_template_donley2007_apj_660_167.txt", "r") 
#
#lines_arp220 = template_file_arp220.readlines()
#template_file_arp220.close()
#
#response_curve_arp220=[]
#
#for i in range(0, len(lines_arp220)):
#    separated_lines_arp220=lines_arp220[i].split() 
#    response_curve_arp220.append(separated_lines_arp220)
#
#
#response_curve_arp220 = np.array(response_curve_arp220)
#wavelength_micron_arp220=response_curve_arp220[:,0]
#flux_mjy_arp220=response_curve_arp220[:,1]
#
#
#wavelength_micron_arp220 = np.array(wavelength_micron_arp220)
#flux_mjy_arp220 = np.array(flux_mjy_arp220)
#wavelength_micron_arp220=np.array([float(i) for i in wavelength_micron_arp220])  #(wavelength in angstrom)
#flux_mjy_arp220=np.array([float(i) for i in flux_mjy_arp220])  #(wavelength in angstrom)
#
##print wavelength_micron_arp220
##print flux_mjy_arp220
#
#
#

#---------------Estimate L_IR-------------------------------

wavelength_micron_redshifted_arp220=[]
flux_mjy_redshifted_arp220=[]
for j in range(0, len(wavelength_micron_arp220)):

    wavelength_micron_redshifted_arp220.append(wavelength_micron_arp220[j]*(1.0+z_want_arp220[l]))
    flux_mjy_redshifted_arp220.append(flux_mjy_arp220[j])


wavelength_micron_redshifted_arp220 = np.array(wavelength_micron_redshifted_arp220)
wavelength_micron_redshifted_arp220=np.array([float(i) for i in wavelength_micron_redshifted_arp220])  #(wavelength in angstrom)
flux_mjy_redshifted_arp220 = np.array(flux_mjy_redshifted_arp220)
flux_mjy_redshifted_arp220=np.array([float(i) for i in flux_mjy_redshifted_arp220])  #(wavelength in angstrom)

#print wavelength_micron_redshifted_arp220
#print flux_mjy_redshifted_arp220
#print 'xxxxxxxxxxxxxxxxxxxxxxxxxxxx'
#print z_want[i]  #z of this sample
#print a_want[i]  #a of this sample

#------------------normalise to fixed flux density (8 mJy at observed 850 micron)

mult_factor = 8.0/np.interp(850.0, wavelength_micron_redshifted_arp220, flux_mjy_redshifted_arp220)
#print mult_factor
#print 'ffffffffffffff'

flux_mjy_redshifted_arp220 = flux_mjy_redshifted_arp220*mult_factor

#-------------------------------------------------------------------------------


wavelength_micron_redshifted_cut_arp220=[]
flux_mjy_redshifted_cut_arp220=[]
for k in range(0, len(wavelength_micron_redshifted_arp220)):
    if wavelength_micron_redshifted_arp220[k] >= 8.0*(1.0+z_want_arp220[l]) and wavelength_micron_redshifted_arp220[k] <= 1000.0*(1.0+z_want_arp220[l]):
        wavelength_micron_redshifted_cut_arp220.append(wavelength_micron_redshifted_arp220[k])
        flux_mjy_redshifted_cut_arp220.append(flux_mjy_redshifted_arp220[k])

#print wavelength_micron_redshifted_cut_arp220
#print flux_mjy_redshifted_cut_arp220
#print np.interp(850.0, wavelength_micron_redshifted_cut_arp220, flux_mjy_redshifted_cut_arp220)
#print 'cccccccccccccccccccccc'    

#plt.plot(wavelength_micron_redshifted_cut_arp220, flux_mjy_redshifted_cut_arp220)
#plt.xscale('log')
#plt.yscale('log')
#plt.show()

#------Unit conversion------------(micron-->m-->Hz-->reverse; mJy-->Jy-->W*m-2*Hz-1-->reverse-->*4pi*luminosity_distance)
wavelength_sed_arp220=[]
flux_sed_arp220=[]
for m in range(0, len(wavelength_micron_redshifted_cut_arp220)):
    wavelength_sed_arp220.append(scipy.constants.c/(wavelength_micron_redshifted_cut_arp220[m]*1e-6))  #now in Hz
    #flux_sed.append(1e-26*flux_mjy_redshifted_cut[m]*1e-3*4.0*math.pi*pow(1e6*3.0857e16*cosmo.luminosity_distance(z_want[i]).value, 2.0))  #now in W/Hz
    flux_sed_arp220.append(1e-26*flux_mjy_redshifted_cut_arp220[m]*1e-3*4.0*math.pi*pow(1e6*3.0857e16*cd.luminosity_distance(z_want_arp220[l], **cosmo), 2.0))  #now in W/Hz

wavelength_sed_arp220=wavelength_sed_arp220[::-1]  #reverse so get increasing in Hz  
flux_sed_arp220=flux_sed_arp220[::-1]

#print wavelength_sed[0:50]
#print flux_sed[0:50]

L_area_arp220 = np.trapz(flux_sed_arp220, wavelength_sed_arp220)
L_area_solar_arp220 = L_area_arp220/3.828e26  #in unit of solar luminosity

L_area_solar_arp220 = np.log10(L_area_solar_arp220)  #make it log to plot in corner diagram

#print 'qqqqqqqqqqqqqqqqqqqqq' 
print L_area_solar_arp220
L_area_solar_arp220_loop.append(L_area_solar_arp220)

#L_want_arp220.append(L_area_solar_arp220)

#SFR_loop_arp220=2.8e-44*L_area_arp220*1e7  #in unit of M_solar/yr (constant from https://ned.ipac.caltech.edu/level5/Sept12/Calzetti/Calzetti1_2.html, L_area need to be in erg/s)
#SFR_loop_arp220=3.89e-44*L_area_arp220*1e7  #in unit of M_solar/yr (constant from Kennedicutt & Evans 2012, L_area need to be in erg/s)
#SFR_loop_arp220=np.log10(SFR_loop_arp220)  #make it log to plot in corner diagram

#print SFR_loop
#print 'HHHHHHHHHHHHHHHHHHHHHHHHH'
    
#SFR_want_arp220.append(SFR_loop_arp220)


