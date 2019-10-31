import numpy as np
import corner
import pyfits
import math
import matplotlib.pyplot as plt
#from astropy.io import fits
import scipy.constants
#from astropy.cosmology import FlatLambdaCDM
#import astropy.units as u
#cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)
import cosmolopy.distance as cd
cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 'omega_k_0':0.0, 'h':0.72}





z_want_arp220 = np.linspace(0.0, 6.0, num=50)  #====================================================Redshift for the SED
hfls3_switch = 1


#--------------import template first -----------------------------

#template_file_arp220 = open("./arp220_template_donley2007_apj_660_167.txt", "r") 
#template_file_arp220 = open("./aless_average_seds.dat", "r") 
#template_file_arp220 = open("./cosmic_eyelash_sed_josh.txt", "r") 
template_file_arp220 = open("./HFLS3.txt", "r") 

lines_arp220 = template_file_arp220.readlines()
template_file_arp220.close()

response_curve_arp220=[]

for i in range(0, len(lines_arp220)):
    separated_lines_arp220=lines_arp220[i].split() 
    response_curve_arp220.append(separated_lines_arp220)


response_curve_arp220 = np.array(response_curve_arp220)
wavelength_micron_arp220=response_curve_arp220[:,0]
flux_mjy_arp220=response_curve_arp220[:,1]


wavelength_micron_arp220 = np.array(wavelength_micron_arp220)
flux_mjy_arp220 = np.array(flux_mjy_arp220)
wavelength_micron_arp220=np.array([float(i) for i in wavelength_micron_arp220])  #(wavelength in angstrom)
flux_mjy_arp220=np.array([float(i) for i in flux_mjy_arp220])  #(wavelength in angstrom)

if hfls3_switch == 1:
    wavelength_micron_arp220=wavelength_micron_arp220/7.34  #redshift back to restframe (z=6.34)
    wavelength_micron_arp220=wavelength_micron_arp220*1.0e6  #convert from m to micron
    flux_mjy_arp220=flux_mjy_arp220*1000.0  #convert from Jy to mJy





L_area_solar_arp220_loop=[]

for l in range(0, len(z_want_arp220)):

    execfile("./sed_test_190925.py")


plt.scatter([2.86, 3.35], [12.85, 13.09], color='r', s=20.0)
plt.errorbar([2.86, 3.35], [12.85, 13.09], xerr=[0.96, 1.09], yerr=[0.22, 0.23], color='r')
plt.plot(z_want_arp220, L_area_solar_arp220_loop, linewidth=1.5, color='blue')
plt.xlabel('z')
plt.ylabel('log10(L_IR)')
plt.rc('font', size=20)
plt.show()


mock_catalogue_z=np.column_stack((z_want_arp220, L_area_solar_arp220_loop))
np.savetxt('./seds_results/HFLS3_LIR.cat', mock_catalogue_z, delimiter=' ') #================================================================================HERE




