#import numpy as np
#import corner
#import pyfits
#import math
#import matplotlib.pyplot as plt
#from astropy.io import fits
import scipy.constants
#from astropy.cosmology import FlatLambdaCDM
#import astropy.units as u
#cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)
import cosmolopy.distance as cd
cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 'omega_k_0':0.0, 'h':0.72}
#import random



#marginalized_array = np.load('./marginalized/samples_values_all_Bootes1_Bootes1_04.npy')  #=============================================HERE

#print type(marginalized_array)
#print np.shape(marginalized_array)

#print len(marginalized_array[:,0])





list_sources_random = np.random.randint(low=0, high=len(aless_array[:,0]), size=200)  #====================HERE: change random sample number (size) (same for all templates, so just use len(aless_array).
#print list_sources_random


z_want_aless=[]
a_want_aless=[]

for j in range(0, len(list_sources_random)):
    z_want_aless.append(aless_array[list_sources_random[j],0])
    a_want_aless.append(aless_array[list_sources_random[j],1])


z_want_arp220=[]
a_want_arp220=[]

for j in range(0, len(list_sources_random)):
    z_want_arp220.append(arp220_array[list_sources_random[j],0])
    a_want_arp220.append(arp220_array[list_sources_random[j],1])


z_want_hfls3=[]
a_want_hfls3=[]

for j in range(0, len(list_sources_random)):
    z_want_hfls3.append(hfls3_array[list_sources_random[j],0])
    a_want_hfls3.append(hfls3_array[list_sources_random[j],1])


z_want_eyelash=[]
a_want_eyelash=[]

for j in range(0, len(list_sources_random)):
    z_want_eyelash.append(eyelash_array[list_sources_random[j],0])
    a_want_eyelash.append(eyelash_array[list_sources_random[j],1])



#print z_want
#print a_want
#print '33333333333333333333'



#----------------template imported before---------------------

#-------------------------------------------------------------



L_want_aless=[]
SFR_want_aless=[]

L_want_arp220=[]
SFR_want_arp220=[]

L_want_hfls3=[]
SFR_want_hfls3=[]

L_want_eyelash=[]
SFR_want_eyelash=[]



marginalized_array_random_aless = []
marginalized_array_random_arp220 = []
marginalized_array_random_hfls3 = []
marginalized_array_random_eyelash = []


#FIR_wavelength = np.logspace(np.log10(8.0), 3.0, num=100, endpoint=True, base=10.0)

for i in range(0, len(z_want_aless)):     #loop over all samples (randomly selected down: should be same for each template)
#for i in range(0, 2):
    print i

    #---------------I. aless-------------------------------
    wavelength_micron_redshifted_aless=[]
    flux_mjy_redshifted_aless=[]
    for j in range(0, len(wavelength_micron_aless)):

        wavelength_micron_redshifted_aless.append(wavelength_micron_aless[j]*(1.0+z_want_aless[i]))
        flux_mjy_redshifted_aless.append(flux_mjy_aless[j]*pow(10.0,a_want_aless[i]))

    #print z_want[i]  #z of this sample
    #print a_want[i]  #a of this sample

    wavelength_micron_redshifted_cut_aless=[]
    flux_mjy_redshifted_cut_aless=[]
    for k in range(0, len(wavelength_micron_redshifted_aless)):
        if wavelength_micron_redshifted_aless[k] >= 8.0*(1.0+z_want_aless[i]) and wavelength_micron_redshifted_aless[k] <= 1000.0*(1.0+z_want_aless[i]):
            wavelength_micron_redshifted_cut_aless.append(wavelength_micron_redshifted_aless[k])
            flux_mjy_redshifted_cut_aless.append(flux_mjy_redshifted_aless[k])

    #print wavelength_micron_redshifted_cut[0:50]
    #print flux_mjy_redshifted_cut[0:50]
    

    #------Unit conversion------------(micron-->m-->Hz-->reverse; mJy-->Jy-->W*m-2*Hz-1-->reverse-->*4pi*luminosity_distance)
    wavelength_sed_aless=[]
    flux_sed_aless=[]
    for m in range(0, len(wavelength_micron_redshifted_cut_aless)):
        wavelength_sed_aless.append(scipy.constants.c/(wavelength_micron_redshifted_cut_aless[m]*1e-6))  #now in Hz
        #flux_sed.append(1e-26*flux_mjy_redshifted_cut[m]*1e-3*4.0*math.pi*pow(1e6*3.0857e16*cosmo.luminosity_distance(z_want[i]).value, 2.0))  #now in W/Hz
        flux_sed_aless.append(1e-26*flux_mjy_redshifted_cut_aless[m]*1e-3*4.0*math.pi*pow(1e6*3.0857e16*cd.luminosity_distance(z_want_aless[i], **cosmo), 2.0))  #now in W/Hz

    wavelength_sed_aless=wavelength_sed_aless[::-1]  #reverse so get increasing in Hz  
    flux_sed_aless=flux_sed_aless[::-1]

    #print wavelength_sed[0:50]
    #print flux_sed[0:50]

    L_area_aless = np.trapz(flux_sed_aless, wavelength_sed_aless)
    L_area_solar_aless = L_area_aless/3.828e26  #in unit of solar luminosity

    L_area_solar_aless = np.log10(L_area_solar_aless)  #make it log to plot in corner diagram
 
    #print L_area_solar

    L_want_aless.append(L_area_solar_aless)

    #SFR_loop_aless=2.8e-44*L_area_aless*1e7  #in unit of M_solar/yr (constant from https://ned.ipac.caltech.edu/level5/Sept12/Calzetti/Calzetti1_2.html, L_area need to be in erg/s)
    SFR_loop_aless=3.89e-44*L_area_aless*1e7  #in unit of M_solar/yr (constant from Kennedicutt & Evans 2012, L_area need to be in erg/s)
    SFR_loop_aless=np.log10(SFR_loop_aless)  #make it log to plot in corner diagram

    #print SFR_loop
    #print 'HHHHHHHHHHHHHHHHHHHHHHHHH'
    
    SFR_want_aless.append(SFR_loop_aless)





    #---------------II. arp220-------------------------------

    wavelength_micron_redshifted_arp220=[]
    flux_mjy_redshifted_arp220=[]
    for j in range(0, len(wavelength_micron_arp220)):

        wavelength_micron_redshifted_arp220.append(wavelength_micron_arp220[j]*(1.0+z_want_arp220[i]))
        flux_mjy_redshifted_arp220.append(flux_mjy_arp220[j]*pow(10.0,a_want_arp220[i]))

    #print z_want[i]  #z of this sample
    #print a_want[i]  #a of this sample

    wavelength_micron_redshifted_cut_arp220=[]
    flux_mjy_redshifted_cut_arp220=[]
    for k in range(0, len(wavelength_micron_redshifted_arp220)):
        if wavelength_micron_redshifted_arp220[k] >= 8.0*(1.0+z_want_arp220[i]) and wavelength_micron_redshifted_arp220[k] <= 1000.0*(1.0+z_want_arp220[i]):
            wavelength_micron_redshifted_cut_arp220.append(wavelength_micron_redshifted_arp220[k])
            flux_mjy_redshifted_cut_arp220.append(flux_mjy_redshifted_arp220[k])

    #print wavelength_micron_redshifted_cut[0:50]
    #print flux_mjy_redshifted_cut[0:50]
    

    #------Unit conversion------------(micron-->m-->Hz-->reverse; mJy-->Jy-->W*m-2*Hz-1-->reverse-->*4pi*luminosity_distance)
    wavelength_sed_arp220=[]
    flux_sed_arp220=[]
    for m in range(0, len(wavelength_micron_redshifted_cut_arp220)):
        wavelength_sed_arp220.append(scipy.constants.c/(wavelength_micron_redshifted_cut_arp220[m]*1e-6))  #now in Hz
        #flux_sed.append(1e-26*flux_mjy_redshifted_cut[m]*1e-3*4.0*math.pi*pow(1e6*3.0857e16*cosmo.luminosity_distance(z_want[i]).value, 2.0))  #now in W/Hz
        flux_sed_arp220.append(1e-26*flux_mjy_redshifted_cut_arp220[m]*1e-3*4.0*math.pi*pow(1e6*3.0857e16*cd.luminosity_distance(z_want_arp220[i], **cosmo), 2.0))  #now in W/Hz

    wavelength_sed_arp220=wavelength_sed_arp220[::-1]  #reverse so get increasing in Hz  
    flux_sed_arp220=flux_sed_arp220[::-1]

    #print wavelength_sed[0:50]
    #print flux_sed[0:50]

    L_area_arp220 = np.trapz(flux_sed_arp220, wavelength_sed_arp220)
    L_area_solar_arp220 = L_area_arp220/3.828e26  #in unit of solar luminosity

    L_area_solar_arp220 = np.log10(L_area_solar_arp220)  #make it log to plot in corner diagram
 
    #print L_area_solar

    L_want_arp220.append(L_area_solar_arp220)

    #SFR_loop_arp220=2.8e-44*L_area_arp220*1e7  #in unit of M_solar/yr (constant from https://ned.ipac.caltech.edu/level5/Sept12/Calzetti/Calzetti1_2.html, L_area need to be in erg/s)
    SFR_loop_arp220=3.89e-44*L_area_arp220*1e7  #in unit of M_solar/yr (constant from Kennedicutt & Evans 2012, L_area need to be in erg/s)
    SFR_loop_arp220=np.log10(SFR_loop_arp220)  #make it log to plot in corner diagram

    #print SFR_loop
    #print 'HHHHHHHHHHHHHHHHHHHHHHHHH'
    
    SFR_want_arp220.append(SFR_loop_arp220)

    #---------------III. hfls3-------------------------------


    wavelength_micron_redshifted_hfls3=[]
    flux_mjy_redshifted_hfls3=[]
    for j in range(0, len(wavelength_micron_hfls3)):

        wavelength_micron_redshifted_hfls3.append(wavelength_micron_hfls3[j]*(1.0+z_want_hfls3[i]))
        flux_mjy_redshifted_hfls3.append(flux_mjy_hfls3[j]*pow(10.0,a_want_hfls3[i]))

    #print z_want[i]  #z of this sample
    #print a_want[i]  #a of this sample

    wavelength_micron_redshifted_cut_hfls3=[]
    flux_mjy_redshifted_cut_hfls3=[]
    for k in range(0, len(wavelength_micron_redshifted_hfls3)):
        if wavelength_micron_redshifted_hfls3[k] >= 8.0*(1.0+z_want_hfls3[i]) and wavelength_micron_redshifted_hfls3[k] <= 1000.0*(1.0+z_want_hfls3[i]):
            wavelength_micron_redshifted_cut_hfls3.append(wavelength_micron_redshifted_hfls3[k])
            flux_mjy_redshifted_cut_hfls3.append(flux_mjy_redshifted_hfls3[k])

    #print wavelength_micron_redshifted_cut[0:50]
    #print flux_mjy_redshifted_cut[0:50]
    

    #------Unit conversion------------(micron-->m-->Hz-->reverse; mJy-->Jy-->W*m-2*Hz-1-->reverse-->*4pi*luminosity_distance)
    wavelength_sed_hfls3=[]
    flux_sed_hfls3=[]
    for m in range(0, len(wavelength_micron_redshifted_cut_hfls3)):
        wavelength_sed_hfls3.append(scipy.constants.c/(wavelength_micron_redshifted_cut_hfls3[m]*1e-6))  #now in Hz
        #flux_sed.append(1e-26*flux_mjy_redshifted_cut[m]*1e-3*4.0*math.pi*pow(1e6*3.0857e16*cosmo.luminosity_distance(z_want[i]).value, 2.0))  #now in W/Hz
        flux_sed_hfls3.append(1e-26*flux_mjy_redshifted_cut_hfls3[m]*1e-3*4.0*math.pi*pow(1e6*3.0857e16*cd.luminosity_distance(z_want_hfls3[i], **cosmo), 2.0))  #now in W/Hz

    wavelength_sed_hfls3=wavelength_sed_hfls3[::-1]  #reverse so get increasing in Hz  
    flux_sed_hfls3=flux_sed_hfls3[::-1]

    #print wavelength_sed[0:50]
    #print flux_sed[0:50]

    L_area_hfls3 = np.trapz(flux_sed_hfls3, wavelength_sed_hfls3)
    L_area_solar_hfls3 = L_area_hfls3/3.828e26  #in unit of solar luminosity

    L_area_solar_hfls3 = np.log10(L_area_solar_hfls3)  #make it log to plot in corner diagram
 
    #print L_area_solar

    L_want_hfls3.append(L_area_solar_hfls3)

    #SFR_loop_hfls3=2.8e-44*L_area_hfls3*1e7  #in unit of M_solar/yr (constant from https://ned.ipac.caltech.edu/level5/Sept12/Calzetti/Calzetti1_2.html, L_area need to be in erg/s)
    SFR_loop_hfls3=3.89e-44*L_area_hfls3*1e7  #in unit of M_solar/yr (constant from Kennedicutt & Evans 2012, L_area need to be in erg/s)
    SFR_loop_hfls3=np.log10(SFR_loop_hfls3)  #make it log to plot in corner diagram

    #print SFR_loop
    #print 'HHHHHHHHHHHHHHHHHHHHHHHHH'
    
    SFR_want_hfls3.append(SFR_loop_hfls3)

    #---------------IV. eyelash-------------------------------


    wavelength_micron_redshifted_eyelash=[]
    flux_mjy_redshifted_eyelash=[]
    for j in range(0, len(wavelength_micron_eyelash)):

        wavelength_micron_redshifted_eyelash.append(wavelength_micron_eyelash[j]*(1.0+z_want_eyelash[i]))
        flux_mjy_redshifted_eyelash.append(flux_mjy_eyelash[j]*pow(10.0,a_want_eyelash[i]))

    #print z_want[i]  #z of this sample
    #print a_want[i]  #a of this sample

    wavelength_micron_redshifted_cut_eyelash=[]
    flux_mjy_redshifted_cut_eyelash=[]
    for k in range(0, len(wavelength_micron_redshifted_eyelash)):
        if wavelength_micron_redshifted_eyelash[k] >= 8.0*(1.0+z_want_eyelash[i]) and wavelength_micron_redshifted_eyelash[k] <= 1000.0*(1.0+z_want_eyelash[i]):
            wavelength_micron_redshifted_cut_eyelash.append(wavelength_micron_redshifted_eyelash[k])
            flux_mjy_redshifted_cut_eyelash.append(flux_mjy_redshifted_eyelash[k])

    #print wavelength_micron_redshifted_cut[0:50]
    #print flux_mjy_redshifted_cut[0:50]
    

    #------Unit conversion------------(micron-->m-->Hz-->reverse; mJy-->Jy-->W*m-2*Hz-1-->reverse-->*4pi*luminosity_distance)
    wavelength_sed_eyelash=[]
    flux_sed_eyelash=[]
    for m in range(0, len(wavelength_micron_redshifted_cut_eyelash)):
        wavelength_sed_eyelash.append(scipy.constants.c/(wavelength_micron_redshifted_cut_eyelash[m]*1e-6))  #now in Hz
        #flux_sed.append(1e-26*flux_mjy_redshifted_cut[m]*1e-3*4.0*math.pi*pow(1e6*3.0857e16*cosmo.luminosity_distance(z_want[i]).value, 2.0))  #now in W/Hz
        flux_sed_eyelash.append(1e-26*flux_mjy_redshifted_cut_eyelash[m]*1e-3*4.0*math.pi*pow(1e6*3.0857e16*cd.luminosity_distance(z_want_eyelash[i], **cosmo), 2.0))  #now in W/Hz

    wavelength_sed_eyelash=wavelength_sed_eyelash[::-1]  #reverse so get increasing in Hz  
    flux_sed_eyelash=flux_sed_eyelash[::-1]

    #print wavelength_sed[0:50]
    #print flux_sed[0:50]

    L_area_eyelash = np.trapz(flux_sed_eyelash, wavelength_sed_eyelash)
    L_area_solar_eyelash = L_area_eyelash/3.828e26  #in unit of solar luminosity

    L_area_solar_eyelash = np.log10(L_area_solar_eyelash)  #make it log to plot in corner diagram
 
    #print L_area_solar

    L_want_eyelash.append(L_area_solar_eyelash)

    #SFR_loop_eyelash=2.8e-44*L_area_eyelash*1e7  #in unit of M_solar/yr (constant from https://ned.ipac.caltech.edu/level5/Sept12/Calzetti/Calzetti1_2.html, L_area need to be in erg/s)
    SFR_loop_eyelash=3.89e-44*L_area_eyelash*1e7  #in unit of M_solar/yr (constant from Kennedicutt & Evans 2012, L_area need to be in erg/s)
    SFR_loop_eyelash=np.log10(SFR_loop_eyelash)  #make it log to plot in corner diagram

    #print SFR_loop
    #print 'HHHHHHHHHHHHHHHHHHHHHHHHH'
    
    SFR_want_eyelash.append(SFR_loop_eyelash)


    #----------------------------------------------------------------












marginalized_array_random_aless = np.c_[np.array(z_want_aless),np.array(a_want_aless), np.array(L_want_aless), np.array(SFR_want_aless)]
marginalized_array_random_arp220 = np.c_[np.array(z_want_arp220),np.array(a_want_arp220), np.array(L_want_arp220), np.array(SFR_want_arp220)]
marginalized_array_random_hfls3 = np.c_[np.array(z_want_hfls3),np.array(a_want_hfls3), np.array(L_want_hfls3), np.array(SFR_want_hfls3)]
marginalized_array_random_eyelash = np.c_[np.array(z_want_eyelash),np.array(a_want_eyelash), np.array(L_want_eyelash), np.array(SFR_want_eyelash)]


#print marginalized_array_random_aless
print type(marginalized_array_random_aless)
print np.shape(marginalized_array_random_aless)
print len(marginalized_array_random_aless[:,0])
print 'xxxxxxxxxxxxxxxxxx'


lower_limit_L_aless = corner.quantile(marginalized_array_random_aless[:,2], q=[0.16], weights=None)
fitted_median_L_aless = corner.quantile(marginalized_array_random_aless[:,2], q=[0.5], weights=None)
higher_limit_L_aless = corner.quantile(marginalized_array_random_aless[:,2], q=[0.84], weights=None)
#print lower_limit_L_aless[0]
#print fitted_median_L_aless[0]
#print higher_limit_L_aless[0]

L_distribution_all_aless.append(fitted_median_L_aless[0])
L_distribution_all_uperr_aless.append(higher_limit_L_aless[0]-fitted_median_L_aless[0])
L_distribution_all_lowerr_aless.append(fitted_median_L_aless[0]-lower_limit_L_aless[0])




lower_limit_L_arp220 = corner.quantile(marginalized_array_random_arp220[:,2], q=[0.16], weights=None)
fitted_median_L_arp220= corner.quantile(marginalized_array_random_arp220[:,2], q=[0.5], weights=None)
higher_limit_L_arp220 = corner.quantile(marginalized_array_random_arp220[:,2], q=[0.84], weights=None)
L_distribution_all_arp220.append(fitted_median_L_arp220[0])
L_distribution_all_uperr_arp220.append(higher_limit_L_arp220[0]-fitted_median_L_arp220[0])
L_distribution_all_lowerr_arp220.append(fitted_median_L_arp220[0]-lower_limit_L_arp220[0])

lower_limit_L_hfls3 = corner.quantile(marginalized_array_random_hfls3[:,2], q=[0.16], weights=None)
fitted_median_L_hfls3 = corner.quantile(marginalized_array_random_hfls3[:,2], q=[0.5], weights=None)
higher_limit_L_hfls3 = corner.quantile(marginalized_array_random_hfls3[:,2], q=[0.84], weights=None)
L_distribution_all_hfls3.append(fitted_median_L_hfls3[0])
L_distribution_all_uperr_hfls3.append(higher_limit_L_hfls3[0]-fitted_median_L_hfls3[0])
L_distribution_all_lowerr_hfls3.append(fitted_median_L_hfls3[0]-lower_limit_L_hfls3[0])

lower_limit_L_eyelash = corner.quantile(marginalized_array_random_eyelash[:,2], q=[0.16], weights=None)
fitted_median_L_eyelash = corner.quantile(marginalized_array_random_eyelash[:,2], q=[0.5], weights=None)
higher_limit_L_eyelash = corner.quantile(marginalized_array_random_eyelash[:,2], q=[0.84], weights=None)
L_distribution_all_eyelash.append(fitted_median_L_eyelash[0])
L_distribution_all_uperr_eyelash.append(higher_limit_L_eyelash[0]-fitted_median_L_eyelash[0])
L_distribution_all_lowerr_eyelash.append(fitted_median_L_eyelash[0]-lower_limit_L_eyelash[0])








#print np.shape(L_want)

#np.save('./marginalized/samples_luminosity_all.npy', L_want)
#print L_want

#fig2 = corner.corner(marginalized_array_random_aless, labels=["$z$", "$log(a)$", r"$log(L(L_{\odot}))$", r"$log(SFR(M_{\odot}/yr))$"],
#                    quantiles=[0.16, 0.5, 0.84],
#                    show_titles=True, title_kwargs={"fontsize": 12})
        
#fig2.savefig("./marginalized/test_"+cluster_list[l]+"_"+source_name[x]+".png")   #=========================================================HERE
#plt.close(fig)


#print 'SSSSSSSSSSSSSSSSSSSSSSSTOPPPPPPPPPPPPPPPPPPPPPPPPP'

