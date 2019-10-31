import numpy as np
import corner
import pyfits
import math
import matplotlib.pyplot as plt
#from astropy.io import fits


cluster_list=['totalSFR']
#cluster_list=['Bootes1', 'EGS', 'Lockman', 'G12', 'NGP1', 'NGP2', 'NGP3', 'NGP4', 'NGP5', 'NGP6', 'NGP7', 'NGP8', 'NGP9', 'PCL1002', 'S2CLS_random']
template_name=['aless', 'arp220', 'hfls3', 'eyelash'] 

#=================================I. Marginalize over the four templates====================================================================================

z_distribution_all=[]
z_distribution_all_uperr=[]
z_distribution_all_lowerr=[]

a_distribution_all=[]
a_distribution_all_uperr=[]
a_distribution_all_lowerr=[]

L_distribution_all=[]
L_distribution_all_uperr=[]
L_distribution_all_lowerr=[]

SFR_distribution_all=[]
SFR_distribution_all_uperr=[]
SFR_distribution_all_lowerr=[]

map_size_all=0.0


#--------------import template first -----------------------------
template_file_aless = open("../aless_average_seds.dat", "r") 

lines_aless = template_file_aless.readlines()
template_file_aless.close()

response_curve_aless=[]

for i in range(0, len(lines_aless)):
    separated_lines_aless=lines_aless[i].split() 
    response_curve_aless.append(separated_lines_aless)


response_curve_aless = np.array(response_curve_aless)
wavelength_micron_aless=response_curve_aless[:,0]
flux_mjy_aless=response_curve_aless[:,1]


wavelength_micron_aless = np.array(wavelength_micron_aless)
flux_mjy_aless = np.array(flux_mjy_aless)
wavelength_micron_aless=np.array([float(i) for i in wavelength_micron_aless])  #(wavelength in angstrom)
flux_mjy_aless=np.array([float(i) for i in flux_mjy_aless])  #(wavelength in angstrom)




template_file_arp220 = open("../arp220_template_donley2007_apj_660_167.txt", "r") 

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





template_file_hfls3 = open("../HFLS3.txt", "r") 

lines_hfls3 = template_file_hfls3.readlines()
template_file_hfls3.close()

response_curve_hfls3=[]

for i in range(0, len(lines_hfls3)):
    separated_lines_hfls3=lines_hfls3[i].split() 
    response_curve_hfls3.append(separated_lines_hfls3)


response_curve_hfls3 = np.array(response_curve_hfls3)
wavelength_micron_hfls3=response_curve_hfls3[:,0]
flux_mjy_hfls3=response_curve_hfls3[:,1]


wavelength_micron_hfls3 = np.array(wavelength_micron_hfls3)
flux_mjy_hfls3 = np.array(flux_mjy_hfls3)
wavelength_micron_hfls3=np.array([float(i) for i in wavelength_micron_hfls3])  #(wavelength in angstrom)
flux_mjy_hfls3=np.array([float(i) for i in flux_mjy_hfls3])  #(wavelength in angstrom)


wavelength_micron_hfls3=wavelength_micron_hfls3/7.34  #redshift back to restframe (z=6.34)
wavelength_micron_hfls3=wavelength_micron_hfls3*1.0e6  #convert from m to micron
flux_mjy_hfls3=flux_mjy_hfls3*1000.0  #convert from Jy to mJy











template_file_eyelash = open("../cosmic_eyelash_sed_josh.txt", "r") 

lines_eyelash = template_file_eyelash.readlines()
template_file_eyelash.close()

response_curve_eyelash=[]

for i in range(0, len(lines_eyelash)):
    separated_lines_eyelash=lines_eyelash[i].split() 
    response_curve_eyelash.append(separated_lines_eyelash)


response_curve_eyelash = np.array(response_curve_eyelash)
wavelength_micron_eyelash=response_curve_eyelash[:,0]
flux_mjy_eyelash=response_curve_eyelash[:,1]


wavelength_micron_eyelash = np.array(wavelength_micron_eyelash)
flux_mjy_eyelash = np.array(flux_mjy_eyelash)
wavelength_micron_eyelash=np.array([float(i) for i in wavelength_micron_eyelash])  #(wavelength in angstrom)
flux_mjy_eyelash=np.array([float(i) for i in flux_mjy_eyelash])  #(wavelength in angstrom)



#-----------------------------------------------------------


L_distribution_all_aless=[]
L_distribution_all_uperr_aless=[]
L_distribution_all_lowerr_aless=[]

L_distribution_all_arp220=[]
L_distribution_all_uperr_arp220=[]
L_distribution_all_lowerr_arp220=[]

L_distribution_all_hfls3=[]
L_distribution_all_uperr_hfls3=[]
L_distribution_all_lowerr_hfls3=[]

L_distribution_all_eyelash=[]
L_distribution_all_uperr_eyelash=[]
L_distribution_all_lowerr_eyelash=[]




for l in range(0, len(cluster_list)):    #=================================Loop over clusters




    ##----------------------------Estimate the map size in deg^2 (using variance map!!)----------------------------------------
    #
    #if cluster_list[l]=='S2CLS' or cluster_list[l]=='S2CLS_random' or cluster_list[l]=='PCL1002':
    #    hdulist = pyfits.open('../../'+cluster_list[l]+'/'+cluster_list[l]+'_rms_170919.fits')      #==================================================================================================================HERE Orz
    #    #w = wcs.WCS(hdulist[0].header, hdulist)
    #    NAXIS1=hdulist[0].header['NAXIS1']
    #    NAXIS2=hdulist[0].header['NAXIS2']
    #    CDELT1=hdulist[0].header['CDELT1']
    #    CDELT2=hdulist[0].header['CDELT2']
    #    #CDELT1=-5.5555555555555E-4
    #    #CDELT2=0.00055555555555555
    #    scidata = hdulist[0].data  #================variance: unit of pW^2
    #    hdulist.close()
    #
    #else:
    #    hdulist = pyfits.open('../../'+cluster_list[l]+'/'+cluster_list[l]+'_flux_170919.fits')      #================================================================================================HERE
    #    #w = wcs.WCS(hdulist[1].header, hdulist)
    #    NAXIS1=hdulist[1].header['NAXIS1']
    #    NAXIS2=hdulist[1].header['NAXIS2']
    #    CDELT1=hdulist[1].header['CDELT1']
    #    CDELT2=hdulist[1].header['CDELT2']
    #    scidata = hdulist[1].data
    #    hdulist.close()
    #
    #
    #pix_value_array=[]
    #
    #for i in range(0, NAXIS2):
    #    for j in range(0, NAXIS1):
    #        if cluster_list[l]=='S2CLS' or cluster_list[l]=='S2CLS_random' or cluster_list[l]=='PCL1002':
    #            pix_value_array.append(scidata[i,j])   #Orz  ([i,j] for S2CLS; [0,i,j] for ordinary fields)
    #        else:
    #            pix_value_array.append(scidata[0,i,j])   #Orz  ([i,j] for S2CLS; [0,i,j] for ordinary fields)
    #
    #pix_value_array = np.array(pix_value_array)
    #pix_value_array=np.array([float(i) for i in pix_value_array])  #(wavelength in angstrom)
    #
    #pix_value_array = pix_value_array[~np.isnan(pix_value_array)]
    #
    ##print len(pix_value_array), 'pixels'
    #
    #area_map=len(pix_value_array)*abs(CDELT1)*abs(CDELT2)  #area in deg^2 of a cluster map
    #
    #map_size_all=map_size_all+area_map
    #



    #---------------Start magninalizing------------------------------------------




    print cluster_list[l]
    names_wht = pyfits.open('../../'+cluster_list[l]+'/'+cluster_list[l]+'_herschel_match_3p5sigma_cat_edgesourcedelete_allscuba2source_herschelflux.fits')    #===============================HERE
    names_wht_data = names_wht[1].data
    flux_350_observed_pre= names_wht_data.field('F350')
    source_name_pre= names_wht_data.field('name')
    SN_ratio_pre= names_wht_data.field('SN_ratio_1')
    SN_ratio_confusion_pre= names_wht_data.field('SN_ratio_confusion')

    flux_350_observed_pre = np.array(flux_350_observed_pre)
    flux_350_observed_pre=np.array([float(i) for i in flux_350_observed_pre])
    SN_ratio_pre = np.array(SN_ratio_pre)
    SN_ratio_pre=np.array([float(i) for i in SN_ratio_pre])
    SN_ratio_confusion_pre = np.array(SN_ratio_confusion_pre)
    SN_ratio_confusion_pre=np.array([float(i) for i in SN_ratio_confusion_pre])

    flux_350_observed=[]
    source_name=[]

    for i in range(0, len(flux_350_observed_pre)):

        if  math.isnan(flux_350_observed_pre[i]) == False and SN_ratio_confusion_pre[i]>=3.5:

            flux_350_observed.append(flux_350_observed_pre[i])
            source_name.append(source_name_pre[i])

    flux_350_observed = np.array(flux_350_observed)
    flux_350_observed=np.array([float(i) for i in flux_350_observed])

    if cluster_list[l]=='G12' or cluster_list[l]=='NGP1' or cluster_list[l]=='NGP2' or cluster_list[l]=='NGP3' or cluster_list[l]=='NGP4' or cluster_list[l]=='NGP5' or cluster_list[l]=='NGP6' or cluster_list[l]=='NGP7' or cluster_list[l]=='NGP8' or cluster_list[l]=='NGP9':
        flux_350_observed=flux_350_observed*1000.0

    #print len(source_name)
    #print '1111111111111111111111111111'



    for x in range(0, len(source_name)):    #===============================Loop over sources
    #for x in range(2, 4):    #===============================Loop over sources
       
        marginalized_array=[]
        print source_name[x]

        aless_array = np.load('./'+cluster_list[l]+'/'+str(source_name[x])+'_mcmc_values_aless.npy')
        arp220_array = np.load('./'+cluster_list[l]+'/'+str(source_name[x])+'_mcmc_values_arp220.npy')
        hfls3_array = np.load('./'+cluster_list[l]+'/'+str(source_name[x])+'_mcmc_values_hfls3.npy')
        eyelash_array = np.load('./'+cluster_list[l]+'/'+str(source_name[x])+'_mcmc_values_eyelash.npy')



        #------------Calculate luminosity, dust temperature, SFR for each source, before marginalize over templates!-------------------

        #execfile("./luminosity_etc_181114.py")
        #execfile("./luminosity_etc_190408.py")
        execfile("./luminosity_etc_190715.py")











        marginalized_array = np.concatenate((marginalized_array_random_aless, marginalized_array_random_arp220, marginalized_array_random_hfls3, marginalized_array_random_eyelash), axis=0)

        print np.shape(marginalized_array)
        #print marginalized_array
        print '#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'


        np.save('./marginalized_totalSFR/samples_values_all_'+str(cluster_list[l])+'_'+str(source_name[x])+'.npy', marginalized_array)
        #fig = corner.corner(marginalized_array, labels=["$z$", "$log(a)$", r"$log(L(L_{\odot}))$", r"$log(SFR(M_{\odot}/yr))$"],
        #                    quantiles=[0.16, 0.5, 0.84],
        #                    show_titles=True, title_kwargs={"fontsize": 12})
        
        #fig.savefig("./marginalized/triangle_marginalized_"+cluster_list[l]+"_"+str(source_name[x])+".png")  
        #plt.close(fig)


        lower_limit_z = corner.quantile(marginalized_array[:,0], q=[0.16], weights=None)
        fitted_median_z = corner.quantile(marginalized_array[:,0], q=[0.5], weights=None)
        higher_limit_z = corner.quantile(marginalized_array[:,0], q=[0.84], weights=None)
        z_distribution_all.append(fitted_median_z[0])
        z_distribution_all_uperr.append(higher_limit_z[0]-fitted_median_z[0])
        z_distribution_all_lowerr.append(fitted_median_z[0]-lower_limit_z[0])

        #print '2222222222222222222222222222'
        #print source_name[x]
        #print fitted_median_z
        #print lower_limit_z
        #print higher_limit_z


        lower_limit_a = corner.quantile(marginalized_array[:,1], q=[0.16], weights=None)
        fitted_median_a = corner.quantile(marginalized_array[:,1], q=[0.5], weights=None)
        higher_limit_a = corner.quantile(marginalized_array[:,1], q=[0.84], weights=None)
        a_distribution_all.append(fitted_median_a[0])
        a_distribution_all_uperr.append(higher_limit_a[0]-fitted_median_a[0])
        a_distribution_all_lowerr.append(fitted_median_a[0]-lower_limit_a[0])

        #print fitted_median_a
        #print lower_limit_a
        #print higher_limit_a


        #==================================II. Calculate Luminosity and SFR=================================================================
        #marginalized_array_random=[]

        #np.save('./marginalized/samples_values_all_'+str(cluster_list[l])+'_'+str(source_name[x])+'.npy', marginalized_array)
        #execfile("./luminosity_etc_181113.py")

        lower_limit_L = corner.quantile(marginalized_array[:,2], q=[0.16], weights=None)
        fitted_median_L = corner.quantile(marginalized_array[:,2], q=[0.5], weights=None)
        higher_limit_L = corner.quantile(marginalized_array[:,2], q=[0.84], weights=None)
        L_distribution_all.append(fitted_median_L[0])
        L_distribution_all_uperr.append(higher_limit_L[0]-fitted_median_L[0])
        L_distribution_all_lowerr.append(fitted_median_L[0]-lower_limit_L[0])
        #print fitted_median_L
        #print lower_limit_L
        #print higher_limit_L

        lower_limit_SFR = corner.quantile(marginalized_array[:,3], q=[0.16], weights=None)
        fitted_median_SFR = corner.quantile(marginalized_array[:,3], q=[0.5], weights=None)
        higher_limit_SFR = corner.quantile(marginalized_array[:,3], q=[0.84], weights=None)
        SFR_distribution_all.append(fitted_median_SFR[0])
        SFR_distribution_all_uperr.append(higher_limit_SFR[0]-fitted_median_SFR[0])
        SFR_distribution_all_lowerr.append(fitted_median_SFR[0]-lower_limit_SFR[0])
        #print fitted_median_SFR
        #print lower_limit_SFR
        #print higher_limit_SFR











mock_catalogue_aless=np.column_stack((L_distribution_all_aless, L_distribution_all_uperr_aless, L_distribution_all_lowerr_aless))
np.savetxt('./marginalized_totalSFR/luminosity_allsource_aless.cat', mock_catalogue_aless, delimiter=' ') #================================================================================HERE

mock_catalogue_arp220=np.column_stack((L_distribution_all_arp220, L_distribution_all_uperr_arp220, L_distribution_all_lowerr_arp220))
np.savetxt('./marginalized_totalSFR/luminosity_allsource_arp220.cat', mock_catalogue_arp220, delimiter=' ') #================================================================================HERE

mock_catalogue_hfls3=np.column_stack((L_distribution_all_hfls3, L_distribution_all_uperr_hfls3, L_distribution_all_lowerr_hfls3))
np.savetxt('./marginalized_totalSFR/luminosity_allsource_hfls3.cat', mock_catalogue_hfls3, delimiter=' ') #================================================================================HERE

mock_catalogue_eyelash=np.column_stack((L_distribution_all_eyelash, L_distribution_all_uperr_eyelash, L_distribution_all_lowerr_eyelash))
np.savetxt('./marginalized_totalSFR/luminosity_allsource_eyelash.cat', mock_catalogue_eyelash, delimiter=' ') #================================================================================HERE





#=====================================III. Plot all z, L, SFR distribution of all sources===================================================
#np.save('./z_distribution_all.npy', samples)
#np.save('./L_distribution_all.npy', samples)
#np.save('./SFR_distribution_all.npy', samples)


mock_catalogue_z=np.column_stack((z_distribution_all, z_distribution_all_uperr, z_distribution_all_lowerr))
np.savetxt('./marginalized_totalSFR/allsource_z.cat', mock_catalogue_z, delimiter=' ') #================================================================================HERE

mock_catalogue_L=np.column_stack((L_distribution_all, L_distribution_all_uperr, L_distribution_all_lowerr))
np.savetxt('./marginalized_totalSFR/allsource_L.cat', mock_catalogue_L, delimiter=' ') #================================================================================HERE

mock_catalogue_SFR=np.column_stack((SFR_distribution_all, SFR_distribution_all_uperr, SFR_distribution_all_lowerr))
np.savetxt('./marginalized_totalSFR/allsource_SFR.cat', mock_catalogue_SFR, delimiter=' ') #================================================================================HERE



print '========================Reaching this line is fine==============================='











