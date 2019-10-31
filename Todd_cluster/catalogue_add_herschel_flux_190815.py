import pyfits
import matplotlib.pyplot as plt
import numpy as np
import math
import glob
import os.path
from astropy.io import fits
import astropy.wcs as wcs
from astropy.table import Table

cluster_list=['Planck18p194', 'Planck18p735', 'Planck24p194', 'PLCK_DU_G045.7-41.2', 'PLCK_DU_G059.1-67.1', 'PLCK_DU_G073.4-57.5', 'PLCK_G006.1+61.8', 'PLCK_G009.8+72.6', 'PLCK_G056.7+62.6', 'PLCK_G068.3+31.9', 'PLCK_G075.1+33.2', 'PLCK_G077.7+32.6', 'PLCK_G078.9+48.2', 'PLCK_G082.5+38.4', 'PLCK_G083.3+51.0', 'PLCK_G091.9+43.0', 'PLCK_G093.6+55.9', 'PLCK_G132.9-76.0', 'PLCK_G144.1+81.0', 'PLCK_G160.7+41.0', 'PLCK_G162.1-59.3', 'PLCK_G165.8+45.3', 'PLCK_G173.8+59.3', 'PLCK_G177.0+35.9', 'PLCK_G179.3+50.7', 'PLCK_G186.3-72.7', 'PLCK_G186.6+66.7', 'PLCK_G188.6-68.9', 'PLCK_G191.3+62.0', 'PLCK_G191.8-83.4', 'PLCK_G201.1+50.7', 'PLCK_G213.0+65.9', 'PLCK_G223.9+41.2', 'PLCK_G328.9+71.4', 'PLCK_G49.6-42.9', 'PLCK_G84.0-71.5', 'PLCK_HZ_G038.0-51.5', 'PLCK_HZ_G067.2-63.8', 'PLCK_HZ_G103.1-73.6', 'PLCK_HZ_G106.8-83.3', 'PLCK_HZ_G119.4-76.6', 'PLCK_HZ_G132.6-81.1', 'PLCK_HZ_G171.1-78.7', 'PLCK_HZ_G173.9+57.0', 'PLCK_HZ_G176.6+59.0', 'PLCK_HZ_G214.1+48.3']



for l in range(0, len(cluster_list)):  #-----------------loop over clusters
    
    print cluster_list[l]

    hdulist_250 = fits.open('../'+cluster_list[l]+'/250herschel.fits')      #========================================================================HERE
    w_flux_250 = wcs.WCS(hdulist_250[1].header, hdulist_250)
    scidata_flux_250 = hdulist_250[1].data
    w_noise_250 = wcs.WCS(hdulist_250[2].header, hdulist_250)
    scidata_noise_250 = hdulist_250[2].data
    hdulist_250.close()

    hdulist_350 = fits.open('../'+cluster_list[l]+'/350herschel.fits')      #========================================================================HERE
    w_flux_350 = wcs.WCS(hdulist_350[1].header, hdulist_350)
    scidata_flux_350 = hdulist_350[1].data
    w_noise_350 = wcs.WCS(hdulist_350[2].header, hdulist_350)
    scidata_noise_350 = hdulist_350[2].data
    hdulist_350.close()

    hdulist_500 = fits.open('../'+cluster_list[l]+'/500herschel.fits')      #========================================================================HERE
    w_flux_500 = wcs.WCS(hdulist_500[1].header, hdulist_500)
    scidata_flux_500 = hdulist_500[1].data
    w_noise_500 = wcs.WCS(hdulist_500[2].header, hdulist_500)
    scidata_noise_500 = hdulist_500[2].data
    hdulist_500.close()


    check_Herschel_table = Table.read('./edgesourcedelete_deboost_190812/'+cluster_list[l]+'_3p5sigma_deboost.cat', format='ascii')    


    flux_250_pre=[]
    error_250_pre=[]
    flux_350_pre=[]
    error_350_pre=[]
    flux_500_pre=[]
    error_500_pre=[]
    for ss in range(0, len(check_Herschel_table['RA'])):  #-------------------------loop over sources

        pix_x_250, pix_y_250 = w_flux_250.all_world2pix(check_Herschel_table['RA'][ss], check_Herschel_table['Dec'][ss], 1)
        pix_x_noise_250, pix_y_noise_250 = w_noise_250.all_world2pix(check_Herschel_table['RA'][ss], check_Herschel_table['Dec'][ss], 1)
        flux_250_pre.append(1000.0*scidata_flux_250[int(pix_y_250-1), int(pix_x_250-1)])  #note that x and y are reversed!!! x1000 to convert to mJy
        error_250_pre.append(1000.0*scidata_noise_250[int(pix_y_noise_250-1), int(pix_x_noise_250-1)])

        pix_x_350, pix_y_350 = w_flux_350.all_world2pix(check_Herschel_table['RA'][ss], check_Herschel_table['Dec'][ss], 1)
        pix_x_noise_350, pix_y_noise_350 = w_noise_350.all_world2pix(check_Herschel_table['RA'][ss], check_Herschel_table['Dec'][ss], 1)
        flux_350_pre.append(1000.0*scidata_flux_350[int(pix_y_350-1), int(pix_x_350-1)])  #note that x and y are reversed!!! x1000 to convert to mJy
        error_350_pre.append(1000.0*scidata_noise_350[int(pix_y_noise_350-1), int(pix_x_noise_350-1)])

        pix_x_500, pix_y_500 = w_flux_500.all_world2pix(check_Herschel_table['RA'][ss], check_Herschel_table['Dec'][ss], 1)
        pix_x_noise_500, pix_y_noise_500 = w_noise_500.all_world2pix(check_Herschel_table['RA'][ss], check_Herschel_table['Dec'][ss], 1)
        flux_500_pre.append(1000.0*scidata_flux_500[int(pix_y_500-1), int(pix_x_500-1)])  #note that x and y are reversed!!! x1000 to convert to mJy
        error_500_pre.append(1000.0*scidata_noise_500[int(pix_y_noise_500-1), int(pix_x_noise_500-1)])

    print flux_250_pre
    print error_250_pre
    print flux_350_pre
    print error_350_pre
    print flux_500_pre
    print error_500_pre
    

    #if math.isnan(flux_250_pre) == True:
    #    print 'yes!'


    #check_Herschel_table = Table.read('./edgesourcedelete_deboost_190812/'+cluster_list[l]+'_3p5sigma_deboost.cat', format='ascii')    


    check_Herschel_table['Herschel_250_mJy'] = flux_250_pre
    check_Herschel_table['Herschel_250err_mJy'] = error_250_pre
    check_Herschel_table['Herschel_350_mJy'] = flux_350_pre
    check_Herschel_table['Herschel_350err_mJy'] = error_350_pre
    check_Herschel_table['Herschel_500_mJy'] = flux_500_pre
    check_Herschel_table['Herschel_500err_mJy'] = error_500_pre


    check_Herschel_table.write('./edgesourcedelete_deboost_190812/'+cluster_list[l]+'_herschel_flux.fits', format='fits')   









