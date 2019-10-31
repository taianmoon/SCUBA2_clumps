#import matplotlib.pyplot as plt
import numpy as np
#import math
#import matplotlib.cm as cm
from astropy.io import fits
import astropy.wcs as wcs
#from scipy.interpolate import interp1d


cluster_list=['Planck18p194', 'Planck18p735', 'Planck24p194', 'PLCK_DU_G045.7-41.2', 'PLCK_DU_G059.1-67.1', 'PLCK_DU_G073.4-57.5', 'PLCK_G006.1+61.8', 'PLCK_G009.8+72.6', 'PLCK_G056.7+62.6', 'PLCK_G068.3+31.9', 'PLCK_G075.1+33.2', 'PLCK_G077.7+32.6', 'PLCK_G078.9+48.2', 'PLCK_G082.5+38.4', 'PLCK_G083.3+51.0', 'PLCK_G091.9+43.0', 'PLCK_G093.6+55.9', 'PLCK_G132.9-76.0', 'PLCK_G144.1+81.0', 'PLCK_G160.7+41.0', 'PLCK_G162.1-59.3', 'PLCK_G165.8+45.3', 'PLCK_G173.8+59.3', 'PLCK_G177.0+35.9', 'PLCK_G179.3+50.7', 'PLCK_G186.3-72.7', 'PLCK_G186.6+66.7', 'PLCK_G188.6-68.9', 'PLCK_G191.3+62.0', 'PLCK_G191.8-83.4', 'PLCK_G201.1+50.7', 'PLCK_G213.0+65.9', 'PLCK_G223.9+41.2', 'PLCK_G328.9+71.4', 'PLCK_G49.6-42.9', 'PLCK_G84.0-71.5', 'PLCK_HZ_G038.0-51.5', 'PLCK_HZ_G067.2-63.8', 'PLCK_HZ_G103.1-73.6', 'PLCK_HZ_G106.8-83.3', 'PLCK_HZ_G119.4-76.6', 'PLCK_HZ_G132.6-81.1', 'PLCK_HZ_G171.1-78.7', 'PLCK_HZ_G173.9+57.0', 'PLCK_HZ_G176.6+59.0', 'PLCK_HZ_G214.1+48.3']



area_array=[]
area_array_cum=0.0

for l in range(0, len(cluster_list)):

    #print cluster_list[l]


    hdulist = fits.open('./'+cluster_list[l]+'/mf.fits')      #================================================================================HERE
    w = wcs.WCS(hdulist[1].header, hdulist)
    NAXIS1=hdulist[1].header['NAXIS1']
    NAXIS2=hdulist[1].header['NAXIS2']
    CDELT1=hdulist[1].header['CDELT1']
    CDELT2=hdulist[1].header['CDELT2']
    scidata = hdulist[1].data
    hdulist.close()


    pix_value_array=[]

    for i in range(0, NAXIS2):
        for j in range(0, NAXIS1):
            pix_value_array.append(scidata[0,i,j])

    pix_value_array = np.array(pix_value_array)
    pix_value_array=np.array([float(i) for i in pix_value_array])  #(wavelength in angstrom)

    pix_value_array = pix_value_array[~np.isnan(pix_value_array)]

    #print len(pix_value_array), 'pixels'

    area_map=len(pix_value_array)*abs(CDELT1)*abs(CDELT2)

    #print 'the total area from flux map (deg^2) is', area_map, 'degree^2'
    print area_map


    area_array.append(area_map)
    area_array_cum=area_array_cum+area_map

area_array = np.array(area_array)
area_array=np.array([float(i) for i in area_array])  #(wavelength in angstrom)

#print area_array
print 'the total area for all clusters is (deg^2):'
print area_array_cum

