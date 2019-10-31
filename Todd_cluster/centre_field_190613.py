#import matplotlib.pyplot as plt
import numpy as np
import math
#import matplotlib.cm as cm
from astropy.io import fits
import astropy.wcs as wcs


#cut_radius = 350.0  #==============HERE: in arcsec



cluster_list=['Planck18p194', 'Planck18p735', 'Planck24p194', 'PLCK_DU_G045.7-41.2', 'PLCK_DU_G059.1-67.1', 'PLCK_DU_G073.4-57.5', 'PLCK_G006.1+61.8', 'PLCK_G009.8+72.6', 'PLCK_G056.7+62.6', 'PLCK_G068.3+31.9', 'PLCK_G075.1+33.2', 'PLCK_G077.7+32.6', 'PLCK_G078.9+48.2', 'PLCK_G082.5+38.4', 'PLCK_G083.3+51.0', 'PLCK_G091.9+43.0', 'PLCK_G093.6+55.9', 'PLCK_G132.9-76.0', 'PLCK_G144.1+81.0', 'PLCK_G160.7+41.0', 'PLCK_G162.1-59.3', 'PLCK_G165.8+45.3', 'PLCK_G173.8+59.3', 'PLCK_G177.0+35.9', 'PLCK_G179.3+50.7', 'PLCK_G186.3-72.7', 'PLCK_G186.6+66.7', 'PLCK_G188.6-68.9', 'PLCK_G191.3+62.0', 'PLCK_G191.8-83.4', 'PLCK_G201.1+50.7', 'PLCK_G213.0+65.9', 'PLCK_G223.9+41.2', 'PLCK_G328.9+71.4', 'PLCK_G49.6-42.9', 'PLCK_G84.0-71.5', 'PLCK_HZ_G038.0-51.5', 'PLCK_HZ_G067.2-63.8', 'PLCK_HZ_G103.1-73.6', 'PLCK_HZ_G106.8-83.3', 'PLCK_HZ_G119.4-76.6', 'PLCK_HZ_G132.6-81.1', 'PLCK_HZ_G171.1-78.7', 'PLCK_HZ_G173.9+57.0', 'PLCK_HZ_G176.6+59.0', 'PLCK_HZ_G214.1+48.3']




centre_array_RA=[]
centre_array_Dec=[]
cum_no=[]

for l in range (0, len(cluster_list)):  #=================loop over realization times
    
    print l
    cum_no.append(l)
    print cluster_list[l]
    hdulist = fits.open('./'+cluster_list[l]+'/mf_crop.fits')      #==========================================================================================================HERE
    w_flux = wcs.WCS(hdulist[0].header, hdulist)
    w_variance = wcs.WCS(hdulist[1].header, hdulist)

    CRVAL1_flux=hdulist[0].header['CRVAL1']
    CRVAL2_flux=hdulist[0].header['CRVAL2']

    centre_array_RA.append(CRVAL1_flux)
    centre_array_Dec.append(CRVAL2_flux)


print centre_array_RA
print centre_array_Dec



mock_catalogue=np.column_stack((cum_no, centre_array_RA, centre_array_Dec)) 
np.savetxt('./centre_field.cat', mock_catalogue, delimiter=' ', header='no RA Dec')
