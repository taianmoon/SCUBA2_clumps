##import matplotlib.pyplot as plt
#import numpy as np
#import math
##import matplotlib.cm as cm
#from astropy.io import fits
#import astropy.wcs as wcs



#cut_radius = 350.0  #==============HERE: in arcsec



hdulist = fits.open('./'+cluster_list[l]+'/mf.fits')      #==========================================================================================================HERE
w_flux = wcs.WCS(hdulist[0].header, hdulist)
w_variance = wcs.WCS(hdulist[1].header, hdulist)
NAXIS1_flux=hdulist[0].header['NAXIS1']
NAXIS2_flux=hdulist[0].header['NAXIS2']
CDELT1_flux=hdulist[0].header['CDELT1']
CDELT2_flux=hdulist[0].header['CDELT2']
CRPIX1_flux=hdulist[0].header['CRPIX1']
CRPIX2_flux=hdulist[0].header['CRPIX2']
CRVAL1_flux=hdulist[0].header['CRVAL1']
CRVAL2_flux=hdulist[0].header['CRVAL2']
NAXIS1_variance=hdulist[1].header['NAXIS1']
NAXIS2_variance=hdulist[1].header['NAXIS2']
CDELT1_variance=hdulist[1].header['CDELT1']
CDELT2_variance=hdulist[1].header['CDELT2']
CRPIX1_variance=hdulist[1].header['CRPIX1']
CRPIX2_variance=hdulist[1].header['CRPIX2']
CRVAL1_variance=hdulist[1].header['CRVAL1']
CRVAL2_variance=hdulist[1].header['CRVAL2']
scidata_flux = hdulist[0].data
header_flux = hdulist[0].header
scidata_variance = hdulist[1].data
header_variance = hdulist[1].header





for i in range(0, NAXIS2_flux):
    for j in range(0, NAXIS1_flux):
        if np.sqrt(pow(j-(CRPIX1_flux-1.0), 2)+pow(i-(CRPIX2_flux-1.0), 2)) >= (cut_radius/3600.0)/abs(CDELT1_flux):
            scidata_flux[0, i, j] = 'nan'

#print scidata_flux[0, 100, :]

for i in range(0, NAXIS2_variance):
    for j in range(0, NAXIS1_variance):
        if np.sqrt(pow(j-(CRPIX1_variance-1.0), 2)+pow(i-(CRPIX2_variance-1.0), 2)) >= (cut_radius/3600.0)/abs(CDELT1_variance):
            scidata_variance[0, i, j] = 'nan'

#print scidata_variance[0, 100, :]




hdulist.writeto('./'+cluster_list[l]+'/mf_crop.fits')




hdulist.close()













hdulist_sn = fits.open('./'+cluster_list[l]+'/mfsnr.fits')      #==========================================================================================================HERE
w_flux_sn = wcs.WCS(hdulist_sn[0].header, hdulist_sn)
w_variance_sn = wcs.WCS(hdulist_sn[1].header, hdulist_sn)
NAXIS1_flux_sn=hdulist_sn[0].header['NAXIS1']
NAXIS2_flux_sn=hdulist_sn[0].header['NAXIS2']
CDELT1_flux_sn=hdulist_sn[0].header['CDELT1']
CDELT2_flux_sn=hdulist_sn[0].header['CDELT2']
CRPIX1_flux_sn=hdulist_sn[0].header['CRPIX1']
CRPIX2_flux_sn=hdulist_sn[0].header['CRPIX2']
CRVAL1_flux_sn=hdulist_sn[0].header['CRVAL1']
CRVAL2_flux_sn=hdulist_sn[0].header['CRVAL2']
NAXIS1_variance_sn=hdulist_sn[1].header['NAXIS1']
NAXIS2_variance_sn=hdulist_sn[1].header['NAXIS2']
CDELT1_variance_sn=hdulist_sn[1].header['CDELT1']
CDELT2_variance_sn=hdulist_sn[1].header['CDELT2']
CRPIX1_variance_sn=hdulist_sn[1].header['CRPIX1']
CRPIX2_variance_sn=hdulist_sn[1].header['CRPIX2']
CRVAL1_variance_sn=hdulist_sn[1].header['CRVAL1']
CRVAL2_variance_sn=hdulist_sn[1].header['CRVAL2']
scidata_flux_sn = hdulist_sn[0].data
header_flux_sn = hdulist_sn[0].header
scidata_variance_sn = hdulist_sn[1].data
header_variance_sn = hdulist_sn[1].header



#print scidata_flux[0,37,142]
#print scidata_variance[0,37,142]





for i in range(0, NAXIS2_flux_sn):
    for j in range(0, NAXIS1_flux_sn):
        if np.sqrt(pow(j-(CRPIX1_flux_sn-1.0), 2)+pow(i-(CRPIX2_flux_sn-1.0), 2)) >= (cut_radius/3600.0)/abs(CDELT1_flux_sn):
            scidata_flux_sn[0, i, j] = 'nan'

#print scidata_flux_sn[0, 100, :]

for i in range(0, NAXIS2_variance_sn):
    for j in range(0, NAXIS1_variance_sn):
        if np.sqrt(pow(j-(CRPIX1_variance_sn-1.0), 2)+pow(i-(CRPIX2_variance_sn-1.0), 2)) >= (cut_radius/3600.0)/abs(CDELT1_variance_sn):
            scidata_variance_sn[0, i, j] = 'nan'

#print scidata_variance_sn[0, 100, :]


hdulist_sn.writeto('./'+cluster_list[l]+'/mfsnr_crop.fits')


hdulist_sn.close()














