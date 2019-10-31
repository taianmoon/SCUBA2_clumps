import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm
from astropy.io import fits
import astropy.wcs as wcs
from scipy.interpolate import interp1d


cluster_list=['Bootes1', 'EGS', 'G12', 'Lockman', 'NGP1', 'NGP2', 'NGP3', 'NGP4', 'NGP5', 'NGP6', 'NGP7', 'NGP8', 'NGP9', 'S2CLS']

#cluster_list=['Planck18p194']

binsize=2  #=================================================================================HERE: bin size 1 or 2?


if binsize==1:

    geach_flux_cum=[3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5]  #mJy
    geach_cum_number_cum=[1012.3, 508.0, 271.9, 151.8, 85.3, 47.1, 26.4, 14.5, 8.7, 5.5, 3.2, 2.4, 1.8]  #/deg^2
    geach_cum_number_err_up_cum=[19.6, 12.3, 8.5, 6.2, 4.7, 3.6, 2.8, 2.2, 1.8, 1.5, 1.2, 1.1, 1.0] #Poisson error
    geach_cum_number_err_low_cum=[19.2, 12.0, 8.2, 6.0, 4.4, 3.3, 2.5, 1.9, 1.5, 1.2, 0.9, 0.8, 0.7]


    geach_flux_diff=[3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5]  #mJy
    geach_diff_number_diff=[451.0, 204.4, 102.6, 56.1, 32.5, 18.0, 9.8, 5.8, 3.4, 2.1, 0.8, 0.5, 0.3]  #/deg^2 /mJy
    geach_diff_number_err_up_diff=[17.1, 9.3, 6.0, 4.3, 3.2, 2.5, 1.9, 1.5, 1.2, 1.1, 0.8, 0.7, 0.6] #Poisson error
    geach_diff_number_err_low_diff=[16.4, 8.9, 5.7, 4.0, 2.9, 2.2, 1.6, 1.2, 0.9, 0.7, 0.4, 0.3, 0.2]


if binsize==2:

    geach_flux_cum=[4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0]  #mJy
    geach_cum_number_cum=[760.15, 211.85, 66.2, 20.45, 7.1, 2.8, 1.8]  #/deg^2
    geach_cum_number_err_up_cum=[11.57, 5.26, 2.96, 1.78, 1.17, 0.81, 1.0] #Poisson error
    geach_cum_number_err_low_cum=[11.32, 5.08, 2.75, 1.57, 0.96, 0.60, 0.7]


    geach_flux_diff=[4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0]  #mJy
    geach_diff_number_diff=[327.7, 79.35, 25.25, 7.8, 2.75, 0.65, 0.3]  #/deg^2 /mJy
    geach_diff_number_err_up_diff=[9.73, 3.69, 2.03, 1.21, 0.81, 0.53, 0.6] #Poisson error
    geach_diff_number_err_low_diff=[9.33, 3.48, 1.82, 1.0, 0.57, 0.25, 0.2]






plt.close('all')
plt.figure()

for l in range(0, len(cluster_list)):

    print cluster_list[l]
    execfile("./overdensity_calculation_180509.py")



#plt.tick_params(width=2, length=16, which='major')
#plt.tick_params(width=2, length=5, which='minor')

#plt.xscale('log',nonposy='clip')
#plt.yscale('log',nonposy='clip')

#plt.legend(loc=4)
#plt.grid()
#plt.xlabel('Flux (mJy)')  
#plt.set_title(str(source_name[i]))
#plt.ylabel('overdensity (sigma)')  
#plt.rc('font', size=15)
#plt.show()






