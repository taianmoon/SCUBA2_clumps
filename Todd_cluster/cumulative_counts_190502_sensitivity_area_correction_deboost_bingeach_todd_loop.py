import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm
from astropy.io import fits
import astropy.wcs as wcs
from scipy.interpolate import interp1d



#cluster_list=['Planck18p194', 'Planck18p735', 'Planck24p194', 'PLCK_DU_G045.7-41.2', 'PLCK_DU_G059.1-67.1', 'PLCK_DU_G073.4-57.5', 'PLCK_G006.1+61.8', 'PLCK_G009.8+72.6', 'PLCK_G056.7+62.6', 'PLCK_G068.3+31.9', 'PLCK_G075.1+33.2', 'PLCK_G077.7+32.6', 'PLCK_G078.9+48.2', 'PLCK_G082.5+38.4', 'PLCK_G083.3+51.0']

#cluster_list=['PLCK_G091.9+43.0', 'PLCK_G093.6+55.9', 'PLCK_G132.9-76.0', 'PLCK_G144.1+81.0', 'PLCK_G160.7+41.0', 'PLCK_G162.1-59.3', 'PLCK_G165.8+45.3', 'PLCK_G173.8+59.3', 'PLCK_G177.0+35.9', 'PLCK_G179.3+50.7', 'PLCK_G186.3-72.7', 'PLCK_G186.6+66.7', 'PLCK_G188.6-68.9', 'PLCK_G191.3+62.0', 'PLCK_G191.8-83.4']

#cluster_list=['PLCK_G201.1+50.7', 'PLCK_G213.0+65.9', 'PLCK_G223.9+41.2', 'PLCK_G328.9+71.4', 'PLCK_G49.6-42.9', 'PLCK_G84.0-71.5', 'PLCK_HZ_G038.0-51.5', 'PLCK_HZ_G067.2-63.8', 'PLCK_HZ_G103.1-73.6', 'PLCK_HZ_G106.8-83.3', 'PLCK_HZ_G119.4-76.6', 'PLCK_HZ_G132.6-81.1', 'PLCK_HZ_G171.1-78.7', 'PLCK_HZ_G173.9+57.0', 'PLCK_HZ_G176.6+59.0']

cluster_list=['PLCK_HZ_G214.1+48.3']




#cluster_list=['Bootes1', 'EGS', 'Lockman', 'G12', 'NGP1', 'NGP2', 'NGP3', 'NGP4', 'NGP5', 'NGP6', 'NGP7', 'NGP8', 'NGP9', 'S2CLS']
#cluster_list=['Bootes1', 'EGS', 'G12']
#cluster_list=['Bootes1', 'EGS', 'Lockman', 'G12', 'NGP1', 'NGP2', 'NGP3', 'NGP4', 'NGP5', 'NGP6', 'NGP7', 'NGP8',  'NGP9']


binsize=2  #=================================================================================HERE: bin size 1 or 2?
sn_want=3.5  #========================================================================HERE: change threshold for sn ratio




#====================================================Adding points from literature============================================

if binsize==1:
    geach_flux=[3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5]  #mJy
    geach_cum_number=[1012.3, 508.0, 271.9, 151.8, 85.3, 47.1, 26.4, 14.5, 8.7, 5.5, 3.2, 2.4, 1.8]  #/deg^2
    geach_cum_number_err_up=[19.6, 12.3, 8.5, 6.2, 4.7, 3.6, 2.8, 2.2, 1.8, 1.5, 1.2, 1.1, 1.0] #Poisson error
    geach_cum_number_err_low=[19.2, 12.0, 8.2, 6.0, 4.4, 3.3, 2.5, 1.9, 1.5, 1.2, 0.9, 0.8, 0.7]


if binsize==2:
    geach_flux=[4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0]  #mJy
    geach_cum_number=[760.15, 211.85, 66.2, 20.45, 7.1, 2.8, 1.8]  #/deg^2
    geach_cum_number_err_up=[11.57, 5.26, 2.96, 1.78, 1.17, 0.81, 1.0] #Poisson error
    geach_cum_number_err_low=[11.32, 5.08, 2.75, 1.57, 0.96, 0.60, 0.7]

asymmetric_error = [geach_cum_number_err_low, geach_cum_number_err_up]



#figa = plt.figure()

#ax =  [[ax1, ax2, ax3], [ax4, ax5, ax6], [ax7, ax8, ax9], [ax10, ax11, ax12], [ax13, ax14, ax15]]
figa, ax_all =plt.subplots(3,5, figsize=(16.0, 10.0), sharex=True, sharey=True)
#figa, ax = plt.subplots(3,5, sharex=True, sharey=True)

for l in range(0, len(cluster_list)):


    print cluster_list[l]
    #execfile("./cumulative_counts_180524_sensitivity_area_correction_deboost_bingeach.py")
    #execfile("./cumulative_counts_180831_sensitivity_area_correction_deboost_bingeach.py")
    #execfile("./cumulative_counts_181127_sensitivity_area_correction_deboost_bingeach.py")
    #execfile("./cumulative_counts_190226_sensitivity_area_correction_deboost_bingeach.py")
    execfile("./cumulative_counts_190329_sensitivity_area_correction_deboost_bingeach.py")


plt.show()
