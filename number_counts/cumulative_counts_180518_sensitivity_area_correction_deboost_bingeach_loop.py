import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm
from astropy.io import fits
import astropy.wcs as wcs
from scipy.interpolate import interp1d


#cluster_list=['Bootes1', 'EGS', 'Lockman', 'G12', 'NGP1', 'NGP2', 'NGP3', 'NGP4', 'NGP5', 'NGP6', 'NGP7', 'NGP8', 'NGP9', 'S2CLS']
#cluster_list=['NGP9']
cluster_list=['Bootes1', 'EGS', 'Lockman', 'G12', 'NGP1', 'NGP2', 'NGP3', 'NGP4', 'NGP5', 'NGP6', 'NGP7', 'NGP8',  'NGP9']


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
    #execfile("./cumulative_counts_190329_sensitivity_area_correction_deboost_bingeach.py")
    execfile("./cumulative_counts_190710_sensitivity_area_correction_deboost_bingeach.py")




figa.delaxes(ax_all[2][3])
figa.delaxes(ax_all[2][4])

plt.show()
