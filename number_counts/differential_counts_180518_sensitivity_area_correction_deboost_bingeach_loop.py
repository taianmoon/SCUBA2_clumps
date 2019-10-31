import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm
from astropy.io import fits
import astropy.wcs as wcs
from scipy.interpolate import interp1d



#cluster_list=['Bootes1', 'EGS', 'Lockman', 'G12', 'NGP1', 'NGP2', 'NGP3', 'NGP4', 'NGP5', 'NGP6', 'NGP7', 'NGP8', 'NGP9', 'S2CLS']
cluster_list=['Bootes1', 'EGS', 'Lockman', 'G12', 'NGP1', 'NGP2', 'NGP3', 'NGP4', 'NGP5', 'NGP6', 'NGP7', 'NGP8', 'NGP9']
#cluster_list=['Bootes1']


binsize=2  #=================================================================================HERE: bin size 1 or 2?
sn_want=3.5  #========================================================================HERE: change threshold for sn ratio







#====================================================Adding points from literature============================================

#geach_flux=[3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5]  #mJy

if binsize==1:
    geach_diff_number=[451.0, 204.4, 102.6, 56.1, 32.5, 18.0, 9.8, 5.8, 3.4, 2.1, 0.8, 0.5, 0.3]  #/deg^2 /mJy
    geach_diff_number_err_up=[17.1, 9.3, 6.0, 4.3, 3.2, 2.5, 1.9, 1.5, 1.2, 1.1, 0.8, 0.7, 0.6] #Poisson error
    geach_diff_number_err_low=[16.4, 8.9, 5.7, 4.0, 2.9, 2.2, 1.6, 1.2, 0.9, 0.7, 0.4, 0.3, 0.2]

if binsize==2:
    geach_diff_number=[327.7, 79.35, 25.25, 7.8, 2.75, 0.65, 0.3]  #/deg^2 /mJy
    geach_diff_number_err_up=[9.73, 3.69, 2.03, 1.21, 0.81, 0.53, 0.6] #Poisson error
    geach_diff_number_err_low=[9.33, 3.48, 1.82, 1.0, 0.57, 0.25, 0.2]


asymmetric_error = [geach_diff_number_err_low, geach_diff_number_err_up]



figz = plt.figure(figsize=(16.0, 10.0))
#figa = plt.figure(figsize=(16.0, 10.0))
#figz, zx_all =plt.subplots(3,5, figsize=(16.0, 10.0), sharex=True, sharey=True)
figa, ax_all =plt.subplots(3,5, figsize=(16.0, 10.0), sharex=True, sharey=True)

for l in range(0, len(cluster_list)):

    print cluster_list[l]
    execfile("./differential_counts_190710_sensitivity_area_correction_deboost_bingeach.py")
    #execfile("./differential_counts_190226_sensitivity_area_correction_deboost_bingeach.py")
    #execfile("./differential_counts_181127_sensitivity_area_correction_deboost_bingeach.py")
    #execfile("./differential_counts_180904_sensitivity_area_correction_deboost_bingeach.py")
    #execfile("./differential_counts_180522_sensitivity_area_correction_deboost_bingeach.py")
    #execfile("./differential_counts_180518_sensitivity_area_correction_deboost_bingeach.py")


plt.show()


