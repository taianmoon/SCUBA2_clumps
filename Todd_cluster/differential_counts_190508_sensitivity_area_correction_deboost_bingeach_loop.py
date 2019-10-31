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
    execfile("./differential_counts_190508_sensitivity_area_correction_deboost_bingeach.py")
    #execfile("./differential_counts_181127_sensitivity_area_correction_deboost_bingeach.py")
    #execfile("./differential_counts_180904_sensitivity_area_correction_deboost_bingeach.py")
    #execfile("./differential_counts_180522_sensitivity_area_correction_deboost_bingeach.py")
    #execfile("./differential_counts_180518_sensitivity_area_correction_deboost_bingeach.py")


plt.show()


