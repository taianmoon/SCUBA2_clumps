from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm
from astropy.io import fits
import astropy.wcs as wcs
from scipy.interpolate import interp1d
import random
import os
import statistics



cluster_list=['PLCK_HZ_G214.1+48.3']  #can't do more clusters because of code constraints #=========================================================================================HERE
sensitivity_switch=1 #===========================HERE: 1--confine to central region, 2--divide the sensitivity area curve
realization_times=1500  #=========================================================================================HERE
radius_sources=(350.0-30.0)/60.0  #=========================================================================================HERE:arcmin, radius
number_of_sources=10  #=========================================================================================HERE: number of mock sources in the map




random_flux_mJy=[]#mJy
random_flux_mJy_notboost=[2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0]#mJy
random_flux_mJy_int=[2, 4, 6, 8, 10, 12, 14, 16, 18, 20]#mJy



random_flux=np.array([2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0]) #mJy for Todd's clusters
#print random_flux


#sn_ratio_average=[3.85, 3.63, 5.01, 3.99, 4.32, 4.57, 5.36, 4.19, 4.18, 4.16] #from 2.0 to 20.0 mJy
sn_ratio_average=[3.85, 3.22, 3.61, 3.92, 4.33, 4.41, 5.17, 5.07, 8.17, 5.12] #from 2.0 to 20.0 mJy #values for Todd's clusters

for sa in range (0, len(sn_ratio_average)):
    random_flux_mJy.append(random_flux_mJy_notboost[sa]/(1.0+0.2*pow(sn_ratio_average[sa]/5.0,-2.3)))

print random_flux_mJy
print ';;;;;;;;;;;;;;;;;;;;;'

completeness_level_2=[]
completeness_level_4=[]
completeness_level_6=[]
completeness_level_8=[]
completeness_level_10=[]
completeness_level_12=[]
completeness_level_14=[]
completeness_level_16=[]
completeness_level_18=[]
completeness_level_20=[]

completeness_level_median=[]  #from 2 to 20 mJy
completeness_level_mean=[]  #from 2 to 20 mJy
completeness_level_std=[]  #from 2 to 20 mJy
completeness_level_cumulative=[]


completeness_level_cumulative_loop_2=0
completeness_level_cumulative_loop_4=0
completeness_level_cumulative_loop_6=0
completeness_level_cumulative_loop_8=0
completeness_level_cumulative_loop_10=0
completeness_level_cumulative_loop_12=0
completeness_level_cumulative_loop_14=0
completeness_level_cumulative_loop_16=0
completeness_level_cumulative_loop_18=0
completeness_level_cumulative_loop_20=0

for rt in range (0, realization_times):  #=================loop over realization times
    
    print 'current realization: '+str(rt)+' out of '+str(realization_times-1)
    #execfile("./completeness_190226.py")
    execfile("./completeness_190326.py")






completeness_level_median.append(statistics.median(completeness_level_2))
completeness_level_median.append(statistics.median(completeness_level_4))
completeness_level_median.append(statistics.median(completeness_level_6))
completeness_level_median.append(statistics.median(completeness_level_8))
completeness_level_median.append(statistics.median(completeness_level_10))
completeness_level_median.append(statistics.median(completeness_level_12))
completeness_level_median.append(statistics.median(completeness_level_14))
completeness_level_median.append(statistics.median(completeness_level_16))
completeness_level_median.append(statistics.median(completeness_level_18))
completeness_level_median.append(statistics.median(completeness_level_20))

completeness_level_mean.append(statistics.mean(completeness_level_2))
completeness_level_mean.append(statistics.mean(completeness_level_4))
completeness_level_mean.append(statistics.mean(completeness_level_6))
completeness_level_mean.append(statistics.mean(completeness_level_8))
completeness_level_mean.append(statistics.mean(completeness_level_10))
completeness_level_mean.append(statistics.mean(completeness_level_12))
completeness_level_mean.append(statistics.mean(completeness_level_14))
completeness_level_mean.append(statistics.mean(completeness_level_16))
completeness_level_mean.append(statistics.mean(completeness_level_18))
completeness_level_mean.append(statistics.mean(completeness_level_20))

#completeness_level_std.append(statistics.stdev(completeness_level_2))
#completeness_level_std.append(statistics.stdev(completeness_level_4))
#completeness_level_std.append(statistics.stdev(completeness_level_6))
#completeness_level_std.append(statistics.stdev(completeness_level_8))
#completeness_level_std.append(statistics.stdev(completeness_level_10))
#completeness_level_std.append(statistics.stdev(completeness_level_12))
#completeness_level_std.append(statistics.stdev(completeness_level_14))
#completeness_level_std.append(statistics.stdev(completeness_level_16))
#completeness_level_std.append(statistics.stdev(completeness_level_18))
#completeness_level_std.append(statistics.stdev(completeness_level_20))

if abs(completeness_level_cumulative_loop_2) != 0:
    completeness_level_cumulative.append(math.sqrt(abs(completeness_level_cumulative_loop_2))/abs(completeness_level_cumulative_loop_2))
else:
    completeness_level_cumulative.append(1.0)
if abs(completeness_level_cumulative_loop_4) != 0:
    completeness_level_cumulative.append(math.sqrt(abs(completeness_level_cumulative_loop_4))/abs(completeness_level_cumulative_loop_4))
else:
    completeness_level_cumulative.append(1.0)
if abs(completeness_level_cumulative_loop_6) != 0:
    completeness_level_cumulative.append(math.sqrt(abs(completeness_level_cumulative_loop_6))/abs(completeness_level_cumulative_loop_6))
else:
    completeness_level_cumulative.append(1.0)
if abs(completeness_level_cumulative_loop_8) != 0:
    completeness_level_cumulative.append(math.sqrt(abs(completeness_level_cumulative_loop_8))/abs(completeness_level_cumulative_loop_8))
else:
    completeness_level_cumulative.append(1.0)
if abs(completeness_level_cumulative_loop_10) != 0:
    completeness_level_cumulative.append(math.sqrt(abs(completeness_level_cumulative_loop_10))/abs(completeness_level_cumulative_loop_10))
else:
    completeness_level_cumulative.append(1.0)
if abs(completeness_level_cumulative_loop_12) != 0:
    completeness_level_cumulative.append(math.sqrt(abs(completeness_level_cumulative_loop_12))/abs(completeness_level_cumulative_loop_12))
else:
    completeness_level_cumulative.append(1.0)
if abs(completeness_level_cumulative_loop_14) != 0:
    completeness_level_cumulative.append(math.sqrt(abs(completeness_level_cumulative_loop_14))/abs(completeness_level_cumulative_loop_14))
else:
    completeness_level_cumulative.append(1.0)
if abs(completeness_level_cumulative_loop_16) != 0:
    completeness_level_cumulative.append(math.sqrt(abs(completeness_level_cumulative_loop_16))/abs(completeness_level_cumulative_loop_16))
else:
    completeness_level_cumulative.append(1.0)
if abs(completeness_level_cumulative_loop_18) != 0:
    completeness_level_cumulative.append(math.sqrt(abs(completeness_level_cumulative_loop_18))/abs(completeness_level_cumulative_loop_18))
else:
    completeness_level_cumulative.append(1.0)
if abs(completeness_level_cumulative_loop_20) != 0:
    completeness_level_cumulative.append(math.sqrt(abs(completeness_level_cumulative_loop_20))/abs(completeness_level_cumulative_loop_20))
else:
    completeness_level_cumulative.append(0.0)



print completeness_level_median
print completeness_level_mean
#print completeness_level_std
print completeness_level_cumulative



#completeness_catalogue=np.column_stack((random_flux_mJy, completeness_level_median, completeness_level_mean, completeness_level_std))
completeness_catalogue=np.column_stack((random_flux_mJy, completeness_level_median, completeness_level_mean, completeness_level_cumulative))

#np.savetxt('../'+cluster_list[0]+'/'+cluster_list[0]+'_completeness_level.cat', completeness_catalogue, delimiter=' ', header="flux_mJy completeness_level_median completeness_level_mean completeness_level_cumulative") 
np.savetxt('./190911_wholemap/'+cluster_list[0]+'_completeness_level.cat', completeness_catalogue, delimiter=' ', header="flux_mJy completeness_level_median completeness_level_mean completeness_level_cumulative")



#execfile("./interpolate.py")



#=============================Plotting=========================================================

plt.plot(random_flux_mJy, completeness_level_median, '-o', color='blue', linewidth=2.0, markersize=8.0, label='median')
plt.plot(random_flux_mJy, completeness_level_mean, '-o', color='red', linewidth=2.0, markersize=8.0, label='mean')
#plt.errorbar(random_flux_mJy, completeness_level_median, yerr=completeness_level_std, color='b',fmt='o', capsize=5, elinewidth=2,markeredgewidth=2)
#plt.errorbar(random_flux_mJy, completeness_level_mean, yerr=completeness_level_std, color='r',fmt='o', capsize=5, elinewidth=2,markeredgewidth=2)
plt.errorbar(random_flux_mJy, completeness_level_median, yerr=completeness_level_cumulative, color='b',fmt='o', capsize=5, elinewidth=2,markeredgewidth=2)
plt.errorbar(random_flux_mJy, completeness_level_mean, yerr=completeness_level_cumulative, color='r',fmt='o', capsize=5, elinewidth=2,markeredgewidth=2)


plt.grid()

plt.legend(loc=4)
plt.xlabel('Flux density (mJy)')  #==========================================================HERE
plt.title(cluster_list[0])  #==========================================================HERE
plt.ylabel('Completeness')  #==========================================================HERE
plt.rc('font', size=30)


os.system('spd-say "your program has finished"')
plt.show()







