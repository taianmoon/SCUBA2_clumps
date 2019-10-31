#import matplotlib.pyplot as plt
#import numpy as np
##import math
##import matplotlib.cm as cm
##from astropy.io import fits
##import astropy.wcs as wcs
#from scipy.interpolate import interp1d
##import pyfits
##import os
##import subprocess
#import scipy


#===========================z======================================


z_array = np.load('./'+cluster_list[l]+'/'+cluster_list[l]+'_marginalize_all_z.npy')
z_grid=z_array[:,0]
z_prob=z_array[:,1]


print z_grid
print z_prob



z_area = np.trapz(z_prob, z_grid)

print 0.5*z_area
print 0.16*z_area
print 0.84*z_area
print '#####################################################'




cum_area_z=[]
z_grid_cum=[]
z_prob_cum=[]
for i in range(0, len(z_grid)):

    z_grid_cum.append(z_grid[i])
    z_prob_cum.append(z_prob[i])
    cum_area_z.append(np.trapz(z_prob_cum, z_grid_cum))

#print z_grid_cum
#print z_prob_cum
#print cum_area_z


z_median = np.interp(0.5*z_area, cum_area_z, z_grid_cum)
z_16perc = np.interp(0.16*z_area, cum_area_z, z_grid_cum)
z_84perc = np.interp(0.84*z_area, cum_area_z, z_grid_cum)

print z_median
print z_16perc
print z_84perc


#plt.close('all')
#plt.figure()

#plt.plot(z_grid, cum_area_z,  color='blue', linewidth=1.5)
##plt.xscale('log',nonposy='clip')
#plt.yscale('log',nonposy='clip')
#plt.show()
    





#plt.close('all')
#plt.figure()

zx = figz.add_subplot(3,5,l+1)

zx.plot(z_grid, z_prob,  color='blue', linewidth=1.5, label=r'$z = {:.2f}'.format(z_median)+'^{+'+'{:.2f}'.format(z_84perc-z_median)+'}_{-'+'{:.2f}'.format(z_median-z_16perc)+'}$')
zx.axvline(x=z_median, color='k', linestyle='dashed', linewidth=1.0)
zx.axvline(x=z_16perc, color='k', linestyle='dashed', linewidth=1.0)
zx.axvline(x=z_84perc, color='k', linestyle='dashed', linewidth=1.0)
#zx.set_xlim(z_median-1.0, z_median+1.0)
zx.set_xlim(right=6.0)
#zx.set_ylim((1e-15, max(z_prob)*10.0))
#zx.text(z_median+0.9, max(z_prob)*1.0, r'$z = '+str(round(z_median,2))+'^{+'+str(round(z_84perc-z_median,2))+'}_{-'+str(round(z_median-z_16perc,2))+'}$', ha='right', va='top', size=12.0)
#zx.text(z_median-0.9, max(z_prob)*1.0, cluster_list[l], ha='left', va='top', size=12.0)


zx.tick_params(width=2, length=16, which='major')
#zx.tick_params(width=2, length=5, which='minor')




#zx.set_xscale('log',nonposy='clip')
zx.set_yscale('log',nonposy='clip')

#zx.legend(handlelength=0, handletextpad=0, fancybox=True, loc='lower left')
leg2 = zx.legend(handlelength=0, handletextpad=0, fancybox=True, loc='lower left')
print leg2
for item in leg2.legendHandles:
    item.set_visible(False)
leg2.get_frame().set_linewidth(0.0)
leg2.get_frame().set_alpha(0.0)


#zx.grid()
zx.set_xlabel('z')  
zx.set_title(cluster_list[l])
#zx.ylabel('P(z|F), all sources')  
plt.rc('font', size=18)

zx.tick_params(axis='y', which='both', bottom=False, top=False)
zx.set_yticks([])


#figz.savefig('./photoz_z_all.eps')
#figz.show()





#===========================a======================================


a_array = np.load('./'+cluster_list[l]+'/'+cluster_list[l]+'_marginalize_all_a.npy')
a_grid=a_array[:,0]
a_prob=a_array[:,1]


print a_grid
print a_prob



a_area = np.trapz(a_prob, a_grid)

print 0.5*a_area
print 0.16*a_area
print 0.84*a_area
print '#####################################################'




cum_area_a=[]
a_grid_cum=[]
a_prob_cum=[]
for i in range(0, len(a_grid)):

    a_grid_cum.append(a_grid[i])
    a_prob_cum.append(a_prob[i])
    cum_area_a.append(np.trapz(a_prob_cum, a_grid_cum))

#print a_grid_cum
#print a_prob_cum
#print cum_area_a


a_median = np.interp(0.5*a_area, cum_area_a, a_grid_cum)
a_16perc = np.interp(0.16*a_area, cum_area_a, a_grid_cum)
a_84perc = np.interp(0.84*a_area, cum_area_a, a_grid_cum)

print a_median
print a_16perc
print a_84perc


#plt.close('all')
#plt.figure()

#plt.plot(a_grid, cum_area_a,  color='blue', linewidth=1.5)
#plt.xscale('log',nonposy='clip')
#plt.yscale('log',nonposy='clip')
#plt.show()
    





#plt.close('all')
#plt.figure()

ax = figa.add_subplot(3,5,l+1)

ax.plot(a_grid, a_prob,  color='blue', linewidth=1.5, label=r'$a = {:.2f}'.format(a_median)+'^{+'+'{:.2f}'.format(a_84perc-a_median)+'}_{-'+'{:.2f}'.format(a_median-a_16perc)+'}$')
ax.axvline(x=a_median, color='k', linestyle='dashed', linewidth=1.0)
ax.axvline(x=a_16perc, color='k', linestyle='dashed', linewidth=1.0)
ax.axvline(x=a_84perc, color='k', linestyle='dashed', linewidth=1.0)
#ax.set_xlim(a_median-1.0, a_median+1.0)
ax.set_xlim(left=1e-1)
#ax.set_ylim((1e-15, max(a_prob)*10.0))
#ax.text(a_median+0.9, max(a_prob)*1.0, r'$a = '+str(round(a_median,2))+'^{+'+str(round(a_84perc-a_median,2))+'}_{-'+str(round(a_median-a_16perc,2))+'}$', ha='right', va='top', size=12.0)
#ax.text(a_median+0.9, max(a_prob)*1.0, r'$a = {:.2f}'.format(a_median)+'^{'+'{:.2f}'.format(a_84perc-a_median)+'}_{'+'{:.2f}'.format(a_median-a_16perc)+'}$', ha='right', va='top', size=12.0)
#ax.text(a_median-0.9, max(a_prob)*1.0, cluster_list[l], ha='left', va='top', size=12.0)

ax.tick_params(width=2, length=16, which='major')
#ax.tick_params(width=2, length=5, which='minor')

#r'$a = {:.3f}'.format(a_median)


ax.set_xscale('log',nonposy='clip')
ax.set_yscale('log',nonposy='clip')

#ax.legend(loc='lower left')
leg = ax.legend(handlelength=0, handletextpad=0, fancybox=True, loc='lower left')
print leg
print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
for item in leg.legendHandles:
    item.set_visible(False)
leg.get_frame().set_linewidth(0.0)
leg.get_frame().set_alpha(0.0)



#ax.grid()
ax.set_xlabel('a')  
ax.set_title(cluster_list[l])
#ax.set_ylabel('P(z|F), all sources')  
plt.rc('font', size=18)

ax.tick_params(axis='y', which='both', bottom=False, top=False)
ax.set_yticks([])


#figa.savefig('./photoz_a_all.eps')
#figa.show()




