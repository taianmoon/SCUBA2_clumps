import matplotlib.pyplot as plt
import numpy as np
#import math
#import matplotlib.cm as cm
#from astropy.io import fits
#import astropy.wcs as wcs
#from scipy.interpolate import interp1d
#import pyfits
#import os
#import subprocess
#import scipy



#cluster_list=['Bootes1', 'EGS', 'G12', 'Lockman']
#cluster_list=['EGS', 'Lockman']
#cluster_list=['G12', 'NGP1', 'NGP2', 'NGP3', 'NGP4', 'NGP5', 'NGP6', 'NGP7', 'NGP8', 'NGP9']
cluster_list=['Bootes1', 'EGS', 'Lockman', 'G12', 'NGP1', 'NGP2', 'NGP3', 'NGP4', 'NGP5', 'NGP6', 'NGP7', 'NGP8', 'NGP9', 'PCL1002', 'S2CLS_random']



#zx_list=['zx1', 'zx2', 'zx3', 'zx4', 'zx5', 'zx6', 'zx7', 'zx8', 'zx9', 'zx10', 'zx11', 'zx12', 'zx13', ]
#ax_list=['ax1', 'ax2', 'ax3', 'ax4', 'ax5', 'ax6', 'ax7', 'ax8', 'ax9', 'ax10', 'ax11', 'ax12', 'ax13', ]

figz = plt.figure(figsize=(16.0, 10.0))
#zx1 = figz.add_subplot(3,5,1)
#zx2 = figz.add_subplot(3,5,2)
#zx3 = figz.add_subplot(3,5,3)
#zx4 = figz.add_subplot(3,5,4)
#zx5 = figz.add_subplot(3,5,5)
#zx6 = figz.add_subplot(3,5,6)
#zx7 = figz.add_subplot(3,5,7)
#zx8 = figz.add_subplot(3,5,8)
#zx9 = figz.add_subplot(3,5,9)
#zx10 = figz.add_subplot(3,5,10)
#zx11 = figz.add_subplot(3,5,11)
#zx12 = figz.add_subplot(3,5,12)
#zx13 = figz.add_subplot(3,5,13)

figa = plt.figure(figsize=(16.0, 10.0))
#ax1 = figa.add_subplot(3,5,1)
#ax2 = figa.add_subplot(3,5,2)
#ax3 = figa.add_subplot(3,5,3)
#ax4 = figa.add_subplot(3,5,4)
#ax5 = figa.add_subplot(3,5,5)
#ax6 = figa.add_subplot(3,5,6)
#ax7 = figa.add_subplot(3,5,7)
#ax8 = figa.add_subplot(3,5,8)
#ax9 = figa.add_subplot(3,5,9)
#ax10 = figa.add_subplot(3,5,10)
#ax11 = figa.add_subplot(3,5,11)
#ax12 = figa.add_subplot(3,5,12)
#ax13 = figa.add_subplot(3,5,13)


#leg = plt.legend(handlelength=0, handletextpad=0, fancybox=True, loc='lower left')
#for item in leg.legendHandles:
#    item.set_visible(False)
#leg.get_frame().set_linewidth(0.0)
#leg.get_frame().set_alpha(0.75)


for l in range(0, len(cluster_list)):

    print cluster_list[l]
    execfile("./photoz_err.py")



plt.show()
