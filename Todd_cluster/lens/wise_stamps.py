#DO FIRST
#import os
#os.system('export PATH=$PATH:$HOME/Downloads/montage/bin')

import pyfits
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
#import montage_wrapper as montage
from astropy.io import fits
import matplotlib.cm as cm
import aplpy
import glob


##========================================================HERE============================================
#
#input_fits_image_band1='./input/bootes1_ks_2h_20160228_astro_2MASS_0p248rms.fits'  #=====Ks
#input_fits_image_band2='./input/bootes1_j_1h12min_20160226_astro_2MASS_0p182rms.fits'  #=====J
#input_fits_image_band3='./input/bootes1_h_42min_20160226_astro_2MASS_0p220rms.fits'    #=====H
##input_fits_image_band4='./input/gama12_i_band_ACAM_53min_20160226_astro_0p241rms.fits'  #=====I
##input_fits_image_band5='./input/G12v230_IRAC_Mosaic_36.fits'                            #=====3.6 micron
##input_fits_image_band6='./input/G12v230_IRAC_Mosaic_45.fits'                            #=====4.5 micron
#
band4_present=1  #=================================================================1:yes 0:no (I band)
band5_present=0  #=================================================================1:yes 0:no (3.6 band)
band6_present=0  #=================================================================1:yes 0:no (4.5 band)
#
#
#================================Load in fits table===============================================

#names_cutoff = pyfits.open('../lens_candidates_cat.fits')                            #========================================================HERE
#names_cutoff_data = names_cutoff[1].data
#names_new= names_cutoff_data.field('names_new')
#Y_WORLD_K= names_cutoff_data.field('Y_WORLD_K')
#K_mag = names_cutoff_data.field('M_calibrated_K=MAG_APER_K+26.2715655err0.0033601851') #==================================================HERE note to correct column name
#J_mag = names_cutoff_data.field('M_calibrated_J=MAG_APER_J+27.93667err0.0224088083') #==================================================HERE note to correct column name
#H_mag = names_cutoff_data.field('M_calibrated_H=MAG_APER_H+27.467062err0.0230624783') #==================================================HERE note to correct column name
##I_mag = names_cutoff_data.field('M_calibrated_I=MAG_APER_I+32.5366530465err0.0029509302') #==================================================HERE note to correct column name
#NUMBER_K= names_cutoff_data.field('NUMBER_K')
#
#

file_array_W1 = glob.glob("./WISE_images_cluster/W1/*")  #===============HERE
file_array_W2 = glob.glob("./WISE_images_cluster/W2/*")  #===============HERE
file_array_W3 = glob.glob("./WISE_images_cluster/W3/*")  #===============HERE
file_array_W4 = glob.glob("./WISE_images_cluster/W4/*")  #===============HERE


print len(file_array_W1)


##===============================make a webpage showing the cutoff images============================================
#
#
#
#f=open('./cutoff_images.html','w+')
#
#f.write('<!DOCTYPE html>\n')
#f.write('<html>\n')
#f.write('<head>\n')
#f.write('\t<title>cutoff_images_wise</title>\n')
#f.write('</head>\n')
#f.write('<body>\n')
#
#f.write('<table border="1">')
#
#f.write('<tr>')
#
#
#
#f.write('<td>W1</td>')
#f.write('<td>W2</td>')
#f.write('<td>W3</td>')
#if band4_present==1:
#    f.write('<td>W4</td>')
#
#if band5_present==1:
#    f.write('<td>3.6 micron image</td>')
#if band6_present==1:
#    f.write('<td>4.5 micron image</td>')
#
#f.write('</tr>')
#
##==================to be continued in the loop=========================











#==============================do cutoff==========================================================

#for counter_cutoff in range(len(file_array_W1)):
#for counter_cutoff in range(0,10):
#for counter_cutoff in range(10,20):
#for counter_cutoff in range(20,30):
#for counter_cutoff in range(30,40):
for counter_cutoff in range(40,46):

    #try:
    #    montage.mSubimage(input_fits_image_band1, './output/cutoff_ks_'+str(NUMBER_K[counter_cutoff])+'.fits', ra=X_WORLD_K[counter_cutoff], dec=Y_WORLD_K[counter_cutoff], xsize=0.0023)
    #except:
    #    print '==problem in ks image=='
    #try:
    #    montage.mSubimage(input_fits_image_band2, './output/cutoff_j_'+str(NUMBER_K[counter_cutoff])+'.fits', ra=X_WORLD_K[counter_cutoff], dec=Y_WORLD_K[counter_cutoff], xsize=0.0023)
    #except:
    #    print '==problem in j image=='
    #try:
    #    montage.mSubimage(input_fits_image_band3, './output/cutoff_h_'+str(NUMBER_K[counter_cutoff])+'.fits', ra=X_WORLD_K[counter_cutoff], dec=Y_WORLD_K[counter_cutoff], xsize=0.0023)
    #except:
    #    print '==problem in h image=='    
    #if band4_present==1:
    #    try:
    #        montage.mSubimage(input_fits_image_band4, './output/cutoff_i_'+str(NUMBER_K[counter_cutoff])+'.fits', ra=X_WORLD_K[counter_cutoff], dec=Y_WORLD_K[counter_cutoff], xsize=0.0023)
    #    except:
    #        print '==problem in i image=='
    #if band5_present==1:
    #    try:
    #        montage.mSubimage(input_fits_image_band5, './output/cutoff_3p6_'+str(NUMBER_K[counter_cutoff])+'.fits', ra=X_WORLD_K[counter_cutoff], dec=Y_WORLD_K[counter_cutoff], xsize=0.0023)
    #    except:
    #        print '==problem in 3.6 micron image=='
    #if band6_present==1:
    #    try:    
    #        montage.mSubimage(input_fits_image_band6, './output/cutoff_4p5_'+str(NUMBER_K[counter_cutoff])+'.fits', ra=X_WORLD_K[counter_cutoff], dec=Y_WORLD_K[counter_cutoff], xsize=0.0023)
    #    except:
    #        print '==problem in 4.5 micron image=='





    try:
        cutoff_1 = aplpy.FITSFigure(file_array_W1[counter_cutoff])
        cutoff_1.show_colorscale(stretch='log')
        #cutoff_1.show_markers(xw=[X_WORLD_K[counter_cutoff]], yw=[Y_WORLD_K[counter_cutoff]], marker='+', c='w')
        cutoff_1.add_colorbar()
        cutoff_1.add_scalebar(1.0/60.0)
        cutoff_1.scalebar.set_label("1'")
        cutoff_1.scalebar.set_color('white')
        #cutoff_1.add_beam()
        #cutoff_1.beam.set_color('white')
        #cutoff_1.beam.set_hatch('+')
        cutoff_1.set_title(file_array_W1[counter_cutoff])
        cutoff_1.save('./WISE_images_cluster/output/cutoff_W1_'+str(counter_cutoff)+'.jpg')
    except:
        print '============='
    
    try:
        cutoff_2 = aplpy.FITSFigure(file_array_W2[counter_cutoff])
        cutoff_2.show_colorscale(stretch='log')
        #cutoff_2.show_markers(xw=[X_WORLD_K[counter_cutoff]], yw=[Y_WORLD_K[counter_cutoff]], marker='+', c='w')
        cutoff_2.add_colorbar()
        cutoff_2.add_scalebar(1.0/60.0)
        cutoff_2.scalebar.set_label("1'")
        cutoff_2.scalebar.set_color('white')
        #cutoff_1.add_beam()
        #cutoff_1.beam.set_color('white')
        #cutoff_1.beam.set_hatch('+')
        cutoff_2.set_title(file_array_W2[counter_cutoff])
        cutoff_2.save('./WISE_images_cluster/output/cutoff_W2_'+str(counter_cutoff)+'.jpg')
    except:
        print '============='

    try:
        print 'doing W3 now'
        cutoff_3 = aplpy.FITSFigure(file_array_W3[counter_cutoff])
        cutoff_3.show_colorscale(stretch='log')
        #cutoff_3.show_markers(xw=[X_WORLD_K[counter_cutoff]], yw=[Y_WORLD_K[counter_cutoff]], marker='+', c='w')
        cutoff_3.add_colorbar()
        cutoff_3.add_scalebar(1.0/60.0)
        cutoff_3.scalebar.set_label("1'")
        cutoff_3.scalebar.set_color('white')
        #cutoff_1.add_beam()
        #cutoff_1.beam.set_color('white')
        #cutoff_1.beam.set_hatch('+')
        cutoff_3.set_title(file_array_W3[counter_cutoff])
        cutoff_3.save('./WISE_images_cluster/output/cutoff_W3_'+str(counter_cutoff)+'.jpg')
    except:
        print '============='

    try:
        if band4_present==1:
            cutoff_4 = aplpy.FITSFigure(file_array_W4[counter_cutoff])
            cutoff_4.show_colorscale(stretch='log')
            #cutoff_4.show_markers(xw=[X_WORLD_K[counter_cutoff]], yw=[Y_WORLD_K[counter_cutoff]], marker='+', c='w')
            cutoff_4.add_colorbar()
            cutoff_4.add_scalebar(1.0/60.0)
            cutoff_4.scalebar.set_label("1'")
            cutoff_4.scalebar.set_color('white')
            #cutoff_1.add_beam()
            #cutoff_1.beam.set_color('white')
            #cutoff_1.beam.set_hatch('+')
            cutoff_4.set_title(file_array_W4[counter_cutoff])
            cutoff_4.save('./WISE_images_cluster/output/cutoff_W4_'+str(counter_cutoff)+'.jpg')
    except:
        print '============='

    try:
        if band5_present==1:
            cutoff_5 = aplpy.FITSFigure('./output/cutoff_3p6_'+str(NUMBER_K[counter_cutoff])+'.fits')
            cutoff_5.show_colorscale(cmap=cm.gist_heat, smooth=1)
            cutoff_5.show_markers(xw=[X_WORLD_K[counter_cutoff]], yw=[Y_WORLD_K[counter_cutoff]], marker='+', c='w')
            cutoff_5.add_colorbar()
            cutoff_5.add_scalebar(1.0/1800.0)
            cutoff_5.scalebar.set_label('2"')
            cutoff_5.scalebar.set_color('white')
            #cutoff_1.add_beam()
            #cutoff_1.beam.set_color('white')
            #cutoff_1.beam.set_hatch('+')
            cutoff_5.set_title('#'+str(NUMBER_K[counter_cutoff])+';  RA='+str(X_WORLD_K[counter_cutoff])+';  Dec='+str(Y_WORLD_K[counter_cutoff])+';  3.6 micron')
            cutoff_5.save('./output/cutoff_3p6_'+str(NUMBER_K[counter_cutoff])+'.jpg')
    except:
        print '============='

    try:
        if band6_present==1:
            cutoff_6 = aplpy.FITSFigure('./output/cutoff_4p5_'+str(NUMBER_K[counter_cutoff])+'.fits')
            cutoff_6.show_colorscale(cmap=cm.gist_heat, smooth=1)
            cutoff_6.show_markers(xw=[X_WORLD_K[counter_cutoff]], yw=[Y_WORLD_K[counter_cutoff]], marker='+', c='w')
            cutoff_6.add_colorbar()
            cutoff_6.add_scalebar(1.0/1800.0)
            cutoff_6.scalebar.set_label('2"')
            cutoff_6.scalebar.set_color('white')
            #cutoff_1.add_beam()
            #cutoff_1.beam.set_color('white')
            #cutoff_1.beam.set_hatch('+')
            cutoff_6.set_title('#'+str(NUMBER_K[counter_cutoff])+';  RA='+str(X_WORLD_K[counter_cutoff])+';  Dec='+str(Y_WORLD_K[counter_cutoff])+';  4.5 micron')
            cutoff_6.save('./output/cutoff_4p5_'+str(NUMBER_K[counter_cutoff])+'.jpg')
    except:
        print '============='





#    #===============================make a webpage showing the cutoff images (continued)============================================
#
#
#
#
#
#    f.write('<tr>')
#
#
#    f.write('<td>')
#    f.write('<a href="./output/cutoff_W1_'+str(counter_cutoff)+'.jpg" target="_blank"><img src="./output/cutoff_W1_'+str(counter_cutoff)+'.jpg" alt="./output/cutoff_W1_'+str(counter_cutoff)+'.jpg" width="600" height="450"></a>')
#    f.write('</td>')
#
#    f.write('<td>')
#    f.write('<a href="./output/cutoff_W2_'+str(counter_cutoff)+'.jpg" target="_blank"><img src="./output/cutoff_W2_'+str(counter_cutoff)+'.jpg" alt="./output/cutoff_W2_'+str(counter_cutoff)+'.jpg" width="600" height="450"></a>')
#    f.write('</td>')
#
#    f.write('<td>')
#    f.write('<a href="./output/cutoff_W3_'+str(counter_cutoff)+'.jpg" target="_blank"><img src="./output/cutoff_W3_'+str(counter_cutoff)+'.jpg" alt="./output/cutoff_W3_'+str(counter_cutoff)+'.jpg" width="600" height="450"></a>')
#    f.write('</td>')
#
#    if band4_present==1:
#        f.write('<td>')
#        f.write('<a href="./output/cutoff_W4_'+str(counter_cutoff)+'.jpg" target="_blank"><img src="./output/cutoff_W4_'+str(counter_cutoff)+'.jpg" alt="./output/cutoff_W4_'+str(counter_cutoff)+'.jpg" width="600" height="450"></a>')
#        f.write('</td>')
#
#    if band5_present==1:
#        f.write('<td>')
#        f.write('<a href="./output/cutoff_3p6_'+str(NUMBER_K[counter_cutoff])+'.jpg" target="_blank"><img src="./output/cutoff_3p6_'+str(NUMBER_K[counter_cutoff])+'.jpg" alt="./output/cutoff_3p6_'+str(NUMBER_K[counter_cutoff])+'.jpg" width="600" height="450"></a>')
#        f.write('</td>')
#
#    if band6_present==1:
#        f.write('<td>')
#        f.write('<a href="./output/cutoff_4p5_'+str(NUMBER_K[counter_cutoff])+'.jpg" target="_blank"><img src="./output/cutoff_4p5_'+str(NUMBER_K[counter_cutoff])+'.jpg" alt="./output/cutoff_4p5_'+str(NUMBER_K[counter_cutoff])+'.jpg" width="600" height="450"></a>')
#        f.write('</td>')
#
#
#    f.write('</tr>')
#
#
#
#
#
#
#
#f.write('</table>')
#
#f.write('</body>\n')
#f.write('</html>')
#f.close()
#
#













