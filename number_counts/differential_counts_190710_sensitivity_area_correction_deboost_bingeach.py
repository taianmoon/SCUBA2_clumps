#import matplotlib.pyplot as plt
#import numpy as np
#import math
#import matplotlib.cm as cm
#from astropy.io import fits
#import astropy.wcs as wcs
#from scipy.interpolate import interp1d


#=====================================Load in source catalogue====================================================================================

if cluster_list[l]=='S2CLS':
    #template_file = open("./catalogues/S2CLS_3p5sigma.cat", "r")    #==================================================================================================HERE
    template_file = open("../"+cluster_list[l]+"/S2CLS_3p5sigma_deboost.cat", "r")    #==================================================================================================HERE
else:
    #template_file = open("./catalogues/"+cluster_list[l]+"_3p5sigma_edgesourcedelete.cat", "r")    #==================================================================================================HERE
    template_file = open("../"+cluster_list[l]+"/"+cluster_list[l]+"_3p5sigma_edgesourcedelete_deboost.cat", "r")    #================================================================================HERE


lines = template_file.readlines()[1:]
#lines = template_file.readlines()
template_file.close()

response_curve=[]

for i in range(0, len(lines)):
    separated_lines=lines[i].split() 
    response_curve.append(separated_lines)


response_curve = np.array(response_curve)
#flux_mjy_pre=response_curve[:,6]      #=========================================================HERE: 7 for Todd's cluster, 6 for our clusters
#flux_mjy_err_pre=response_curve[:,7]      #=========================================================HERE: 8 for Todd's cluster, 7 for our clusters
#sn_ratio_pre=response_curve[:,4]          #=========================================================HERE: 5 for Todd's cluster, 4 for our clusters
#flux_mjy=response_curve[:,9]
#flux_mjy_err=response_curve[:,10]



if cluster_list[l]=='S2CLS':
    flux_mjy_pre=response_curve[:,8]      #=========================================================HERE: 7 for Todd's cluster, 6 for our clusters
    flux_mjy_err_pre=response_curve[:,10]    #=========================================================HERE: 8 for Todd's cluster, 7 for our clusters
    sn_ratio_pre=response_curve[:,11]     #=========================================================HERE: 5 for Todd's cluster, 4 for our clusters
else:
    flux_mjy_pre=response_curve[:,9]      #=========================================================HERE: 7 for Todd's cluster, 6 for our clusters
    flux_mjy_err_pre=response_curve[:,14]    #=========================================================HERE: 8 for Todd's cluster, 7 for our clusters  #confusion noise removed!!!
    sn_ratio_pre=response_curve[:,5]     #=========================================================HERE: 5 for Todd's cluster, 4 for our clusters
    sn_ratio_confusion_pre=response_curve[:,12]     #=========================================================HERE: 5 for Todd's cluster, 4 for our clusters



flux_mjy_pre = np.array(flux_mjy_pre)
flux_mjy_err_pre = np.array(flux_mjy_err_pre)
sn_ratio_pre = np.array(sn_ratio_pre)
sn_ratio_confusion_pre = np.array(sn_ratio_confusion_pre)
flux_mjy_pre=np.array([float(i) for i in flux_mjy_pre])  #(wavelength in angstrom)
flux_mjy_err_pre=np.array([float(i) for i in flux_mjy_err_pre])  #(wavelength in angstrom)
sn_ratio_pre=np.array([float(i) for i in sn_ratio_pre])  #(wavelength in angstrom)
sn_ratio_confusion_pre=np.array([float(i) for i in sn_ratio_confusion_pre])  #(wavelength in angstrom)


flux_mjy_sn=[]
flux_mjy_err=[]
sn_ratio=[]
for i in range(0, len(flux_mjy_pre)):
    if sn_ratio_confusion_pre[i] >= sn_want: 
        flux_mjy_sn.append(flux_mjy_pre[i])
        flux_mjy_err.append(flux_mjy_err_pre[i])
        sn_ratio.append(sn_ratio_pre[i])



flux_mjy=flux_mjy_sn


flux_mjy_deboost=[]
for k in range(0, len(flux_mjy_sn)):

    #flux_mjy_deboost.append(flux_mjy_sn[k]/(1.0+0.2*pow(sn_ratio[k]/5.0,-2.3)))
    flux_mjy_deboost.append(flux_mjy_sn[k])  #already deboost in cat



flux_mjy = np.array(flux_mjy)
flux_mjy=np.array([float(i) for i in flux_mjy])  #(wavelength in angstrom)
flux_mjy_deboost = np.array(flux_mjy_deboost)
flux_mjy_deboost=np.array([float(i) for i in flux_mjy_deboost])  #(wavelength in angstrom)


print 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
#print flux_mjy
print flux_mjy_deboost
#print flux_mjy_err
#print sn_ratio
#print 'DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD'







#=========================Estimate the map size in deg^2 (using variance map!!)===============================

if cluster_list[l]=='S2CLS':
    hdulist = fits.open('../S2CLS/S2CLS_rms_170919.fits')      #=======================================================================================================================HERE Orz
    w = wcs.WCS(hdulist[0].header, hdulist)
    NAXIS1=hdulist[0].header['NAXIS1']
    NAXIS2=hdulist[0].header['NAXIS2']
    CDELT1=hdulist[0].header['CDELT1']
    CDELT2=hdulist[0].header['CDELT2']
    #CDELT1=-5.5555555555555E-4
    #CDELT2=0.00055555555555555
    scidata = hdulist[0].data  #================variance: unit of pW^2
    hdulist.close()

else:
    hdulist = fits.open('../'+cluster_list[l]+'/'+cluster_list[l]+'_flux_170919.fits')      #==========================================================================================================HERE
    w = wcs.WCS(hdulist[1].header, hdulist)
    NAXIS1=hdulist[1].header['NAXIS1']
    NAXIS2=hdulist[1].header['NAXIS2']
    CDELT1=hdulist[1].header['CDELT1']
    CDELT2=hdulist[1].header['CDELT2']
    scidata = hdulist[1].data
    hdulist.close()







pix_value_array=[]

for i in range(0, NAXIS2):
    for j in range(0, NAXIS1):
        if cluster_list[l]=='S2CLS':
            pix_value_array.append(scidata[i,j])   #Orz  ([i,j] for S2CLS; [0,i,j] for ordinary fields)
        else:
            pix_value_array.append(scidata[0,i,j])   #Orz  ([i,j] for S2CLS; [0,i,j] for ordinary fields)

pix_value_array = np.array(pix_value_array)
pix_value_array=np.array([float(i) for i in pix_value_array])  #(wavelength in angstrom)

pix_value_array = pix_value_array[~np.isnan(pix_value_array)]

#print len(pix_value_array), 'pixels'

area_map=len(pix_value_array)*abs(CDELT1)*abs(CDELT2)

#area_correction_factor=1.0  #======================================================================================================================================HERE
#area_map=area_map*area_correction_factor


#print 'the total area from flux map (deg^2) is', area_map, 'degree^2'





#---------------doing sensitivity map------------------------------------------------

noise_value_array_tmp=[]

for i in range(0, NAXIS2):
    for j in range(0, NAXIS1):
        if cluster_list[l]=='S2CLS':
            noise_value_array_tmp.append(math.sqrt(scidata[i,j]))   #Orz ([i,j] for S2CLS; [0,i,j] for ordinary fields)
        else:
            noise_value_array_tmp.append(math.sqrt(scidata[0,i,j]))   #Orz ([i,j] for S2CLS; [0,i,j] for ordinary fields)

noise_value_array_tmp = np.array(noise_value_array_tmp)
noise_value_array_tmp=np.array([float(i) for i in noise_value_array_tmp])  #(wavelength in angstrom)

noise_value_array_tmp = noise_value_array_tmp[~np.isnan(noise_value_array_tmp)]

if cluster_list[l]=='S2CLS':
    noise_value_array_tmp=noise_value_array_tmp*1.0*1.0  #--> for 2SCLS
else:
    noise_value_array_tmp=noise_value_array_tmp*537.0*1000.0  #--> for ordinary fields   #convert unit form pW to mJy #Orz



noise_value_array=[]
for i in range(0, len(noise_value_array_tmp)):
    noise_value_array.append(math.sqrt(pow(noise_value_array_tmp[i],2.0)+pow(0.7,2.0)))




noise_value_array_sorted= sorted(noise_value_array)


noise_value_array_sorted = np.array(noise_value_array_sorted)
noise_value_array_sorted=np.array([float(i) for i in noise_value_array_sorted])  #(wavelength in angstrom)

#oise_value_array_sorted = noise_value_array_sorted[::-1]

#print noise_value_array_sorted[0:10]

noise_value_array_sorted=noise_value_array_sorted*sn_want #========================================================================HERE: detection threshold 3.5-sigma ?

#print 'the faintest sensitivity (mJy) is ', noise_value_array_sorted[0]

noise_area=[]
noise_pix_counter=0.0
for i in range(0, len(noise_value_array_sorted)):
    noise_pix_counter=noise_pix_counter+1.0
    noise_area.append(noise_pix_counter*abs(CDELT1)*abs(CDELT2))    #area in deg^2

    
#print noise_area[0:10]
#print 'the total cumulative area (deg^2) for sensitivity pixels is', noise_area[-1]





area_source=[]
for k in range(0, len(flux_mjy)):

    #area_source.append(np.interp(flux_mjy[k], noise_value_array_sorted, noise_area))
    area_source.append(np.interp(flux_mjy[k]*(1.0+0.2*pow(sn_ratio[k]/5.0,-2.3)), noise_value_array_sorted, noise_area))

area_source = np.array(area_source)
area_source=np.array([float(i) for i in area_source])  #(wavelength in angstrom)


#print area_source
#print "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"






zx = figz.add_subplot(3,5,l+1)

#zx.scatter(flux_mjy, area_source, color='red')
#zx.errorbar(flux_mjy, area_source, xerr=flux_mjy_err, color='red',fmt='o', capsize=5, elinewidth=2,markeredgewidth=2)
zx.plot(noise_value_array_sorted, noise_area, linewidth=2.0, label=cluster_list[l])
zx.tick_params(width=2, length=16, which='major')
zx.tick_params(width=2, length=5, which='minor')
zx.set_xscale('log',nonposy='clip')
zx.set_yscale('log',nonposy='clip')
zx.grid()

zx.set_xlim([4.0, 30.0])
zx.set_ylim([1e-5, 1e-1])
#zx.set_ylim(bottom=1e-5)


leg = zx.legend(handlelength=0, handletextpad=0, fancybox=True, loc='upper left')
#print leg
#print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
for item in leg.legendHandles:
    item.set_visible(False)
leg.get_frame().set_linewidth(0.0)
leg.get_frame().set_alpha(0.75)


zx.set_xlabel('Sensitivity (mJy)')  
#zx.set_title(cluster_list[l])  #===================================================================================================HERE
zx.set_ylabel(r'Cumulative Area ($deg^{2}$)')  
plt.rc('font', size=12)
#figz.savefig('./diff_sensitivity_all.eps')






#================================Import errorbars from completeness===============================================================

template_file = open("../completeness/"+cluster_list[l]+"/"+cluster_list[l]+"_completeness_level_origflux.cat", "r") 

lines = template_file.readlines()[1:]
#lines = template_file.readlines()
template_file.close()

response_curve=[]

for i in range(0, len(lines)):
    separated_lines=lines[i].split() 
    response_curve.append(separated_lines)


response_curve = np.array(response_curve)

flux_deboosted=response_curve[:,0]      #=========================================================HERE: 7 for Todd's cluster, 6 for our clusters
completeness_level=response_curve[:,1]    #=========================================================HERE: 8 for Todd's cluster, 7 for our clusters
completeness_level_err=response_curve[:,2]    #=========================================================HERE: 8 for Todd's cluster, 7 for our clusters

flux_deboosted = np.array(flux_deboosted)
completeness_level = np.array(completeness_level)
completeness_level_err = np.array(completeness_level_err)
flux_deboosted=np.array([float(i) for i in flux_deboosted])  #(wavelength in angstrom)
completeness_level=np.array([float(i) for i in completeness_level])  #(wavelength in angstrom)
completeness_level_err=np.array([float(i) for i in completeness_level_err])  #(wavelength in angstrom)









#------------------binning------------------------------

bin_number=[]
bin_number_err=[]
correct_area=[]

if binsize==1:
    geach_flux=[3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5]  #mJy
    #geach_flux_widebin=[4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0]
    #geach_flux_widebin=np.arange(4.0, 100.0, 2.0)  #=======================================================HERE: binwidth
    geach_flux_widebin=np.arange(3.5, 100.5, 1.0) #=======================================================HERE: binwidth

    #flux_bin_low=3.0
    #flux_bin_high=4.0

    flux_widebin_low=3.0  #=======================================================HERE: binwidth
    flux_widebin_high=4.0  #=======================================================HERE: binwidth


if binsize==2:
    geach_flux=[4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0]  #mJy
    geach_flux_widebin=np.arange(4.0, 100.0, 2.0) #=======================================================HERE: binwidth
    flux_widebin_low=3.0  #=======================================================HERE: binwidth
    flux_widebin_high=5.0  #=======================================================HERE: binwidth


completeness_err=[]

#for j in range(0, len(geach_flux)):
for j in range(0, len(geach_flux_widebin)):
    counter=0.0
    counter_area=0.0
    counter_per_area=0.0
    counter_per_area_err=0.0
    #print flux_bin_low
    #print flux_widebin_low
    for i in range(0, len(flux_mjy_deboost)):
        #if flux_mjy[i] >= flux_bin_low and flux_mjy[i] < flux_bin_high:
        if flux_mjy_deboost[i] >= flux_widebin_low and flux_mjy_deboost[i] < flux_widebin_high:
            counter=counter+1.0
            counter_area=counter_area + area_source[i]
            per_area=1.0/area_source[i]
            counter_per_area=counter_per_area + per_area
            counter_per_area_err=counter_per_area_err + pow(math.sqrt(1.0)*per_area, 2.0)

    if binsize==1:
        bin_number.append(counter_per_area/1.0)  #=======================================================HERE: binwidth
        bin_number_err.append(math.sqrt(counter_per_area_err)/1.0)  #=======================================================HERE: binwidth
        if counter==0.0:
            correct_area.append(0.0)
        else:
            correct_area.append(counter_area/counter)
        #flux_bin_low=flux_bin_low+1.0
        #flux_bin_high=flux_bin_high+1.0
        flux_widebin_low=flux_widebin_low+1.0  #=======================================================HERE: binwidth
        flux_widebin_high=flux_widebin_high+1.0  #=======================================================HERE: binwidth
        #flux_widebin_low=flux_widebin_low+1.0
        #flux_widebin_high=flux_widebin_high+1.0
    if binsize==2:
        bin_number.append(counter_per_area/2.0)  #=======================================================HERE: binwidth
        bin_number_err.append(math.sqrt(counter_per_area_err)/2.0)  #=======================================================HERE: binwidth
        if counter==0.0:
            correct_area.append(0.0)
        else:
            correct_area.append(counter_area/counter)
        flux_widebin_low=flux_widebin_low+2.0  #=======================================================HERE: binwidth
        flux_widebin_high=flux_widebin_high+2.0  #=======================================================HERE: binwidth

    if 4.0 <= geach_flux_widebin[j] <= 20.0:
        for db in range(0, len(flux_deboosted)):
            if geach_flux_widebin[j] == flux_deboosted[db]:
                completeness_err.append((counter_per_area/2.0)*noise_area[-1]*completeness_level_err[db]/pow(completeness_level[db], 2.0))
    else:
        completeness_err.append(0.0)


#remember to change the scatter x and err when plotting in the end!
#print bin_number
#print bin_number_err
#print correct_area
print completeness_err
print geach_flux_widebin
print "QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ"














#========================Estimate differential number per deg^2==========================

bin_number_perarea=[]
bin_number_perarea_err=[]
for i in range(0, len(bin_number)):
    if correct_area[i]==0.0:
        bin_number_perarea.append(bin_number[i])  #---sensitivity area ============================HERE: total area or sensitivity area?
        bin_number_perarea_err.append(bin_number_err[i])   #---sensitivity area
    else:
        #bin_number_perarea.append(bin_number[i]/area_map)   #-------------------total area
        #bin_number_perarea_err.append(bin_number_err[i]/area_map)   #-------------------total area
        bin_number_perarea.append(bin_number[i]*noise_area[-1])  #---sensitivity area ============================HERE: total area or sensitivity area?
        bin_number_perarea_err.append(bin_number_err[i]*noise_area[-1])   #---sensitivity area





bin_number_perarea = np.array(bin_number_perarea)
bin_number_perarea=np.array([float(i) for i in bin_number_perarea])  #(wavelength in angstrom)
bin_number_perarea_err = np.array(bin_number_perarea_err)
bin_number_perarea_err=np.array([float(i) for i in bin_number_perarea_err])  #(wavelength in angstrom)

#print bin_number_perarea
#print bin_number_perarea_err
#print "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"





#=======================plotting========================================================


if 0<=l<=4:
    ax=ax_all[0,l]
elif 5<=l<=9:
    ax=ax_all[1,l-5]
else:
    ax=ax_all[2,l-10]


#ax = figa.add_subplot(3,5,l+1)


#plt.errorbar(geach_flux, bin_number_perarea, yerr=bin_number_perarea_err, color='blue',fmt='o', capsize=5, elinewidth=2,markeredgewidth=2)
#ax.errorbar(geach_flux_widebin, bin_number_perarea, yerr=bin_number_perarea_err, color='blue',fmt='o', capsize=5, elinewidth=2,markeredgewidth=2) #removed on 180703
ax.errorbar(geach_flux_widebin, bin_number_perarea, yerr=completeness_err, color='blue',fmt='o', capsize=5, elinewidth=2,markeredgewidth=2) #removed on 180703
#plt.plot(geach_flux, bin_number_perarea, label='This work')
#plt.scatter(geach_flux, bin_number_perarea, label='This work')
ax.scatter(geach_flux_widebin, bin_number_perarea, label=cluster_list[l], s=40.0)


geach_diff_number_plot= [rr * noise_area[-1] for rr in geach_diff_number]
geach_diff_number_err_up_plot = [rr2 * noise_area[-1] for rr2 in geach_diff_number_err_up]
geach_diff_number_err_low_plot = [rr3 * noise_area[-1] for rr3 in geach_diff_number_err_low]
asymmetric_error_plot = [geach_diff_number_err_up_plot, geach_diff_number_err_low_plot]

print geach_diff_number_plot
print noise_area[-1]
print '222222222222222222222222222'





ax.plot(geach_flux, geach_diff_number_plot, color='red')
ax.scatter(geach_flux, geach_diff_number_plot, color='red')
ax.errorbar(geach_flux, geach_diff_number_plot, yerr=asymmetric_error_plot, color='red',fmt='o', capsize=5, elinewidth=2,markeredgewidth=2)


ax.tick_params(width=2, length=16, which='major')
ax.tick_params(width=2, length=5, which='minor')

#ax.set_xlim((1.0, 1e2))
ax.set_xlim((3.0, 30.0))
ax.set_ylim((1e-3, 5e1))
#ax.set_ylim((1e-1, 5e3))

ax.set_xscale('log',nonposy='clip')
ax.set_yscale('log',nonposy='clip')


leg2 = ax.legend(handlelength=0, handletextpad=0, fancybox=True, loc='upper right')
#print leg2
#print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
for item in leg2.legendHandles:
    item.set_visible(False)
leg2.get_frame().set_linewidth(0.0)
leg2.get_frame().set_alpha(0.75)


#ax.legend(loc=3)
ax.grid()

if 10<=l<=14:
    ax.set_xlabel('Flux Density (mJy)')  
#ax.set_title(cluster_list[l])  #===================================================================================================================================================================HERE
if l==0 or l==5 or l==10:
    ax.set_ylabel(r'Differential Counts' +'\n'+'(per map size)')  
plt.rc('font', size=18)
#figa.savefig('./diff_counts_all.eps')




#mock_catalogue=np.column_stack((geach_flux_widebin, bin_number_perarea, bin_number_perarea_err))
mock_catalogue=np.column_stack((geach_flux_widebin, bin_number_perarea, completeness_err))


if binsize==1:
    np.savetxt('./diff_counts_'+cluster_list[l]+'_bin1.cat', mock_catalogue, delimiter=' ') #================================================================================HERE
if binsize==2:
    np.savetxt('./diff_counts_'+cluster_list[l]+'_bin2.cat', mock_catalogue, delimiter=' ') #================================================================================HERE





