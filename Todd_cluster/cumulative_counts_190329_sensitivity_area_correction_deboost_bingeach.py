#import matplotlib.pyplot as plt
#import numpy as np
#import math
#import matplotlib.cm as cm
#from astropy.io import fits
#import astropy.wcs as wcs
#from scipy.interpolate import interp1d













#=====================================Load in source catalogue====================================================================================




#if cluster_list[l]=='S2CLS':
#    #template_file = open("./catalogues/S2CLS_3p5sigma.cat", "r")    #==================================================================================================HERE
#    template_file = open("../"+cluster_list[l]+"/S2CLS_3p5sigma_deboost.cat", "r")    #==================================================================================================HERE
#else:
#    #template_file = open("./catalogues/"+cluster_list[l]+"_3p5sigma_edgesourcedelete.cat", "r")    #==================================================================================================HERE
#    template_file = open("../"+cluster_list[l]+"/"+cluster_list[l]+"_3p5sigma_edgesourcedelete_deboost.cat", "r")    
template_file = open("./catalogues/edgesourcedelete_deboost_190812/"+cluster_list[l]+"_3p5sigma_deboost.cat", "r")   
#==================================================================================================HERE

lines = template_file.readlines()[1:]
#lines = template_file.readlines()
template_file.close()

response_curve=[]

for i in range(0, len(lines)):
    separated_lines=lines[i].split() 
    response_curve.append(separated_lines)


response_curve = np.array(response_curve)
#flux_mjy_pre=response_curve[:,6]      #=========================================================HERE: 7 for Todd's cluster, 6 for our clusters
#flux_mjy_err_pre=response_curve[:,7]    #=========================================================HERE: 8 for Todd's cluster, 7 for our clusters
#sn_ratio_pre=response_curve[:,4]     #=========================================================HERE: 5 for Todd's cluster, 4 for our clusters

if cluster_list[l]=='S2CLS':
    flux_mjy_pre=response_curve[:,8]      #=========================================================HERE: 7 for Todd's cluster, 6 for our clusters
    flux_mjy_err_pre=response_curve[:,10]    #=========================================================HERE: 8 for Todd's cluster, 7 for our clusters
    sn_ratio_pre=response_curve[:,11]     #=========================================================HERE: 5 for Todd's cluster, 4 for our clusters
else:
    flux_mjy_pre=response_curve[:,9]      #=========================================================HERE: 7 for Todd's cluster, 6 for our clusters
    flux_mjy_err_pre=response_curve[:,14]    #=========================================================HERE: 8 for Todd's cluster, 7 for our clusters #ncounts analysis doesn't use this error so keep it as the one including confusion noise.
    #sn_ratio_pre=response_curve[:,12]     #=========================================================HERE: 5 for Todd's cluster, 4 for our clusters: this one adds confusion in quadrature, but actually it's already included in the variance map
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


flux_mjy_sorted_pre=[]
flux_mjy_err_sorted=[]
sn_ratio_sorted=[]


if len(flux_mjy_sn) > 0:
    flux_mjy_sorted_pre, flux_mjy_err_sorted, sn_ratio_sorted = (list(t) for t in zip(*sorted(zip(flux_mjy_sn, flux_mjy_err, sn_ratio))))

flux_mjy_sorted_pre=flux_mjy_sorted_pre[::-1]
flux_mjy_err_sorted=flux_mjy_err_sorted[::-1]
sn_ratio_sorted=sn_ratio_sorted[::-1]

#print 'the faintest source has flux (mJy) of ', flux_mjy_sorted_pre[-1]
#print flux_mjy_sorted_pre
#print flux_mjy_err_sorted
#print sn_ratio_sorted


flux_mjy_sorted=[]
for k in range(0, len(flux_mjy_sorted_pre)):

    #flux_mjy_sorted.append(flux_mjy_sorted_pre[k]/(1.0+0.2*pow(sn_ratio_sorted[k]/5.0,-2.3)))
    flux_mjy_sorted.append(flux_mjy_sorted_pre[k])  #already deboost in cat


flux_mjy_sorted = np.array(flux_mjy_sorted)
flux_mjy_sorted=np.array([float(i) for i in flux_mjy_sorted])  #(wavelength in angstrom)


#print 'xxxxxxxxxxxxxxxxxxxxx'
#print flux_mjy_sorted








#===========================produce noise/sensitivity map=============================================================

if cluster_list[l]=='S2CLS':
    hdulist_noise = fits.open('../S2CLS/S2CLS_rms_170919.fits')      #================================================================================HERE Orz
    w_noise = wcs.WCS(hdulist_noise[0].header, hdulist_noise)
    NAXIS1_noise=hdulist_noise[0].header['NAXIS1']
    NAXIS2_noise=hdulist_noise[0].header['NAXIS2']
    CDELT1_noise=hdulist_noise[0].header['CDELT1']
    CDELT2_noise=hdulist_noise[0].header['CDELT2']
    scidata_noise = hdulist_noise[0].data  #================variance: unit of pW^2
    hdulist_noise.close()

else:
    hdulist_noise = fits.open('./'+cluster_list[l]+'/mf_crop.fits')      #================================================================================HERE
    w_noise = wcs.WCS(hdulist_noise[1].header, hdulist_noise)
    NAXIS1_noise=hdulist_noise[1].header['NAXIS1']
    NAXIS2_noise=hdulist_noise[1].header['NAXIS2']
    CDELT1_noise=hdulist_noise[1].header['CDELT1']
    CDELT2_noise=hdulist_noise[1].header['CDELT2']
    scidata_noise = hdulist_noise[1].data  #================variance: unit of pW^2
    hdulist_noise.close()

#print CDELT1_noise
#print CDELT2_noise







#print scidata_noise[0,50,:]


noise_value_array_tmp=[]

for i in range(0, NAXIS2_noise):
    for j in range(0, NAXIS1_noise):
        if cluster_list[l]=='S2CLS':
            noise_value_array_tmp.append(math.sqrt(scidata_noise[i,j]))   #Orz ([i,j] for S2CLS; [0,i,j] for ordinary fields)
        else:
            noise_value_array_tmp.append(math.sqrt(scidata_noise[0,i,j]))   #Orz ([i,j] for S2CLS; [0,i,j] for ordinary fields)

noise_value_array_tmp = np.array(noise_value_array_tmp)
noise_value_array_tmp=np.array([float(i) for i in noise_value_array_tmp])  #(wavelength in angstrom)

noise_value_array_tmp = noise_value_array_tmp[~np.isnan(noise_value_array_tmp)]

#if cluster_list[l]=='S2CLS':
#    noise_value_array_tmp=noise_value_array_tmp*1.0*1.0  #--> for 2SCLS
#else:
#    noise_value_array_tmp=noise_value_array_tmp*537.0*1000.0  #--> for ordinary fields   #convert unit form pW to mJy #Orz

noise_value_array_tmp=noise_value_array_tmp*1.0*1.0  #--> for Todd's clusters

noise_value_array=[]
for i in range(0, len(noise_value_array_tmp)):
    noise_value_array.append(math.sqrt(pow(noise_value_array_tmp[i],2.0)+pow(0.7,2.0)))
    #noise_value_array.append(noise_value_array_tmp[i])
    

noise_value_array_sorted=[]
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
    noise_area.append(noise_pix_counter*abs(CDELT1_noise)*abs(CDELT2_noise))    #area in deg^2


#print noise_area[0:10]
#print 'the total cumulative area (deg^2) for sensitivity pixels is', noise_area[-1]





area_source=[]
for k in range(0, len(flux_mjy_sorted_pre)):

    area_source.append(np.interp(flux_mjy_sorted_pre[k]*(1.0+0.2*pow(sn_ratio_sorted[k]/5.0,-2.3)), noise_value_array_sorted, noise_area)) #Have to use flux "without" boosted, so boost back here!

area_source = np.array(area_source)
area_source=np.array([float(i) for i in area_source])  #(wavelength in angstrom)


print 'xxxxxxxxxxxxxxxxxxxxx'
#print flux_mjy_sorted
print area_source
#print "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"








#plt.plot(noise_value_array_sorted, noise_area, linewidth=2.0)
#plt.tick_params(width=2, length=16, which='major')
#plt.tick_params(width=2, length=5, which='minor')
#plt.xscale('log',nonposy='clip')
#plt.yscale('log',nonposy='clip')
#plt.grid()
#plt.xlabel('sensitivity (mJy)')  
#plt.title(cluster_list[l])  #===================================================================================================HERE
#plt.ylabel('Cumulative Area (deg^2)')  
#plt.rc('font', size=30)
#plt.show()









#==================================import fits image to calculate the map area (not used)================================


#hdulist = fits.open('./NGP9_snr_170919.fits')      #================================================================================HERE
#w = wcs.WCS(hdulist[0].header, hdulist)
#NAXIS1=hdulist[0].header['NAXIS1']
#NAXIS2=hdulist[0].header['NAXIS2']
#CDELT1=hdulist[0].header['CDELT1']
#CDELT2=hdulist[0].header['CDELT2']
#scidata = hdulist[0].data
#hdulist.close()


#pix_value_array=[]

#for i in range(0, NAXIS2):
#    for j in range(0, NAXIS1):
#        pix_value_array.append(scidata[0,i,j])

#pix_value_array = np.array(pix_value_array)
#pix_value_array=np.array([float(i) for i in pix_value_array])  #(wavelength in angstrom)

#pix_value_array = pix_value_array[~np.isnan(pix_value_array)]

##print len(pix_value_array), 'pixels'

#area_map=len(pix_value_array)*abs(CDELT1)*abs(CDELT2)

#print 'the total area from flux map (deg^2) is', area_map, 'degree^2'



##n_cum_perdeg=n_cum/area_map  #cumulative number per deg^2


##print n_cum_perdeg





#================================Import errorbars from completeness================================================================

template_file = open("./"+cluster_list[l]+"/"+cluster_list[l]+"_completeness_level_origflux.cat", "r") 

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
#binsize=2  #=================================================================================HERE: bin size 1 or 2?

bin_flux=[]
bin_number=[]
bin_number_err=[]
correct_area=[]

if binsize==1:
    #geach_flux=[3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5]  #mJy
    geach_flux_widebin=np.arange(100.5, 2.5, -1.0)  #====Here to change binwidth

    flux_widebin_low=100.0  #====Here to change binwidth
    flux_widebin_high=101.0  #====Here to change binwidth


if binsize==2:
    #geach_flux=[3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5]  #mJy
    geach_flux_widebin=np.arange(100.0, 2.0, -2.0)  #====Here to change binwidth

    flux_widebin_low=99.0  #====Here to change binwidth
    flux_widebin_high=101.0  #====Here to change binwidth



for j in range(0, len(geach_flux_widebin)):
    counter=0.0
    counter_per_area=0.0
    counter_per_area_err=0.0
    #print flux_widebin_low
    for i in range(0, len(flux_mjy_sorted)):
        if flux_mjy_sorted[i] >= flux_widebin_low and flux_mjy_sorted[i] < flux_widebin_high:
            counter=counter+1.0
            per_area=1.0/area_source[i]
            counter_per_area=counter_per_area + per_area
            counter_per_area_err=counter_per_area_err + pow(math.sqrt(1.0)*per_area, 2.0)
    if counter > 0.0:
        bin_flux.append(geach_flux_widebin[j])
        bin_number.append(counter_per_area)
        bin_number_err.append(math.sqrt(counter_per_area_err))
    if binsize==1:
        flux_widebin_low=flux_widebin_low-1.0  #====Here to change binwidth
        flux_widebin_high=flux_widebin_high-1.0  #====Here to change binwidth
    if binsize==2:
        flux_widebin_low=flux_widebin_low-2.0  #====Here to change binwidth
        flux_widebin_high=flux_widebin_high-2.0  #====Here to change binwidth

#remember to change the scatter x and err when plotting in the end!
#print bin_number
#print bin_number_err
#print bin_flux
#print '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'


#------------------------------------------------------------------------




##area_per_source=[]
#n_per_source_per_area=[]
#n_per_source_per_area_err=[]
#for i in range(0, len(flux_mjy_sorted)):
#    ##area_per_source.append(np.interp(flux_mjy_sorted[i], noise_value_array_sorted, noise_area))    
#    n_per_source_per_area.append(1.0/np.interp(flux_mjy_sorted[i], noise_value_array_sorted, noise_area))  #---sensitivity area ===============================================HERE: total area or sensitivity area?
#    n_per_source_per_area_err.append(math.sqrt(i+1.0)/np.interp(flux_mjy_sorted[i], noise_value_array_sorted, noise_area))   #---sensitivity area
#    #n_per_source_per_area.append(1.0/noise_area[-1])   #-------------------total area
#    #n_per_source_per_area_err.append(math.sqrt(i+1.0)/noise_area[-1])   #-------------------total area
 

##print n_per_source_per_area
##print n_per_source_per_area_err
##print 'eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee'





n_cum=[]
n_cum_err=[]
completeness_err=[]

for i in range(0, len(bin_flux)):
    counter=0
    add_amount=0.0
    add_amount_err=0.0
    while counter <= i:
        add_amount=add_amount+bin_number[counter]
        add_amount_err=add_amount_err+pow(bin_number_err[counter],2.0)
        counter=counter+1
    #n_cum.append(add_amount/noise_area[-1])  #-------------------total area
    #n_cum_err.append(math.sqrt(add_amount)/noise_area[-1])  #-------------------total area
    n_cum.append(add_amount*noise_area[-1])  #---sensitivity area ===============================================HERE: total area or sensitivity area?
    n_cum_err.append(math.sqrt(add_amount_err)*noise_area[-1])  #---sensitivity area

    for db in range(0, len(flux_deboosted)): 
        if flux_deboosted[db] == bin_flux[i]:
            completeness_err.append(add_amount*noise_area[-1]*completeness_level_err[db]/pow(completeness_level[db], 2.0))

n_cum = np.array(n_cum)
n_cum=np.array([float(i) for i in n_cum])  #(wavelength in angstrom)


print n_cum
print n_cum_err
#print completeness_err

print 'hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh'



#------------------try binning------------------------------

#bin_number=[]
#counter=0.0
#for i in range(0, len(flux_mjy_sorted)):
#    if flux_mjy_sorted[i] >= 5.0 and flux_mjy_sorted[i] < 6.0:
#        counter=counter+1.0
#    bin_number.append(counter)

#print bin_number



#-----------------------------------------------------------
















#n_cum=[]
#n_cum_err=[]
#counter_cum=0.0
#for i in range(0, len(flux_mjy_sorted)):

#    counter_cum=counter_cum+1.0
#    n_cum.append(counter_cum/np.interp(flux_mjy_sorted[i], noise_value_array_sorted, noise_area))
#    n_cum_err.append(math.sqrt(counter_cum/np.interp(flux_mjy_sorted[i], noise_value_array_sorted, noise_area)))
#    #n_cum.append(counter_cum/area_map)
#    #n_cum_err.append(math.sqrt(counter_cum/area_map))





#------test area: testfor binning like Geach---------------------------------------------
#print 'test starts'

#test_flux=[14.5,13.5,12.5,11.5,10.5,9.5,8.5,7.5,6.5,5.5,4.5]
#test_area=[]
#for i in range(0, len(test_flux)):
#    test_area.append(np.interp(test_flux[i], noise_value_array_sorted, noise_area))
#print test_area

#test_n=[1.0,0.0,0.0,4.0,4.0,3.0,2.0,5.0,6.0,3.0,1.0]

#test_n_area=[]
#for i in range(0, len(test_n)):
#    test_n_area.append(test_n[i]/test_area[i])

#print test_n_area

##test_n_cum_area=[]
##for i in range(0, len(test_n_area)):
##    counter=0
##    add_amount=0.0
##    while counter <= i:
##        add_amount=add_amount+test_n_area[counter]
##        counter=counter+1
##    test_n_cum_area.append(add_amount)

##print test_n_cum_area

#test_n_cum=[1.0, 1.0, 1.0, 5.0, 9.0, 12.0, 14.0, 19.0, 25.0, 28.0, 29.0]

#test_n_cum_area=[]
#for i in range(0, len(test_area)):
#    test_n_cum_area.append(test_n_cum[i]/test_area[i])

#print test_flux
#print test_n_cum_area

#print 'test ends'
#----------------------------------------------------------------------------------------

#flux_mjy_sorted = np.array(flux_mjy_sorted)
#flux_mjy_sorted=np.array([float(i) for i in flux_mjy_sorted])  #(wavelength in angstrom)


#n_cum = np.array(n_cum)
#n_cum=np.array([float(i) for i in n_cum])  #(wavelength in angstrom)
#n_cum_err = np.array(n_cum_err)
#n_cum_err=np.array([float(i) for i in n_cum_err])  #(wavelength in angstrom)

##n_cum = n_cum[::-1]
##n_cum_err = n_cum_err[::-1]

#print n_cum
#print n_cum_err










#====================================================Plotting======================================================================

#manual_flux=[11.38795, 10.661, 9.514, 9.39545, 8.75898, 7.7187, 7.60141, 7.15786, 7.12553, 6.84517, 6.75961, 6.58465, 6.5636, 6.21689, 6.01243, 5.909, 5.89506, 5.03845, 4.93665]
#manual_cum_n=[71.6197243914, 143.239448783, 214.859173174, 286.478897565, 358.098621957, 429.718346348, 501.338070739, 572.957795131, 644.577519522, 716.197243914, 787.816968305, 859.436692696, 931.056417088, 1002.67614148, 1074.29586587, 1145.91559026, 1217.53531465, 1289.15503904, 1360.77476344]

#plt.plot(manual_flux, manual_cum_n, color='green', label='NGP9')  #plots for manually putting the numbers in
#plt.scatter(manual_flux, manual_cum_n, color='green')  #plots for manually putting the numbers in


#plt.plot(test_flux, test_n_cum_area, color='green')  #plots for binning like Geach
#plt.scatter(test_flux, test_n_cum_area, color='green') #plots for binning like Geach


##plt.errorbar(flux_mjy_sorted, n_cum, yerr=n_cum_err, color='blue',fmt='o', capsize=5, elinewidth=2,markeredgewidth=2)
#plt.errorbar(flux_mjy_sorted, n_cum, yerr=n_per_source_per_area_err, color='blue',fmt='o', capsize=5, elinewidth=2,markeredgewidth=2)
#plt.plot(flux_mjy_sorted, n_cum, label='This work')
#plt.scatter(flux_mjy_sorted, n_cum)

#ax = figa.add_subplot(3,5,l+1)

if 0<=l<=4:
    ax=ax_all[0,l]
elif 5<=l<=9:
    ax=ax_all[1,l-5]
else:
    ax=ax_all[2,l-10]


geach_cum_number_plot= [rr * noise_area[-1] for rr in geach_cum_number]
geach_cum_number_err_up_plot = [rr2 * noise_area[-1] for rr2 in geach_cum_number_err_up]
geach_cum_number_err_low_plot = [rr3 * noise_area[-1] for rr3 in geach_cum_number_err_low]
asymmetric_error_plot = [geach_cum_number_err_low_plot, geach_cum_number_err_up_plot]

print geach_cum_number_plot
print noise_area[-1]
print '222222222222222222222222222'

ax.plot(geach_flux, geach_cum_number_plot, color='red')
ax.scatter(geach_flux, geach_cum_number_plot, color='red')
ax.errorbar(geach_flux, geach_cum_number_plot, yerr=asymmetric_error_plot, color='red',fmt='o', capsize=5, elinewidth=2,markeredgewidth=2)

print bin_flux
print flux_deboosted
print len(n_cum)
print len(completeness_err)


ax.errorbar(bin_flux, n_cum, yerr=completeness_err, color='blue',fmt='o', capsize=5, elinewidth=2,markeredgewidth=2)
ax.plot(bin_flux, n_cum, label=cluster_list[l])
ax.scatter(bin_flux, n_cum, s=80.0, marker='s')





ax.tick_params(width=2, length=16, which='major')
ax.tick_params(width=2, length=5, which='minor')


ax.set_xscale('log',nonposy='clip')
ax.set_yscale('log',nonposy='clip')

ax.set_xlim((3.0, 30.0))
ax.set_ylim((1e-2, 1e2))


leg = ax.legend(handlelength=0, handletextpad=0, fancybox=True, loc='upper right')
#print leg
#print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
for item in leg.legendHandles:
    item.set_visible(False)
leg.get_frame().set_linewidth(0.0)
leg.get_frame().set_alpha(0.0)


#plt.legend(loc=3)
ax.grid()


if 10<=l<=14:
    ax.set_xlabel('Flux Density (mJy)')  
#ax.title(cluster_list[l])  #===============================================================================================================HERE
if l==0 or l==5 or l==10:
    ax.set_ylabel(r'Cumulative Count' +'\n'+'(per map size)')  
plt.rc('font', size=18)
#plt.show()



mock_catalogue=np.column_stack((bin_flux, n_cum, completeness_err)) 
#mock_catalogue=np.column_stack((bin_flux, n_cum))

mock_catalogue_geach=np.column_stack((geach_flux, geach_cum_number_plot, geach_cum_number_err_low_plot, geach_cum_number_err_up_plot)) 

if binsize==1:
    np.savetxt('./ncounts/cum_190813/cum_counts_'+cluster_list[l]+'_bin1.cat', mock_catalogue, delimiter=' ', header='bin_flux n_cum completeness_err')
    np.savetxt('./ncounts/cum_190813/cum_counts_'+cluster_list[l]+'_bin1.cat', mock_catalogue, delimiter=' ', header='bin_flux n_cum completeness_err')
if binsize==2:
    np.savetxt('./ncounts/cum_190813/cum_counts_'+cluster_list[l]+'_bin2.cat', mock_catalogue, delimiter=' ', header='bin_flux n_cum completeness_err')
    np.savetxt('./ncounts/cum_190813/cum_counts_geach_bin2.cat', mock_catalogue_geach, delimiter=' ', header='geach_flux geach_cum_number_plot geach_cum_number_err_low_plot geach_cum_number_err_up_plot')


