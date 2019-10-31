#import pyfits
#import matplotlib.pyplot as plt
#import numpy as np
#import math
#import glob
#import os.path


#names_wht = pyfits.open('../'+cluster_list[l]+'/'+cluster_list[l]+'_herschel_match_3p5sigma_cat_edgesourcedelete.fits')                            #================================================HERE
names_wht = pyfits.open('../'+cluster_list[l]+'/'+cluster_list[l]+'_herschel_match_3p5sigma_cat_edgesourcedelete_allscuba2source_herschelflux.fits')                            #====================================HERE
names_wht_data = names_wht[1].data

#RA= names_wht_data.field('RA_1')
#Dec= names_wht_data.field('Dec_1')






if cluster_list[l]=='Bootes1':
    error_250_pre= names_wht_data.field('et_f250')    #========================================================HERE
    error_350_pre= names_wht_data.field('et_f350')    #========================================================HERE
    error_500_pre= names_wht_data.field('et_f500')    #========================================================HERE

if cluster_list[l]=='EGS' or cluster_list[l]=='Lockman':
    error_250_pre= names_wht_data.field('et250')    #========================================================HERE
    error_350_pre= names_wht_data.field('et350')    #========================================================HERE
    error_500_pre= names_wht_data.field('et500')    #========================================================HERE

if cluster_list[l]=='G12' or cluster_list[l]=='NGP1' or cluster_list[l]=='NGP2' or cluster_list[l]=='NGP3' or cluster_list[l]=='NGP4' or cluster_list[l]=='NGP5' or cluster_list[l]=='NGP6' or cluster_list[l]=='NGP7' or cluster_list[l]=='NGP8' or cluster_list[l]=='NGP9':
    error_250_pre= names_wht_data.field('e250')    #========================================================HERE
    error_350_pre= names_wht_data.field('e350')    #========================================================HERE
    error_500_pre= names_wht_data.field('e500')    #========================================================HERE



sn_ratio_pre= names_wht_data.field('SN_ratio_1')
flux_250_pre= names_wht_data.field('F250')    #========================================================HERE(Bootes: F250, et_F250; EGS/NGP9: f250, et250; G12/NGP: F250, E250)
#error_250_pre= names_wht_data.field('et_f250')    #========================================================HER
flux_350_pre= names_wht_data.field('F350')    #========================================================HER
#error_350_pre= names_wht_data.field('et_f350')    #========================================================HER
flux_500_pre= names_wht_data.field('F500')    #========================================================HER
#error_500_pre= names_wht_data.field('et_f500')    #========================================================HER
flux_850_pre= names_wht_data.field('flux_mJy_perbeam_deboost')
error_850_pre= names_wht_data.field('flux_error_mJy_deboost_calibration')
sn_ratio_confusion_pre= names_wht_data.field('SN_ratio_confusion')






sn_ratio_pre = np.array(sn_ratio_pre)
sn_ratio_pre=np.array([float(i) for i in sn_ratio_pre])
flux_250_pre = np.array(flux_250_pre)
flux_250_pre=np.array([float(i) for i in flux_250_pre])
error_250_pre = np.array(error_250_pre)
error_250_pre=np.array([float(i) for i in error_250_pre])
flux_350_pre = np.array(flux_350_pre)
flux_350_pre=np.array([float(i) for i in flux_350_pre])
error_350_pre = np.array(error_350_pre)
error_350_pre=np.array([float(i) for i in error_350_pre])
flux_500_pre = np.array(flux_500_pre)
flux_500_pre=np.array([float(i) for i in flux_500_pre])
error_500_pre = np.array(error_500_pre)
error_500_pre=np.array([float(i) for i in error_500_pre])
flux_850_pre = np.array(flux_850_pre)
flux_850_pre=np.array([float(i) for i in flux_850_pre])
error_850_pre = np.array(error_850_pre)
error_850_pre=np.array([float(i) for i in error_850_pre])
sn_ratio_confusion_pre = np.array(sn_ratio_confusion_pre)
sn_ratio_confusion_pre=np.array([float(i) for i in sn_ratio_confusion_pre])




sn_ratio = []
flux_250 = []
error_250 = []
flux_350 = []
error_350 = []
flux_500 = []
error_500 = []
flux_850 = []
error_850 = []

for i in range(0, len(sn_ratio_pre)):
    if sn_ratio_confusion_pre[i] >= sn_want: 

        sn_ratio.append(sn_ratio_pre[i])
        flux_250.append(flux_250_pre[i])
        error_250.append(error_250_pre[i])
        flux_350.append(flux_350_pre[i])
        error_350.append(error_350_pre[i])
        flux_500.append(flux_500_pre[i])
        error_500.append(error_500_pre[i])
        flux_850.append(flux_850_pre[i])
        error_850.append(error_850_pre[i])



sn_ratio = np.array(sn_ratio)
sn_ratio=np.array([float(i) for i in sn_ratio])
flux_250 = np.array(flux_250)
flux_250=np.array([float(i) for i in flux_250])
error_250 = np.array(error_250)
error_250=np.array([float(i) for i in error_250])
flux_350 = np.array(flux_350)
flux_350=np.array([float(i) for i in flux_350])
error_350 = np.array(error_350)
error_350=np.array([float(i) for i in error_350])
flux_500 = np.array(flux_500)
flux_500=np.array([float(i) for i in flux_500])
error_500 = np.array(error_500)
error_500=np.array([float(i) for i in error_500])
flux_850 = np.array(flux_850)
flux_850=np.array([float(i) for i in flux_850])
error_850 = np.array(error_850)
error_850=np.array([float(i) for i in error_850])





#--------------------------------for GAMA and NGP fields--------------------------------#========================================================HER
if cluster_list[l]=='G12' or cluster_list[l]=='NGP1' or cluster_list[l]=='NGP2' or cluster_list[l]=='NGP3' or cluster_list[l]=='NGP4' or cluster_list[l]=='NGP5' or cluster_list[l]=='NGP6' or cluster_list[l]=='NGP7' or cluster_list[l]=='NGP8' or cluster_list[l]=='NGP9':
    flux_850=flux_850/1000.0
    error_850=error_850/1000.0
#----------------------------------------------------------------------------------------



f350_f850=[]
error_f350_f850=[]
f250_f350=[]
error_f250_f350=[]
for i in range(len(flux_250)):

    f350_f850.append(flux_250[i]/flux_850[i])
    error_f350_f850.append(abs(flux_350[i]/flux_850[i])*math.sqrt(pow(error_350[i]/flux_350[i],2)+pow(error_850[i]/flux_850[i],2)))
    f250_f350.append(flux_250[i]/flux_350[i])
    error_f250_f350.append(abs(flux_250[i]/flux_350[i])*math.sqrt(pow(error_250[i]/flux_250[i],2)+pow(error_350[i]/flux_350[i],2)))



print f350_f850
print error_f350_f850
print f250_f350
print error_f250_f350


















#flux_250=np.interp(250.0, wavelength_micron, flux_jy)
#flux_350=np.interp(350.0, wavelength_micron, flux_jy)
#flux_500=np.interp(500.0, wavelength_micron, flux_jy)
#flux_850=np.interp(850.0, wavelength_micron, flux_jy)

#print flux_250
#print flux_350
#print flux_500
#print flux_850


#======================redshift======================================================

flux_250_z=[]
flux_350_z=[]
flux_500_z=[]
flux_850_z=[]
color_350_850_z=[]
color_250_350_z=[]

flux_250_z_aless=[]
flux_350_z_aless=[]
flux_500_z_aless=[]
flux_850_z_aless=[]
color_350_850_z_aless=[]
color_250_350_z_aless=[]

flux_250_z_hfls3=[]
flux_350_z_hfls3=[]
flux_500_z_hfls3=[]
flux_850_z_hfls3=[]
color_350_850_z_hfls3=[]
color_250_350_z_hfls3=[]

flux_250_z_eyelash=[]
flux_350_z_eyelash=[]
flux_500_z_eyelash=[]
flux_850_z_eyelash=[]
color_350_850_z_eyelash=[]
color_250_350_z_eyelash=[]

#wavelength_micron_z=[]

#z_list=np.arange(0.0, 5.0, 0.5)

#print z_list

for i in range(0, len(z_list)):

    wavelength_micron_z=wavelength_micron*(1.0+z_list[i])
    flux_250_z.append(np.interp(250.0, wavelength_micron_z, flux_jy))
    flux_350_z.append(np.interp(350.0, wavelength_micron_z, flux_jy))
    flux_500_z.append(np.interp(500.0, wavelength_micron_z, flux_jy))
    flux_850_z.append(np.interp(850.0, wavelength_micron_z, flux_jy))
    color_350_850_z.append(np.interp(350.0, wavelength_micron_z, flux_jy)/np.interp(850.0, wavelength_micron_z, flux_jy))
    color_250_350_z.append(np.interp(250.0, wavelength_micron_z, flux_jy)/np.interp(350.0, wavelength_micron_z, flux_jy))

    wavelength_micron_z_aless=wavelength_micron_aless*(1.0+z_list[i])
    flux_250_z_aless.append(np.interp(250.0, wavelength_micron_z_aless, flux_jy_aless))
    flux_350_z_aless.append(np.interp(350.0, wavelength_micron_z_aless, flux_jy_aless))
    flux_500_z_aless.append(np.interp(500.0, wavelength_micron_z_aless, flux_jy_aless))
    flux_850_z_aless.append(np.interp(850.0, wavelength_micron_z_aless, flux_jy_aless))
    color_350_850_z_aless.append(np.interp(350.0, wavelength_micron_z_aless, flux_jy_aless)/np.interp(850.0, wavelength_micron_z_aless, flux_jy_aless))
    color_250_350_z_aless.append(np.interp(250.0, wavelength_micron_z_aless, flux_jy_aless)/np.interp(350.0, wavelength_micron_z_aless, flux_jy_aless))

    wavelength_micron_z_hfls3=wavelength_micron_hfls3*(1.0+z_list[i])
    flux_250_z_hfls3.append(np.interp(250.0, wavelength_micron_z_hfls3, flux_jy_hfls3))
    flux_350_z_hfls3.append(np.interp(350.0, wavelength_micron_z_hfls3, flux_jy_hfls3))
    flux_500_z_hfls3.append(np.interp(500.0, wavelength_micron_z_hfls3, flux_jy_hfls3))
    flux_850_z_hfls3.append(np.interp(850.0, wavelength_micron_z_hfls3, flux_jy_hfls3))
    color_350_850_z_hfls3.append(np.interp(350.0, wavelength_micron_z_hfls3, flux_jy_hfls3)/np.interp(850.0, wavelength_micron_z_hfls3, flux_jy_hfls3))
    color_250_350_z_hfls3.append(np.interp(250.0, wavelength_micron_z_hfls3, flux_jy_hfls3)/np.interp(350.0, wavelength_micron_z_hfls3, flux_jy_hfls3))

    wavelength_micron_z_eyelash=wavelength_micron_eyelash*(1.0+z_list[i])
    flux_250_z_eyelash.append(np.interp(250.0, wavelength_micron_z_eyelash, flux_jy_eyelash))
    flux_350_z_eyelash.append(np.interp(350.0, wavelength_micron_z_eyelash, flux_jy_eyelash))
    flux_500_z_eyelash.append(np.interp(500.0, wavelength_micron_z_eyelash, flux_jy_eyelash))
    flux_850_z_eyelash.append(np.interp(850.0, wavelength_micron_z_eyelash, flux_jy_eyelash))
    color_350_850_z_eyelash.append(np.interp(350.0, wavelength_micron_z_eyelash, flux_jy_eyelash)/np.interp(850.0, wavelength_micron_z_eyelash, flux_jy_eyelash))
    color_250_350_z_eyelash.append(np.interp(250.0, wavelength_micron_z_eyelash, flux_jy_eyelash)/np.interp(350.0, wavelength_micron_z_eyelash, flux_jy_eyelash))


print color_350_850_z
print color_250_350_z
print color_350_850_z_aless
print color_250_350_z_aless
print color_350_850_z_hfls3
print color_250_350_z_hfls3
print color_350_850_z_eyelash
print color_250_350_z_eyelash



















#ax = fig.add_subplot(3,5,l+1)


if 0<=l<=4:
    ax=ax_all[0,l]
elif 5<=l<=9:
    ax=ax_all[1,l-5]
else:
    ax=ax_all[2,l-10]


line_labels = ["Arp220", "ALESS", "HFLS3", "Cosmic Eyelash"]

l1=ax.plot(color_250_350_z, color_350_850_z, '-o', color='blue', linewidth=1.0, markersize=4.0, label='Arp220') #
l2=ax.plot(color_250_350_z_aless, color_350_850_z_aless, '-o', color='green', linewidth=1.0, markersize=4.0, label='ALESS') #
l3=ax.plot(color_250_350_z_hfls3, color_350_850_z_hfls3, '-o', color='purple', linewidth=1.0, markersize=4.0, label='HFLS3') #
l4=ax.plot(color_250_350_z_eyelash, color_350_850_z_eyelash, '-o', color='grey', linewidth=1.0, markersize=4.0, label='Cosmic Eyelash') #
#plt.legend(loc=2)
#ax.set_xscale('log')
#ax.set_yscale('log')
#plt.xlabel("F250/F350")
#plt.ylabel("F350/F850")
#plt.title("Arp220")    
#plt.rc('font', size=40) 

#plt.show()

ax.plot(color_250_350_z, color_350_850_z, color='blue', linewidth=0.0) # label=cluster_list[l]


#ax.text(color_250_350_z[0], color_350_850_z[0], 'z=0', fontsize=12)
ax.text(color_250_350_z[-1]+0.1, color_350_850_z[-1], 'z=4.5', fontsize=14)

#ax.text(color_250_350_z_aless[0], color_350_850_z_aless[0], 'z=0', fontsize=14)
ax.text(color_250_350_z_aless[4]+0.1, color_350_850_z_aless[4], 'z=1', fontsize=14)
ax.text(color_250_350_z_aless[8]+0.1, color_350_850_z_aless[8], 'z=2', fontsize=14)
#ax.text(color_250_350_z_aless[-1], color_350_850_z_aless[-1], 'z=4.5', fontsize=12, weight='bold')

#ax.text(color_250_350_z_hfls3[0], color_350_850_z_hfls3[0], 'z=0', fontsize=12)
#ax.text(color_250_350_z_hfls3[-1], color_350_850_z_hfls3[-1], 'z=4.5', fontsize=12, weight='bold')

#ax.text(color_250_350_z_eyelash[0], color_350_850_z_eyelash[0], 'z=0', fontsize=12)
#ax.text(color_250_350_z_eyelash[-1], color_350_850_z_eyelash[-1], 'z=4.5', fontsize=12, weight='bold')


ax.text(0.3, 13.0, cluster_list[l], fontsize=20)

#leg = ax.legend(handlelength=0, handletextpad=0, fancybox=True, loc='upper left')
#print leg
#print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
#for item in leg.legendHandles:
#    item.set_visible(False)
#leg.get_frame().set_linewidth(0.0)
#leg.get_frame().set_alpha(0.75)



#ax.set_xlim((0.2, 4.3))
#ax.set_ylim((0.0, 32.0))

ax.set_xlim((0.2, 2.5))
ax.set_ylim((0.0, 15.0))



ax.errorbar(f250_f350, f350_f850, xerr=error_f250_f350, yerr=error_f350_f850, color='r',fmt='o', capsize=2.5, elinewidth=1,markeredgewidth=1)
ax.scatter(f250_f350, f350_f850, color='r', s=50.0) 


ax.axvspan(0.2, 1.3, ymin=0.0, ymax=6.0/15.0, alpha=0.2, color='red', lw=0)
ax.axvspan(0.2, 1.3, ymin=6.0/15.0, ymax=10.0/15.0, alpha=0.2, color='blue', lw=0)
ax.axvspan(1.3, 1.7, ymin=0.0, ymax=10.0/15.0, alpha=0.2, color='blue', lw=0)



#plt.ylim([-10.0,20.0])
#ax.grid()
if 10<=l<=14:
    ax.set_xlabel(r'$S_{250}/S_{350}$')  #==========================================================HERE
#ax.set_title(cluster_list[l])  #==========================================================HERE
if l==0 or l==5 or l==10:
    ax.set_ylabel(r'$S_{350}/S_{850}$')  #==========================================================HERE
plt.rc('font', size=20)


#fig.legend([l1, l2, l3, l4], labels=line_labels, loc='center right')
#plt.figlegend( [l1, l2, l3, l4], labels=line_labels, loc = 'lower center', ncol=5, labelspacing=0. )

handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc=(0.62, 0.15))



#plt.show()






