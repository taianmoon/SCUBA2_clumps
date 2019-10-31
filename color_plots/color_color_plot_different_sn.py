import pyfits
import matplotlib.pyplot as plt
import numpy as np
import math
import glob
import os.path


names_wht = pyfits.open('../NGP9/NGP9_herschel_match_3p0sigma_cat_edgesourcedelete.fits')                            #========================================================HERE
names_wht_data = names_wht[1].data

#RA= names_wht_data.field('RA_1')
#Dec= names_wht_data.field('Dec_1')

flux_250= names_wht_data.field('F250')    #========================================================HERE(Boote.s: F250, et_F250; EG.S/NG.P9/Lockma.n: f250, et250; G1.2/NG.P: F250, E250)
error_250= names_wht_data.field('e250')    #========================================================HER
flux_350= names_wht_data.field('F350')    #========================================================HER
error_350= names_wht_data.field('e350')    #========================================================HER
flux_500= names_wht_data.field('F500')    #========================================================HER
error_500= names_wht_data.field('e500')    #========================================================HER
flux_850= names_wht_data.field('flux_mJy_perbeam_1')
error_850= names_wht_data.field('flux_error_mJy')
sn_ratio= names_wht_data.field('SN_ratio_1')




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
sn_ratio = np.array(sn_ratio)
sn_ratio=np.array([float(i) for i in sn_ratio])

#--------------------------------for GAMA and NGP fields--------------------------------#========================================================HER
flux_850=flux_850/1000.0
error_850=error_850/1000.0
#----------------------------------------------------------------------------------------



f350_f850=[]
error_f350_f850=[]
f250_f350=[]
error_f250_f350=[]

f350_f850_3p5sigma=[]
error_f350_f850_3p5sigma=[]
f250_f350_3p5sigma=[]
error_f250_f350_3p5sigma=[]

f350_f850_4p0sigma=[]
error_f350_f850_4p0sigma=[]
f250_f350_4p0sigma=[]
error_f250_f350_4p0sigma=[]

f350_f850_4p5sigma=[]
error_f350_f850_4p5sigma=[]
f250_f350_4p5sigma=[]
error_f250_f350_4p5sigma=[]

f350_f850_5p0sigma=[]
error_f350_f850_5p0sigma=[]
f250_f350_5p0sigma=[]
error_f250_f350_5p0sigma=[]

for i in range(len(flux_250)):

    if sn_ratio[i] >= 3.0:
        f350_f850.append(flux_250[i]/flux_850[i])
        error_f350_f850.append(abs(flux_350[i]/flux_850[i])*math.sqrt(pow(error_350[i]/flux_350[i],2)+pow(error_850[i]/flux_850[i],2)))
        f250_f350.append(flux_250[i]/flux_350[i])
        error_f250_f350.append(abs(flux_250[i]/flux_350[i])*math.sqrt(pow(error_250[i]/flux_250[i],2)+pow(error_350[i]/flux_350[i],2)))

    if sn_ratio[i] >= 3.5:
        f350_f850_3p5sigma.append(flux_250[i]/flux_850[i])
        error_f350_f850_3p5sigma.append(abs(flux_350[i]/flux_850[i])*math.sqrt(pow(error_350[i]/flux_350[i],2)+pow(error_850[i]/flux_850[i],2)))
        f250_f350_3p5sigma.append(flux_250[i]/flux_350[i])
        error_f250_f350_3p5sigma.append(abs(flux_250[i]/flux_350[i])*math.sqrt(pow(error_250[i]/flux_250[i],2)+pow(error_350[i]/flux_350[i],2)))

    if sn_ratio[i] >= 4.0:
        f350_f850_4p0sigma.append(flux_250[i]/flux_850[i])
        error_f350_f850_4p0sigma.append(abs(flux_350[i]/flux_850[i])*math.sqrt(pow(error_350[i]/flux_350[i],2)+pow(error_850[i]/flux_850[i],2)))
        f250_f350_4p0sigma.append(flux_250[i]/flux_350[i])
        error_f250_f350_4p0sigma.append(abs(flux_250[i]/flux_350[i])*math.sqrt(pow(error_250[i]/flux_250[i],2)+pow(error_350[i]/flux_350[i],2)))

    if sn_ratio[i] >= 4.5:
        f350_f850_4p5sigma.append(flux_250[i]/flux_850[i])
        error_f350_f850_4p5sigma.append(abs(flux_350[i]/flux_850[i])*math.sqrt(pow(error_350[i]/flux_350[i],2)+pow(error_850[i]/flux_850[i],2)))
        f250_f350_4p5sigma.append(flux_250[i]/flux_350[i])
        error_f250_f350_4p5sigma.append(abs(flux_250[i]/flux_350[i])*math.sqrt(pow(error_250[i]/flux_250[i],2)+pow(error_350[i]/flux_350[i],2)))

    if sn_ratio[i] >= 5.0:
        f350_f850_5p0sigma.append(flux_250[i]/flux_850[i])
        error_f350_f850_5p0sigma.append(abs(flux_350[i]/flux_850[i])*math.sqrt(pow(error_350[i]/flux_350[i],2)+pow(error_850[i]/flux_850[i],2)))
        f250_f350_5p0sigma.append(flux_250[i]/flux_350[i])
        error_f250_f350_5p0sigma.append(abs(flux_250[i]/flux_350[i])*math.sqrt(pow(error_250[i]/flux_250[i],2)+pow(error_350[i]/flux_350[i],2)))

print f350_f850
print error_f350_f850
print f250_f350
print error_f250_f350










#==================Add templates in=================================================================================================

#++++++++++++++++++++++++++++++Apr220++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template_file = open("./arp220_template_donley2007_apj_660_167.txt", "r")    #========================================================HERE

lines = template_file.readlines()
template_file.close()

response_curve=[]

for i in range(0, len(lines)):
    separated_lines=lines[i].split() 
    response_curve.append(separated_lines)


response_curve = np.array(response_curve)
wavelength_micron=response_curve[:,0]
flux_jy=response_curve[:,1]


wavelength_micron = np.array(wavelength_micron)
flux_jy = np.array(flux_jy)
wavelength_micron=np.array([float(i) for i in wavelength_micron])  #(wavelength in angstrom)
flux_jy=np.array([float(i) for i in flux_jy])  #(wavelength in angstrom)



#++++++++++++++++++++++++++++++++++++++ALESS+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template_file_aless = open("./aless_average_seds.dat", "r")    #========================================================HERE

lines_aless = template_file_aless.readlines()
template_file_aless.close()

response_curve_aless=[]

for i in range(0, len(lines_aless)):
    separated_lines_aless=lines_aless[i].split() 
    response_curve_aless.append(separated_lines_aless)


response_curve_aless = np.array(response_curve_aless)
wavelength_micron_aless=response_curve_aless[:,0]
flux_jy_aless=response_curve_aless[:,1]


wavelength_micron_aless = np.array(wavelength_micron_aless)
flux_jy_aless = np.array(flux_jy_aless)
wavelength_micron_aless=np.array([float(i) for i in wavelength_micron_aless])  #(wavelength in angstrom)
flux_jy_aless=np.array([float(i) for i in flux_jy_aless])  #(wavelength in angstrom)



#++++++++++++++++++++++++++++++++++HFLS3+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


template_file_hfls3 = open("./HFLS3.txt", "r")    #========================================================HERE

lines_hfls3 = template_file_hfls3.readlines()
template_file_hfls3.close()

response_curve_hfls3=[]

for i in range(0, len(lines_hfls3)):
    separated_lines_hfls3=lines_hfls3[i].split() 
    response_curve_hfls3.append(separated_lines_hfls3)


response_curve_hfls3 = np.array(response_curve_hfls3)
wavelength_micron_hfls3=response_curve_hfls3[:,0]
flux_jy_hfls3=response_curve_hfls3[:,1]


wavelength_micron_hfls3 = np.array(wavelength_micron_hfls3)
flux_jy_hfls3 = np.array(flux_jy_hfls3)
wavelength_micron_hfls3=np.array([float(i) for i in wavelength_micron_hfls3])  #(wavelength in angstrom)
flux_jy_hfls3=np.array([float(i) for i in flux_jy_hfls3])  #(wavelength in angstrom)


wavelength_micron_hfls3=wavelength_micron_hfls3/7.34  #redshift back to restframe (z=6.34)
wavelength_micron_hfls3=wavelength_micron_hfls3*1.0e6  #convert from m to micron
flux_jy_hfls3=flux_jy_hfls3*1000.0  #convert from Jy to mJy




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

#wavelength_micron_z=[]

z_list=np.arange(0.0, 5.0, 0.5)

print z_list

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


print color_350_850_z
print color_250_350_z
print color_350_850_z_aless
print color_250_350_z_aless
print color_350_850_z_hfls3
print color_250_350_z_hfls3







#================================Plotting=================================================================================================================================

fig = plt.figure(figsize=(8.0, 5.0))

ax1 = fig.add_subplot(321)
ax2 = fig.add_subplot(322)
ax3 = fig.add_subplot(323)
ax4 = fig.add_subplot(324)
ax5 = fig.add_subplot(325)
ax6 = fig.add_subplot(326)





ax1.plot(color_250_350_z, color_350_850_z, '-o', color='c', linewidth=2.0, markersize=8.0, label='Arp220')
ax1.plot(color_250_350_z_aless, color_350_850_z_aless, '-o', color='y', linewidth=2.0, markersize=8.0, label='ALESS')
ax1.plot(color_250_350_z_hfls3, color_350_850_z_hfls3, '-o', color='purple', linewidth=2.0, markersize=8.0, label='HFLS3')

#plt.xscale('log')
#plt.yscale('log')
#plt.xlabel("F250/F350")
#plt.ylabel("F350/F850")
#plt.title("Arp220")    
#plt.rc('font', size=40) 

#plt.show()




ax1.text(color_250_350_z[0], color_350_850_z[0], 'z=0', fontsize=20, weight='bold')
ax1.text(color_250_350_z[-1], color_350_850_z[-1], 'z=4.5', fontsize=20, weight='bold')

ax1.text(color_250_350_z_aless[0], color_350_850_z_aless[0], 'z=0', fontsize=20, weight='bold')
#plt.text(color_250_350_z_aless[-1], color_350_850_z_aless[-1], 'z=4.5', fontsize=20, weight='bold')

ax1.text(color_250_350_z_hfls3[0], color_350_850_z_hfls3[0], 'z=0', fontsize=20, weight='bold')
#plt.text(color_250_350_z_hfls3[-1], color_350_850_z_hfls3[-1], 'z=4.5', fontsize=20, weight='bold')



ax2.plot(color_250_350_z, color_350_850_z, '-o', color='c', linewidth=2.0, markersize=8.0)
ax2.plot(color_250_350_z_aless, color_350_850_z_aless, '-o', color='y', linewidth=2.0, markersize=8.0)
ax2.plot(color_250_350_z_hfls3, color_350_850_z_hfls3, '-o', color='purple', linewidth=2.0, markersize=8.0)

#ax2.text(color_250_350_z[0], color_350_850_z[0], 'z=0', fontsize=20, weight='bold')
ax2.text(color_250_350_z[-1], color_350_850_z[-1], 'z=4.5', fontsize=20, weight='bold')
#ax2.text(color_250_350_z_aless[0], color_350_850_z_aless[0], 'z=0', fontsize=20, weight='bold')
#ax2.text(color_250_350_z_aless[-1], color_350_850_z_aless[-1], 'z=4.5', fontsize=20, weight='bold')
#ax2.text(color_250_350_z_hfls3[0], color_350_850_z_hfls3[0], 'z=0', fontsize=20, weight='bold')
#ax2.text(color_250_350_z_hfls3[-1], color_350_850_z_hfls3[-1], 'z=4.5', fontsize=20, weight='bold')



ax3.plot(color_250_350_z, color_350_850_z, '-o', color='c', linewidth=2.0, markersize=8.0,)
ax3.plot(color_250_350_z_aless, color_350_850_z_aless, '-o', color='y', linewidth=2.0, markersize=8.0)
ax3.plot(color_250_350_z_hfls3, color_350_850_z_hfls3, '-o', color='purple', linewidth=2.0, markersize=8.0)

#ax3.text(color_250_350_z[0], color_350_850_z[0], 'z=0', fontsize=20, weight='bold')
ax3.text(color_250_350_z[-1], color_350_850_z[-1], 'z=4.5', fontsize=20, weight='bold')
#ax3.text(color_250_350_z_aless[0], color_350_850_z_aless[0], 'z=0', fontsize=20, weight='bold')
#ax3.text(color_250_350_z_aless[-1], color_350_850_z_aless[-1], 'z=4.5', fontsize=20, weight='bold')
#ax3.text(color_250_350_z_hfls3[0], color_350_850_z_hfls3[0], 'z=0', fontsize=20, weight='bold')
#ax3.text(color_250_350_z_hfls3[-1], color_350_850_z_hfls3[-1], 'z=4.5', fontsize=20, weight='bold')



ax4.plot(color_250_350_z, color_350_850_z, '-o', color='c', linewidth=2.0, markersize=8.0)
ax4.plot(color_250_350_z_aless, color_350_850_z_aless, '-o', color='y', linewidth=2.0, markersize=8.0)
ax4.plot(color_250_350_z_hfls3, color_350_850_z_hfls3, '-o', color='purple', linewidth=2.0, markersize=8.0)

#ax4.text(color_250_350_z[0], color_350_850_z[0], 'z=0', fontsize=20, weight='bold')
ax4.text(color_250_350_z[-1], color_350_850_z[-1], 'z=4.5', fontsize=20, weight='bold')
#ax4.text(color_250_350_z_aless[0], color_350_850_z_aless[0], 'z=0', fontsize=20, weight='bold')
#ax4.text(color_250_350_z_aless[-1], color_350_850_z_aless[-1], 'z=4.5', fontsize=20, weight='bold')
#ax4.text(color_250_350_z_hfls3[0], color_350_850_z_hfls3[0], 'z=0', fontsize=20, weight='bold')
#ax4.text(color_250_350_z_hfls3[-1], color_350_850_z_hfls3[-1], 'z=4.5', fontsize=20, weight='bold')



ax5.plot(color_250_350_z, color_350_850_z, '-o', color='c', linewidth=2.0, markersize=8.0)
ax5.plot(color_250_350_z_aless, color_350_850_z_aless, '-o', color='y', linewidth=2.0, markersize=8.0)
ax5.plot(color_250_350_z_hfls3, color_350_850_z_hfls3, '-o', color='purple', linewidth=2.0, markersize=8.0)

#ax5.text(color_250_350_z[0], color_350_850_z[0], 'z=0', fontsize=20, weight='bold')
ax5.text(color_250_350_z[-1], color_350_850_z[-1], 'z=4.5', fontsize=20, weight='bold')
#ax5.text(color_250_350_z_aless[0], color_350_850_z_aless[0], 'z=0', fontsize=20, weight='bold')
#ax5.text(color_250_350_z_aless[-1], color_350_850_z_aless[-1], 'z=4.5', fontsize=20, weight='bold')
#ax5.text(color_250_350_z_hfls3[0], color_350_850_z_hfls3[0], 'z=0', fontsize=20, weight='bold')
#ax5.text(color_250_350_z_hfls3[-1], color_350_850_z_hfls3[-1], 'z=4.5', fontsize=20, weight='bold')




ax6.plot(color_250_350_z, color_350_850_z, '-o', color='c', linewidth=2.0, markersize=8.0)
ax6.plot(color_250_350_z_aless, color_350_850_z_aless, '-o', color='y', linewidth=2.0, markersize=8.0)
ax6.plot(color_250_350_z_hfls3, color_350_850_z_hfls3, '-o', color='purple', linewidth=2.0, markersize=8.0)

#ax6.text(color_250_350_z[0], color_350_850_z[0], 'z=0', fontsize=20, weight='bold')
ax6.text(color_250_350_z[-1], color_350_850_z[-1], 'z=4.5', fontsize=20, weight='bold')
#ax6.text(color_250_350_z_aless[0], color_350_850_z_aless[0], 'z=0', fontsize=20, weight='bold')
#ax6.text(color_250_350_z_aless[-1], color_350_850_z_aless[-1], 'z=4.5', fontsize=20, weight='bold')
#ax6.text(color_250_350_z_hfls3[0], color_350_850_z_hfls3[0], 'z=0', fontsize=20, weight='bold')
#ax6.text(color_250_350_z_hfls3[-1], color_350_850_z_hfls3[-1], 'z=4.5', fontsize=20, weight='bold')








ax1.errorbar(f250_f350, f350_f850, xerr=error_f250_f350, yerr=error_f350_f850, color='r',fmt='o', capsize=5, elinewidth=2,markeredgewidth=2)
ax1.scatter(f250_f350, f350_f850, color='r', s=100.0) 
ax1.errorbar(f250_f350_3p5sigma, f350_f850_3p5sigma, xerr=error_f250_f350_3p5sigma, yerr=error_f350_f850_3p5sigma,fmt='o', color='r', capsize=5, elinewidth=2,markeredgewidth=2)
ax1.scatter(f250_f350_3p5sigma, f350_f850_3p5sigma, color='r', s=100.0) 
ax1.errorbar(f250_f350_4p0sigma, f350_f850_4p0sigma, xerr=error_f250_f350_4p0sigma, yerr=error_f350_f850_4p0sigma,fmt='o', color='r', capsize=5, elinewidth=2,markeredgewidth=2)
ax1.scatter(f250_f350_4p0sigma, f350_f850_4p0sigma, color='r', s=100.0) 
ax1.errorbar(f250_f350_4p5sigma, f350_f850_4p5sigma, xerr=error_f250_f350_4p5sigma, yerr=error_f350_f850_4p5sigma,fmt='o', color='r', capsize=5, elinewidth=2,markeredgewidth=2)
ax1.scatter(f250_f350_4p5sigma, f350_f850_4p5sigma, color='r', s=100.0) 
ax1.errorbar(f250_f350_5p0sigma, f350_f850_5p0sigma, xerr=error_f250_f350_5p0sigma, yerr=error_f350_f850_5p0sigma,fmt='o', color='r', capsize=5, elinewidth=2,markeredgewidth=2)
ax1.scatter(f250_f350_5p0sigma, f350_f850_5p0sigma, color='r', s=100.0) 



ax2.errorbar(f250_f350, f350_f850, xerr=error_f250_f350, yerr=error_f350_f850, color='r',fmt='o', capsize=5, elinewidth=2,markeredgewidth=2)
ax2.scatter(f250_f350, f350_f850, color='r', s=100.0, label='3.0 sigma') 
ax3.errorbar(f250_f350_3p5sigma, f350_f850_3p5sigma, xerr=error_f250_f350_3p5sigma, yerr=error_f350_f850_3p5sigma,fmt='o', color='r', capsize=5, elinewidth=2,markeredgewidth=2)
ax3.scatter(f250_f350_3p5sigma, f350_f850_3p5sigma, color='r', s=100.0, label='3.5 sigma') 
ax4.errorbar(f250_f350_4p0sigma, f350_f850_4p0sigma, xerr=error_f250_f350_4p0sigma, yerr=error_f350_f850_4p0sigma,fmt='o', color='r', capsize=5, elinewidth=2,markeredgewidth=2)
ax4.scatter(f250_f350_4p0sigma, f350_f850_4p0sigma, color='r', s=100.0, label='4.0 sigma') 
ax5.errorbar(f250_f350_4p5sigma, f350_f850_4p5sigma, xerr=error_f250_f350_4p5sigma, yerr=error_f350_f850_4p5sigma,fmt='o', color='r', capsize=5, elinewidth=2,markeredgewidth=2)
ax5.scatter(f250_f350_4p5sigma, f350_f850_4p5sigma, color='r', s=100.0, label='4.5 sigma') 
ax6.errorbar(f250_f350_5p0sigma, f350_f850_5p0sigma, xerr=error_f250_f350_5p0sigma, yerr=error_f350_f850_5p0sigma,fmt='o', color='r', capsize=5, elinewidth=2,markeredgewidth=2)
ax6.scatter(f250_f350_5p0sigma, f350_f850_5p0sigma, color='r', s=100.0, label='5.0 sigma') 

#ax1.legend(bbox_to_anchor=(2.3, 1), loc=2, borderaxespad=0.)

ax1.legend(loc=2)
ax2.legend(loc=2)
ax3.legend(loc=2)
ax4.legend(loc=2)
ax5.legend(loc=2)
ax6.legend(loc=2)

#plt.ylim([-10.0,20.0])
ax1.grid()
ax2.grid()
ax3.grid()
ax4.grid()
ax5.grid()
ax6.grid()

ax1.set_xlabel('F250/F350')  #==========================================================HERE
ax2.set_xlabel('F250/F350')  #==========================================================HERE
ax3.set_xlabel('F250/F350')  #==========================================================HERE
ax4.set_xlabel('F250/F350')  #==========================================================HERE
ax5.set_xlabel('F250/F350')  #==========================================================HERE
ax6.set_xlabel('F250/F350')  #==========================================================HERE

ax1.set_title('NGP9')  #==========================================================HERE

ax1.set_ylabel('F350/F850')  #==========================================================HERE
ax2.set_ylabel('F350/F850')  #==========================================================HERE
ax3.set_ylabel('F350/F850')  #==========================================================HERE
ax4.set_ylabel('F350/F850')  #==========================================================HERE
ax5.set_ylabel('F350/F850')  #==========================================================HERE
ax6.set_ylabel('F350/F850')  #==========================================================HERE

plt.rc('font', size=20)
plt.show()






