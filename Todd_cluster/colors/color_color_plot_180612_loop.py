import pyfits
import matplotlib.pyplot as plt
import numpy as np
import math
import glob
import os.path


cluster_list=['Planck18p194', 'Planck18p735', 'Planck24p194', 'PLCK_DU_G045.7-41.2', 'PLCK_DU_G059.1-67.1', 'PLCK_DU_G073.4-57.5', 'PLCK_G006.1+61.8', 'PLCK_G009.8+72.6', 'PLCK_G056.7+62.6', 'PLCK_G068.3+31.9', 'PLCK_G075.1+33.2', 'PLCK_G077.7+32.6', 'PLCK_G078.9+48.2', 'PLCK_G082.5+38.4', 'PLCK_G083.3+51.0', 'PLCK_G091.9+43.0', 'PLCK_G093.6+55.9', 'PLCK_G132.9-76.0', 'PLCK_G144.1+81.0', 'PLCK_G160.7+41.0', 'PLCK_G162.1-59.3', 'PLCK_G165.8+45.3', 'PLCK_G173.8+59.3', 'PLCK_G177.0+35.9', 'PLCK_G179.3+50.7', 'PLCK_G186.3-72.7', 'PLCK_G186.6+66.7', 'PLCK_G188.6-68.9', 'PLCK_G191.3+62.0', 'PLCK_G191.8-83.4', 'PLCK_G201.1+50.7', 'PLCK_G213.0+65.9', 'PLCK_G223.9+41.2', 'PLCK_G328.9+71.4', 'PLCK_G49.6-42.9', 'PLCK_G84.0-71.5', 'PLCK_HZ_G038.0-51.5', 'PLCK_HZ_G067.2-63.8', 'PLCK_HZ_G103.1-73.6', 'PLCK_HZ_G106.8-83.3', 'PLCK_HZ_G119.4-76.6', 'PLCK_HZ_G132.6-81.1', 'PLCK_HZ_G171.1-78.7', 'PLCK_HZ_G173.9+57.0', 'PLCK_HZ_G176.6+59.0', 'PLCK_HZ_G214.1+48.3']

#three categories
#cluster_list=['PLCK_G165.8+45.3', 'PLCK_G082.5+38.4', 'PLCK_G49.6-42.9']


overd_cat=[2, 1, 2, 2, 3, 1, 3, 3, 3, 3, 3, 3, 1, 3, 3, 3, 1, 1, 2, 3, 1, 3, 1, 2, 3, 2, 3, 3, 3, 3, 2, 3, 2, 3, 1, 2, 3, 3, 2, 3, 2, 1, 3, 1, 3, 3]


sn_want=3.5  #========================================================================HERE: change threshold for sn ratio
flux_want=0.0  #========================================================================HERE: change threshold for flux (mJy) for checking lens
flux_threshold=10.0  #========================================================================HERE: flux threshold for plotting different colors


#fig, ax_all =plt.subplots(3,5, figsize=(16.0, 10.0), sharex=True, sharey=True)
fig, ax_all =plt.subplots(1,3, figsize=(16.0, 10.0), sharex=True, sharey=True)
#fig = plt.figure(figsize=(16.0, 10.0))

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



#++++++++++++++++++++++++++++++++++++++Cosmic Eyelash+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template_file_eyelash = open("./cosmic_eyelash_sed_josh.txt", "r")    #========================================================HERE

lines_eyelash = template_file_eyelash.readlines()
template_file_eyelash.close()

response_curve_eyelash=[]

for i in range(0, len(lines_eyelash)):
    separated_lines_eyelash=lines_eyelash[i].split() 
    response_curve_eyelash.append(separated_lines_eyelash)


response_curve_eyelash = np.array(response_curve_eyelash)
wavelength_micron_eyelash=response_curve_eyelash[:,0]
flux_jy_eyelash=response_curve_eyelash[:,1]


wavelength_micron_eyelash = np.array(wavelength_micron_eyelash)
flux_jy_eyelash = np.array(flux_jy_eyelash)
wavelength_micron_eyelash=np.array([float(i) for i in wavelength_micron_eyelash])  #(wavelength in angstrom)
flux_jy_eyelash=np.array([float(i) for i in flux_jy_eyelash])  #(wavelength in angstrom)



#===================================================================================================================================

#z_list=np.arange(0.0, 5.0, 0.5)
z_list=np.arange(0.0, 4.75, 0.25)  #================================================================:HERE








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


#
#
#plt.plot(color_250_350_z, color_350_850_z, '-o', color='blue', linewidth=1.0, markersize=4.0, label='Arp220') #
#plt.plot(color_250_350_z_aless, color_350_850_z_aless, '-o', color='green', linewidth=1.0, markersize=4.0, label='ALESS') #
#plt.plot(color_250_350_z_hfls3, color_350_850_z_hfls3, '-o', color='purple', linewidth=1.0, markersize=4.0, label='HFLS3') #
#plt.plot(color_250_350_z_eyelash, color_350_850_z_eyelash, '-o', color='grey', linewidth=1.0, markersize=4.0, label='Cosmic Eyelash') #
##plt.legend(loc=2)
##ax.set_xscale('log')
##ax.set_yscale('log')
##plt.xlabel("F250/F350")
##plt.ylabel("F350/F850")
##plt.title("Arp220")    
##plt.rc('font', size=40) 
#
##plt.show()
#
##ax.plot(color_250_350_z, color_350_850_z, color='blue', linewidth=0.0) # label=cluster_list[l]
#
##ax.text(color_250_350_z[0], color_350_850_z[0], 'z=0', fontsize=12)
#plt.text(color_250_350_z[-1]+0.1, color_350_850_z[-1], 'z=4.5', fontsize=14)
#
##ax.text(color_250_350_z_aless[0], color_350_850_z_aless[0], 'z=0', fontsize=14)
#plt.text(color_250_350_z_aless[4]+0.1, color_350_850_z_aless[4], 'z=1', fontsize=14)
#plt.text(color_250_350_z_aless[8]+0.1, color_350_850_z_aless[8], 'z=2', fontsize=14)
##ax.text(color_250_350_z_aless[-1], color_350_850_z_aless[-1], 'z=4.5', fontsize=12, weight='bold')
#
##ax.text(color_250_350_z_hfls3[0], color_350_850_z_hfls3[0], 'z=0', fontsize=12)
##ax.text(color_250_350_z_hfls3[-1], color_350_850_z_hfls3[-1], 'z=4.5', fontsize=12, weight='bold')
#
##ax.text(color_250_350_z_eyelash[0], color_350_850_z_eyelash[0], 'z=0', fontsize=12)
##ax.text(color_250_350_z_eyelash[-1], color_350_850_z_eyelash[-1], 'z=4.5', fontsize=12, weight='bold')
#
#plt.xlim((0.2, 2.5))
#plt.ylim((0.0, 15.0))
#
#













for l in range(0, len(cluster_list)):

    print cluster_list[l]
    #execfile("./color_color_plot_180517.py")
    #execfile("./color_color_plot_180905.py")
    #execfile("./color_color_plot_181127.py")
    #execfile("./color_color_plot_190304.py")
    #execfile("./color_color_plot_190612.py")
    #execfile("./color_color_plot_190625.py")
    execfile("./color_color_plot_190815.py")


#plt.scatter([], [], c=col, s=s, edgecolors='none', label='Blue: >'+str(flux_threshold)+' mJy') 
#plt.legend(loc=2)

plt.show()


