import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm
from astropy.io import fits
import astropy.wcs as wcs
from scipy.interpolate import interp1d


#=====================================Load in source catalogue====================================================================================




template_file = open("./catalogues/NGP9_3p5sigma_edgesourcedelete.cat", "r")    #==================================================================================================HERE

lines = template_file.readlines()
template_file.close()

response_curve=[]

for i in range(0, len(lines)):
    separated_lines=lines[i].split() 
    response_curve.append(separated_lines)


response_curve = np.array(response_curve)
flux_mjy=response_curve[:,6]
flux_mjy_err=response_curve[:,7]
sn_ratio=response_curve[:,4]


flux_mjy = np.array(flux_mjy)
flux_mjy_err = np.array(flux_mjy_err)
sn_ratio = np.array(sn_ratio)
flux_mjy=np.array([float(i) for i in flux_mjy])  #(wavelength in angstrom)
flux_mjy_err=np.array([float(i) for i in flux_mjy_err])  #(wavelength in angstrom)
sn_ratio=np.array([float(i) for i in sn_ratio])  #(wavelength in angstrom)




#-------------import other S/N's-------------------------------

template_file_3p0 = open("./catalogues/NGP9_3p0sigma_edgesourcedelete.cat", "r")    #========================================================HERE

lines_3p0 = template_file_3p0.readlines()
template_file_3p0.close()

response_curve_3p0=[]

for i in range(0, len(lines_3p0)):
    separated_lines_3p0=lines_3p0[i].split() 
    response_curve_3p0.append(separated_lines_3p0)


response_curve_3p0 = np.array(response_curve_3p0)
flux_mjy_3p0=response_curve_3p0[:,6]
flux_mjy_err_3p0=response_curve_3p0[:,7]


flux_mjy_3p0 = np.array(flux_mjy_3p0)
flux_mjy_err_3p0 = np.array(flux_mjy_err_3p0)
flux_mjy_3p0=np.array([float(i) for i in flux_mjy_3p0])  #(wavelength in angstrom)
flux_mjy_err_3p0=np.array([float(i) for i in flux_mjy_err_3p0])  #(wavelength in angstrom)

#flux_mjy_sorted_3p0, flux_mjy_err_sorted_3p0= (list(t) for t in zip(*sorted(zip(flux_mjy_3p0, flux_mjy_err_3p0))))

#flux_mjy_sorted_3p0=flux_mjy_sorted_3p0[::-1]
#flux_mjy_err_sorted_3p0=flux_mjy_err_sorted_3p0[::-1]

print flux_mjy_3p0



flux_mjy_4p0=[]
for i in range(0, len(sn_ratio)):
    if sn_ratio[i] >= 4.0:
        flux_mjy_4p0.append(flux_mjy[i])

print flux_mjy_4p0



flux_mjy_4p5=[]
for i in range(0, len(sn_ratio)):
    if sn_ratio[i] >= 4.5:
        flux_mjy_4p5.append(flux_mjy[i])

print flux_mjy_4p5



flux_mjy_5p0=[]
for i in range(0, len(sn_ratio)):
    if sn_ratio[i] >= 5.0:
        flux_mjy_5p0.append(flux_mjy[i])

print flux_mjy_5p0









#------------------binning------------------------------
#geach_flux=[3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5]  #mJy
geach_flux=[3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5]  #mJy
#geach_flux_widebin=[4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0]
geach_flux_widebin=np.arange(4.0, 100.0, 2.0)


bin_number=[]
bin_number_err=[]

bin_number_3p0=[]
bin_number_err_3p0=[]

bin_number_4p0=[]
bin_number_err_4p0=[]

bin_number_4p5=[]
bin_number_err_4p5=[]

bin_number_5p0=[]
bin_number_err_5p0=[]





#flux_bin_low=3.0
#flux_bin_high=4.0
flux_widebin_low=3.0
flux_widebin_high=5.0


#for j in range(0, len(geach_flux)):
for j in range(0, len(geach_flux_widebin)):
    counter=0.0
    #print flux_bin_low
    print flux_widebin_low
    for i in range(0, len(flux_mjy)):
        #if flux_mjy[i] >= flux_bin_low and flux_mjy[i] < flux_bin_high:
        if flux_mjy[i] >= flux_widebin_low and flux_mjy[i] < flux_widebin_high:
            counter=counter+1.0
    bin_number.append(counter)
    bin_number_err.append(math.sqrt(counter))
    #flux_bin_low=flux_bin_low+1.0
    #flux_bin_high=flux_bin_high+1.0
    flux_widebin_low=flux_widebin_low+2.0
    flux_widebin_high=flux_widebin_high+2.0

#remember to change the scatter x and err when plotting in the end!
print bin_number
print bin_number_err



#-----other SN's---

flux_widebin_low=3.0
flux_widebin_high=5.0
for j in range(0, len(geach_flux_widebin)):
    counter=0.0
    #print flux_widebin_low
    for i in range(0, len(flux_mjy_3p0)):
        if flux_mjy_3p0[i] >= flux_widebin_low and flux_mjy_3p0[i] < flux_widebin_high:
            counter=counter+1.0
    bin_number_3p0.append(counter)
    bin_number_err_3p0.append(math.sqrt(counter))
    flux_widebin_low=flux_widebin_low+2.0
    flux_widebin_high=flux_widebin_high+2.0

#remember to change the scatter x and err when plotting in the end!
print bin_number_3p0
print bin_number_err_3p0



flux_widebin_low=3.0
flux_widebin_high=5.0
for j in range(0, len(geach_flux_widebin)):
    counter=0.0
    #print flux_widebin_low
    for i in range(0, len(flux_mjy_4p0)):
        if flux_mjy_4p0[i] >= flux_widebin_low and flux_mjy_4p0[i] < flux_widebin_high:
            counter=counter+1.0
    bin_number_4p0.append(counter)
    bin_number_err_4p0.append(math.sqrt(counter))
    flux_widebin_low=flux_widebin_low+2.0
    flux_widebin_high=flux_widebin_high+2.0

#remember to change the scatter x and err when plotting in the end!
print bin_number_4p0
print bin_number_err_4p0


flux_widebin_low=3.0
flux_widebin_high=5.0
for j in range(0, len(geach_flux_widebin)):
    counter=0.0
    #print flux_widebin_low
    for i in range(0, len(flux_mjy_4p5)):
        if flux_mjy_4p5[i] >= flux_widebin_low and flux_mjy_4p5[i] < flux_widebin_high:
            counter=counter+1.0
    bin_number_4p5.append(counter)
    bin_number_err_4p5.append(math.sqrt(counter))
    flux_widebin_low=flux_widebin_low+2.0
    flux_widebin_high=flux_widebin_high+2.0

#remember to change the scatter x and err when plotting in the end!
print bin_number_4p5
print bin_number_err_4p5


flux_widebin_low=3.0
flux_widebin_high=5.0
for j in range(0, len(geach_flux_widebin)):
    counter=0.0
    #print flux_widebin_low
    for i in range(0, len(flux_mjy_5p0)):
        if flux_mjy_5p0[i] >= flux_widebin_low and flux_mjy_5p0[i] < flux_widebin_high:
            counter=counter+1.0
    bin_number_5p0.append(counter)
    bin_number_err_5p0.append(math.sqrt(counter))
    flux_widebin_low=flux_widebin_low+2.0
    flux_widebin_high=flux_widebin_high+2.0

#remember to change the scatter x and err when plotting in the end!
print bin_number_5p0
print bin_number_err_5p0



#=========================Estimate the map size in deg^2===============================



hdulist = fits.open('../NGP9/NGP9_snr_170919.fits')      #========================================================================================================================HERE
w = wcs.WCS(hdulist[0].header, hdulist)
NAXIS1=hdulist[0].header['NAXIS1']
NAXIS2=hdulist[0].header['NAXIS2']
CDELT1=hdulist[0].header['CDELT1']
CDELT2=hdulist[0].header['CDELT2']
scidata = hdulist[0].data
hdulist.close()



#hdulist = fits.open('../S2COSMOS/S2CLS_rms_171120.fits')      #=======================================================================================================================HERE Orz
#w = wcs.WCS(hdulist[0].header, hdulist)
#NAXIS1=hdulist[0].header['NAXIS1']
#NAXIS2=hdulist[0].header['NAXIS2']
#CDELT1=hdulist[0].header['CDELT1']
#CDELT2=hdulist[0].header['CDELT2']
#scidata = hdulist[0].data  #================variance: unit of pW^2
#hdulist.close()




pix_value_array=[]

for i in range(0, NAXIS2):
    for j in range(0, NAXIS1):
        pix_value_array.append(scidata[0,i,j])   #Orz  ([i,j] for S2CLS; [0,i,j] for ordinary fields)

pix_value_array = np.array(pix_value_array)
pix_value_array=np.array([float(i) for i in pix_value_array])  #(wavelength in angstrom)

pix_value_array = pix_value_array[~np.isnan(pix_value_array)]

print len(pix_value_array), 'pixels'

area_map=len(pix_value_array)*abs(CDELT1)*abs(CDELT2)

#area_correction_factor=1.0  #======================================================================================================================================HERE
#area_map=area_map*area_correction_factor


print 'the total area from flux map (deg^2) is', area_map, 'degree^2'


#========================Estimate differential number per deg^2==========================

bin_number_perarea=[]
bin_number_perarea_err=[]
for i in range(0, len(bin_number)):
    bin_number_perarea.append(bin_number[i]/area_map)
    bin_number_perarea_err.append(bin_number_err[i]/area_map)


bin_number_perarea = np.array(bin_number_perarea)
bin_number_perarea=np.array([float(i) for i in bin_number_perarea])  #(wavelength in angstrom)
bin_number_perarea_err = np.array(bin_number_perarea_err)
bin_number_perarea_err=np.array([float(i) for i in bin_number_perarea_err])  #(wavelength in angstrom)

print bin_number_perarea
print bin_number_perarea_err


#-----other SN's---

bin_number_perarea_3p0=[]
bin_number_perarea_err_3p0=[]
for i in range(0, len(bin_number_3p0)):
    bin_number_perarea_3p0.append(bin_number_3p0[i]/area_map)
    bin_number_perarea_err_3p0.append(bin_number_err_3p0[i]/area_map)


bin_number_perarea_3p0 = np.array(bin_number_perarea_3p0)
bin_number_perarea_3p0=np.array([float(i) for i in bin_number_perarea_3p0])  #(wavelength in angstrom)
bin_number_perarea_err_3p0 = np.array(bin_number_perarea_err_3p0)
bin_number_perarea_err_3p0=np.array([float(i) for i in bin_number_perarea_err_3p0])  #(wavelength in angstrom)

print bin_number_perarea_3p0
print bin_number_perarea_err_3p0



bin_number_perarea_4p0=[]
bin_number_perarea_err_4p0=[]
for i in range(0, len(bin_number_4p0)):
    bin_number_perarea_4p0.append(bin_number_4p0[i]/area_map)
    bin_number_perarea_err_4p0.append(bin_number_err_4p0[i]/area_map)


bin_number_perarea_4p0 = np.array(bin_number_perarea_4p0)
bin_number_perarea_4p0=np.array([float(i) for i in bin_number_perarea_4p0])  #(wavelength in angstrom)
bin_number_perarea_err_4p0 = np.array(bin_number_perarea_err_4p0)
bin_number_perarea_err_4p0=np.array([float(i) for i in bin_number_perarea_err_4p0])  #(wavelength in angstrom)

print bin_number_perarea_4p0
print bin_number_perarea_err_4p0



bin_number_perarea_4p5=[]
bin_number_perarea_err_4p5=[]
for i in range(0, len(bin_number_4p5)):
    bin_number_perarea_4p5.append(bin_number_4p5[i]/area_map)
    bin_number_perarea_err_4p5.append(bin_number_err_4p5[i]/area_map)


bin_number_perarea_4p5 = np.array(bin_number_perarea_4p5)
bin_number_perarea_4p5=np.array([float(i) for i in bin_number_perarea_4p5])  #(wavelength in angstrom)
bin_number_perarea_err_4p5 = np.array(bin_number_perarea_err_4p5)
bin_number_perarea_err_4p5=np.array([float(i) for i in bin_number_perarea_err_4p5])  #(wavelength in angstrom)

print bin_number_perarea_4p5
print bin_number_perarea_err_4p5



bin_number_perarea_5p0=[]
bin_number_perarea_err_5p0=[]
for i in range(0, len(bin_number_5p0)):
    bin_number_perarea_5p0.append(bin_number_5p0[i]/area_map)
    bin_number_perarea_err_5p0.append(bin_number_err_5p0[i]/area_map)


bin_number_perarea_5p0 = np.array(bin_number_perarea_5p0)
bin_number_perarea_5p0=np.array([float(i) for i in bin_number_perarea_5p0])  #(wavelength in angstrom)
bin_number_perarea_err_5p0 = np.array(bin_number_perarea_err_5p0)
bin_number_perarea_err_5p0=np.array([float(i) for i in bin_number_perarea_err_5p0])  #(wavelength in angstrom)

print bin_number_perarea_5p0
print bin_number_perarea_err_5p0


#====================================================Adding points from literature============================================

#geach_flux=[3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5]  #mJy
geach_diff_number=[451.0, 204.4, 102.6, 56.1, 32.5, 18.0, 9.8, 5.8, 3.4, 2.1, 0.8, 0.5, 0.3]  #/deg^2 /mJy
geach_diff_number_err_up=[17.1, 9.3, 6.0, 4.3, 3.2, 2.5, 1.9, 1.5, 1.2, 1.1, 0.8, 0.7, 0.6] #Poisson error
geach_diff_number_err_low=[16.4, 8.9, 5.7, 4.0, 2.9, 2.2, 1.6, 1.2, 0.9, 0.7, 0.4, 0.3, 0.2]

asymmetric_error = [geach_diff_number_err_low, geach_diff_number_err_up]




#=======================plotting========================================================

plt.errorbar(geach_flux_widebin, bin_number_perarea_3p0, yerr=bin_number_perarea_err_3p0, color='blue',fmt='o', capsize=5, elinewidth=2,markeredgewidth=2)
plt.plot(geach_flux_widebin, bin_number_perarea_3p0, label='>3.0 sigma', color='blue', linewidth=2.0)
plt.scatter(geach_flux_widebin, bin_number_perarea_3p0, color='blue')



#plt.errorbar(geach_flux, bin_number_perarea, yerr=bin_number_perarea_err, color='blue',fmt='o', capsize=5, elinewidth=2,markeredgewidth=2)
plt.errorbar(geach_flux_widebin, bin_number_perarea, yerr=bin_number_perarea_err, color='g',fmt='o', capsize=5, elinewidth=2,markeredgewidth=2)
plt.plot(geach_flux_widebin, bin_number_perarea, label='>3.5 sigma', color='g', linewidth=2.0)
#plt.scatter(geach_flux, bin_number_perarea, label='This work')
plt.scatter(geach_flux_widebin, bin_number_perarea, color='g')


plt.errorbar(geach_flux_widebin, bin_number_perarea_4p0, yerr=bin_number_perarea_err_4p0, color='c',fmt='o', capsize=5, elinewidth=2,markeredgewidth=2)
plt.plot(geach_flux_widebin, bin_number_perarea_4p0, label='>4.0 sigma', color='c', linewidth=2.0)
plt.scatter(geach_flux_widebin, bin_number_perarea_4p0, color='c')

plt.errorbar(geach_flux_widebin, bin_number_perarea_4p5, yerr=bin_number_perarea_err_4p5, color='m',fmt='o', capsize=5, elinewidth=2,markeredgewidth=2)
plt.plot(geach_flux_widebin, bin_number_perarea_4p5, label='>4.5 sigma', color='m', linewidth=2.0)
plt.scatter(geach_flux_widebin, bin_number_perarea_4p5, color='m')

plt.errorbar(geach_flux_widebin, bin_number_perarea_5p0, yerr=bin_number_perarea_err_5p0, color='y',fmt='o', capsize=5, elinewidth=2,markeredgewidth=2)
plt.plot(geach_flux_widebin, bin_number_perarea_5p0, label='>5.0 sigma', color='y', linewidth=2.0)
plt.scatter(geach_flux_widebin, bin_number_perarea_5p0, color='y')


plt.plot(geach_flux, geach_diff_number, label='Geach 2017', color='red', linewidth=2.0)
plt.scatter(geach_flux, geach_diff_number, color='red')
plt.errorbar(geach_flux, geach_diff_number, yerr=asymmetric_error, color='red',fmt='o', capsize=5, elinewidth=2,markeredgewidth=2)


plt.tick_params(width=2, length=16, which='major')
plt.tick_params(width=2, length=5, which='minor')


plt.xscale('log',nonposy='clip')
plt.yscale('log',nonposy='clip')


plt.legend(loc=3)
plt.grid()
plt.xlabel('flux (mJy)')  
plt.title('NGP9')  #===================================================================================================================================================================HERE
plt.ylabel('N_differential (deg^(-2) mJy^(-1))')  
plt.rc('font', size=30)
plt.show()







