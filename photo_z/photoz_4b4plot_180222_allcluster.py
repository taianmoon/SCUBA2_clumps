import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm
from astropy.io import fits
import astropy.wcs as wcs
from scipy.interpolate import interp1d
import pyfits
import os
import subprocess


#=======================import the catalogues, get flux and err=============================================

names_wht = pyfits.open('./Allcluster_herschel_match_3p5sigma_cat_edgesourcedelete_allscuba2source.fits')                            #========================================================HERE
names_wht_data = names_wht[1].data


flux_250_observed_pre= names_wht_data.field('F250')    #========================================================HERE(Boote.s: F250, et_F250; EG.S/Lockma.n: f250, et250; G1.2/NG.P: F250, E250)
error_250_observed_pre= names_wht_data.field('et_f250')    #========================================================HER
flux_350_observed_pre= names_wht_data.field('F350')    #========================================================HER
error_350_observed_pre= names_wht_data.field('et_f350')    #========================================================HER
flux_500_observed_pre= names_wht_data.field('F500')    #========================================================HER
error_500_observed_pre= names_wht_data.field('et_f500')    #========================================================HER
flux_850_observed_pre= names_wht_data.field('flux_mJy_perbeam_1')
error_850_observed_pre= names_wht_data.field('flux_error_mJy_1')
RA_pre= names_wht_data.field('RA_1')
Dec_pre= names_wht_data.field('Dec_1')
source_name_pre= names_wht_data.field('name')


flux_250_observed_pre = np.array(flux_250_observed_pre)
flux_250_observed_pre=np.array([float(i) for i in flux_250_observed_pre])
error_250_observed_pre = np.array(error_250_observed_pre)
error_250_observed_pre=np.array([float(i) for i in error_250_observed_pre])
flux_350_observed_pre = np.array(flux_350_observed_pre)
flux_350_observed_pre=np.array([float(i) for i in flux_350_observed_pre])
error_350_observed_pre = np.array(error_350_observed_pre)
error_350_observed_pre=np.array([float(i) for i in error_350_observed_pre])
flux_500_observed_pre = np.array(flux_500_observed_pre)
flux_500_observed_pre=np.array([float(i) for i in flux_500_observed_pre])
error_500_observed_pre = np.array(error_500_observed_pre)
error_500_observed_pre=np.array([float(i) for i in error_500_observed_pre])
flux_850_observed_pre = np.array(flux_850_observed_pre)
flux_850_observed_pre=np.array([float(i) for i in flux_850_observed_pre])
error_850_observed_pre = np.array(error_850_observed_pre)
error_850_observed_pre=np.array([float(i) for i in error_850_observed_pre])




flux_250_observed=[]
error_250_observed=[]
flux_350_observed=[]
error_350_observed=[]
flux_500_observed=[]
error_500_observed=[]
flux_850_observed=[]
error_850_observed=[]
RA=[]
Dec=[]
source_name=[]

for i in range(0, len(flux_250_observed_pre)):

    if  math.isnan(flux_250_observed_pre[i]) == False:

        flux_250_observed.append(flux_250_observed_pre[i])
        error_250_observed.append(error_250_observed_pre[i])
        flux_350_observed.append(flux_350_observed_pre[i])
        error_350_observed.append(error_350_observed_pre[i])
        flux_500_observed.append(flux_500_observed_pre[i])
        error_500_observed.append(error_500_observed_pre[i])
        flux_850_observed.append(flux_850_observed_pre[i])
        error_850_observed.append(error_850_observed_pre[i])
        RA.append(RA_pre[i])
        Dec.append(Dec_pre[i])
        source_name.append(source_name_pre[i])


print source_name


flux_250_observed = np.array(flux_250_observed)
flux_250_observed=np.array([float(i) for i in flux_250_observed])
error_250_observed = np.array(error_250_observed)
error_250_observed=np.array([float(i) for i in error_250_observed])
flux_350_observed = np.array(flux_350_observed)
flux_350_observed=np.array([float(i) for i in flux_350_observed])
error_350_observed = np.array(error_350_observed)
error_350_observed=np.array([float(i) for i in error_350_observed])
flux_500_observed = np.array(flux_500_observed)
flux_500_observed=np.array([float(i) for i in flux_500_observed])
error_500_observed = np.array(error_500_observed)
error_500_observed=np.array([float(i) for i in error_500_observed])
flux_850_observed = np.array(flux_850_observed)
flux_850_observed=np.array([float(i) for i in flux_850_observed])
error_850_observed = np.array(error_850_observed)
error_850_observed=np.array([float(i) for i in error_850_observed])



#--------------------------------for GAMA and NGP fields--------------------------------#====================================================================================HERE

flux_250_observed=flux_250_observed*1000.0
error_250_observed=error_250_observed*1000.0
flux_350_observed=flux_350_observed*1000.0
error_350_observed=error_350_observed*1000.0
flux_500_observed=flux_500_observed*1000.0
error_500_observed=error_500_observed*1000.0
#----------------------------------------------------------------------------------------




#=====================import template(s) redshift it, interpolate the model values at the correspoinding bands, and create parameter space (z and L --> delx and dely)==================



template_file = open("./arp220_template_donley2007_apj_660_167.txt", "r")

lines = template_file.readlines()
template_file.close()

response_curve=[]

for i in range(0, len(lines)):
    separated_lines=lines[i].split() 
    response_curve.append(separated_lines)


response_curve = np.array(response_curve)
wavelength_micron=response_curve[:,0]
flux_mjy=response_curve[:,1]


wavelength_micron = np.array(wavelength_micron)
flux_mjy = np.array(flux_mjy)
wavelength_micron=np.array([float(i) for i in wavelength_micron])  #(wavelength in angstrom)
flux_mjy=np.array([float(i) for i in flux_mjy])  #(wavelength in angstrom)



flux_250_z=[]
flux_350_z=[]
flux_500_z=[]
flux_850_z=[]

z_list=np.arange(0.0, 6.0, 0.05)  #========================================================HERE: change parameter grids
dely_list=np.arange(0.0, 10.0, 0.05)  #========================================================HERE: change parameter grids




z_list_space=[]
dely_list_space=[]
for i in range(0, len(z_list)):
    for j in range(0, len(dely_list)):

        wavelength_micron_z = wavelength_micron*(1.0+z_list[i])
        flux_mjy_z = flux_mjy*dely_list[j]
        flux_250_z.append(np.interp(250.0, wavelength_micron_z, flux_mjy_z))
        flux_350_z.append(np.interp(350.0, wavelength_micron_z, flux_mjy_z))
        flux_500_z.append(np.interp(500.0, wavelength_micron_z, flux_mjy_z))
        flux_850_z.append(np.interp(850.0, wavelength_micron_z, flux_mjy_z))
        z_list_space.append(z_list[i])
        dely_list_space.append(dely_list[j])

z_list_space = np.array(z_list_space)
z_list_space=np.array([float(i) for i in z_list_space])  #(wavelength in angstrom)
dely_list_space = np.array(dely_list_space)
dely_list_space=np.array([float(i) for i in dely_list_space])  #(wavelength in angstrom)






execfile("./other_template.py")







#=====================Calculate residual and least sqaure===========================================================


#-----------prerequisite for contourf plot-----------------------
from numpy import linspace, meshgrid
from matplotlib.mlab import griddata
import matplotlib.colors as colors

def grid(x, y, z, resX=100, resY=100):  #=======================================================================================================HERE: change contour grid density
#    "Convert 3 column data to matplotlib grid"
    xi = linspace(min(x), max(x), resX)
    yi = linspace(min(y), max(y), resY)
    Z = griddata(x, y, z, xi, yi)
    X, Y = meshgrid(xi, yi)
    return X, Y, Z
#----------------------------------------------------------------



best_redshift_all=[]
best_dely_all=[]
best_redshift_all_aless=[]
best_dely_all_aless=[]
best_redshift_all_hfls3=[]
best_dely_all_hfls3=[]

best_redshift_all_marginalized=[]
best_dely_all_marginalized=[]

min_chi2_all=[]
min_chi2_all_aless=[]
min_chi2_all_hfls3=[]

max_likelihood_marginalized_all=[]

#likelihood_marginalize_over_z_all=[]
multiply_value_z=[1.0]*len(z_list)
multiply_value_a=[1.0]*len(dely_list)

for i in range(0, len(flux_250_observed)):  #looping different sources

    print "now source "+str(i)+" out of "+str(len(flux_250_observed))

    residual_squared_sum=[]
    residual_squared_sum_aless=[]
    residual_squared_sum_hfls3=[]

    likelihood_sum=[]
    likelihood_sum_aless=[]
    likelihood_sum_hfls3=[]
    likelihood_sum_marginalize=[]

    for j in range(0, len(flux_250_z)):  #looping different parameters
        
        residual_squared_250 = ((flux_250_observed[i]-flux_250_z[j])*(flux_250_observed[i]-flux_250_z[j])/(error_250_observed[i]*error_250_observed[i]))       #=======================Check here!!!! 
        residual_squared_350 = ((flux_350_observed[i]-flux_350_z[j])*(flux_350_observed[i]-flux_350_z[j])/(error_350_observed[i]*error_350_observed[i]))        
        residual_squared_500 = ((flux_500_observed[i]-flux_500_z[j])*(flux_500_observed[i]-flux_500_z[j])/(error_500_observed[i]*error_500_observed[i]))       
        residual_squared_850 = ((flux_850_observed[i]-flux_850_z[j])*(flux_850_observed[i]-flux_850_z[j])/(error_850_observed[i]*error_850_observed[i]))

        residual_squared_250_aless = ((flux_250_observed[i]-flux_250_z_aless[j])*(flux_250_observed[i]-flux_250_z_aless[j])/(error_250_observed[i]*error_250_observed[i]))       #=======================Check here!!!! 
        residual_squared_350_aless = ((flux_350_observed[i]-flux_350_z_aless[j])*(flux_350_observed[i]-flux_350_z_aless[j])/(error_350_observed[i]*error_350_observed[i]))        
        residual_squared_500_aless = ((flux_500_observed[i]-flux_500_z_aless[j])*(flux_500_observed[i]-flux_500_z_aless[j])/(error_500_observed[i]*error_500_observed[i]))       
        residual_squared_850_aless = ((flux_850_observed[i]-flux_850_z_aless[j])*(flux_850_observed[i]-flux_850_z_aless[j])/(error_850_observed[i]*error_850_observed[i]))


        residual_squared_250_hfls3 = ((flux_250_observed[i]-flux_250_z_hfls3[j])*(flux_250_observed[i]-flux_250_z_hfls3[j])/(error_250_observed[i]*error_250_observed[i]))       #=======================Check here!!!! 
        residual_squared_350_hfls3 = ((flux_350_observed[i]-flux_350_z_hfls3[j])*(flux_350_observed[i]-flux_350_z_hfls3[j])/(error_350_observed[i]*error_350_observed[i]))        
        residual_squared_500_hfls3 = ((flux_500_observed[i]-flux_500_z_hfls3[j])*(flux_500_observed[i]-flux_500_z_hfls3[j])/(error_500_observed[i]*error_500_observed[i]))       
        residual_squared_850_hfls3 = ((flux_850_observed[i]-flux_850_z_hfls3[j])*(flux_850_observed[i]-flux_850_z_hfls3[j])/(error_850_observed[i]*error_850_observed[i]))




        residual_squared_sum.append((residual_squared_250+residual_squared_350+residual_squared_500+residual_squared_850))
        residual_squared_sum_aless.append((residual_squared_250_aless+residual_squared_350_aless+residual_squared_500_aless+residual_squared_850_aless))
        residual_squared_sum_hfls3.append((residual_squared_250_hfls3+residual_squared_350_hfls3+residual_squared_500_hfls3+residual_squared_850_hfls3))

        likelihood_sum.append(math.exp(-0.5*(residual_squared_250+residual_squared_350+residual_squared_500+residual_squared_850)))
        likelihood_sum_aless.append(math.exp(-0.5*(residual_squared_250_aless+residual_squared_350_aless+residual_squared_500_aless+residual_squared_850_aless)))
        likelihood_sum_hfls3.append(math.exp(-0.5*(residual_squared_250_hfls3+residual_squared_350_hfls3+residual_squared_500_hfls3+residual_squared_850_hfls3)))
        likelihood_sum_marginalize.append(math.exp(-0.5*(residual_squared_250+residual_squared_350+residual_squared_500+residual_squared_850))+math.exp(-0.5*(residual_squared_250_aless+residual_squared_350_aless+residual_squared_500_aless+residual_squared_850_aless))+math.exp(-0.5*(residual_squared_250_hfls3+residual_squared_350_hfls3+residual_squared_500_hfls3+residual_squared_850_hfls3)))

    print 'Min residual squared using arp220 for this source is ', min(residual_squared_sum)
    print 'Min residual squared using ALESS for this source is ', min(residual_squared_sum_aless)
    print 'Min residual squared using HFLS3 for this source is ', min(residual_squared_sum_hfls3)
    print 'Max likelihood (marginalized) for this source is ', max(likelihood_sum_marginalize)

    min_chi2_all.append(min(residual_squared_sum))
    min_chi2_all_aless.append(min(residual_squared_sum_aless))
    min_chi2_all_hfls3.append(min(residual_squared_sum_hfls3))

    max_likelihood_marginalized_all.append(max(likelihood_sum_marginalize))

    for k in range(0, len(flux_250_z)):
        if residual_squared_sum[k] == min(residual_squared_sum):
            print 'its index (arp220) in parameter space is ', k 
            best_redshift=z_list_space[k]
            best_dely=dely_list_space[k]


    for l in range(0, len(flux_250_z_aless)):
        if residual_squared_sum_aless[l] == min(residual_squared_sum_aless):
            print 'its index (ALESS) in parameter space is ', l 
            best_redshift_aless=z_list_space[l]
            best_dely_aless=dely_list_space[l]



    for m in range(0, len(flux_250_z_hfls3)):
        if residual_squared_sum_hfls3[m] == min(residual_squared_sum_hfls3):
            print 'its index (HFLS3) in parameter space is ', m 
            best_redshift_hfls3=z_list_space[m]
            best_dely_hfls3=dely_list_space[m]


    for n in range(0, len(likelihood_sum_marginalize)):
        if likelihood_sum_marginalize[n] == max(likelihood_sum_marginalize):
            print 'its index (marginalized) in parameter space is ', n 
            best_redshift_marginalized=z_list_space[n]
            best_dely_marginalized=dely_list_space[n]


    best_wavelength_micron=wavelength_micron*(1.0+best_redshift)
    best_flux_mjy = flux_mjy*best_dely

    best_wavelength_micron_aless=wavelength_micron_aless*(1.0+best_redshift_aless)
    best_flux_mjy_aless = flux_mjy_aless*best_dely_aless

    best_wavelength_micron_hfls3=wavelength_micron_hfls3*(1.0+best_redshift_hfls3)
    best_flux_mjy_hfls3 = flux_mjy_hfls3*best_dely_hfls3

    best_redshift_all.append(best_redshift)
    best_dely_all.append(best_dely)

    best_redshift_all_aless.append(best_redshift_aless)
    best_dely_all_aless.append(best_dely_aless)

    best_redshift_all_hfls3.append(best_redshift_hfls3)
    best_dely_all_hfls3.append(best_dely_hfls3)

    best_redshift_all_marginalized.append(best_redshift_marginalized)
    best_dely_all_marginalized.append(best_dely_marginalized)

    print 'best redshift (arp220) for this source is ', best_redshift
    print 'best dely (arp220) for this source is ', best_dely
    print 'best redshift (ALESS) for this source is ', best_redshift_aless
    print 'best dely (ALESS) for this source is ', best_dely_aless
    print 'best redshift (HFLS3) for this source is ', best_redshift_hfls3
    print 'best dely (HFLS3) for this source is ', best_dely_hfls3
    print 'best redshift (marginalized) for this source is ', best_redshift_marginalized
    print 'best dely (marginalized) for this source is ', best_dely_marginalized

    #--------------------------marginalizing over a and z---------------------------------------------------


    likelihood_marginalize_over_z=[]
    for k in range(0, len(z_list)):

        counter_sum=0.0
        for l in range(0, len(z_list_space)): 

            if z_list_space[l]==z_list[k]:
                counter_sum=counter_sum+likelihood_sum_marginalize[l]

        likelihood_marginalize_over_z.append(counter_sum)

    #print 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
    #print likelihood_marginalize_over_z

    likelihood_marginalize_over_a=[]
    for m in range(0, len(dely_list)):

        counter_sum=0.0
        for n in range(0, len(dely_list_space)): 

            if dely_list_space[n]==dely_list[m]:
                counter_sum=counter_sum+likelihood_sum_marginalize[n]

        likelihood_marginalize_over_a.append(counter_sum)

    #print 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
    #print likelihood_marginalize_over_a


    #-------------------multiply marginalized likelihood with previous sources-----------


    for o in range(0, len(z_list)):
        #marginalize_start=marginalize_start*likelihood_marginalize_over_z[o]
        #likelihood_marginalize_over_z_all.append(marginalize_start)


        multiply_value_z[o]=multiply_value_z[o]*likelihood_marginalize_over_z[o]
        #multiply_value[1]=multiply_value*likelihood_marginalize_over_z[1]

    #print 'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz'
    #print multiply_value


    for p in range(0, len(dely_list)):

        multiply_value_a[p]=multiply_value_a[p]*likelihood_marginalize_over_a[p]


    #print 'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz'
    #print multiply_value


    #=========================Plotting=====================================================================

    #fig = plt.figure(figsize=(8.0, 5.0))

    #ax1 = fig.add_subplot(221)

    #ax2 = fig.add_subplot(222)
 
    #ax3 = fig.add_subplot(223)
    #ax4 = fig.add_subplot(224)


    ##1. chi-squared histogram plot

    
    #X, Y, Z = grid(z_list_space, dely_list_space, residual_squared_sum)
    #X_aless, Y_aless, Z_aless = grid(z_list_space, dely_list_space, residual_squared_sum_aless)
    #X_hfls3, Y_hfls3, Z_hfls3 = grid(z_list_space, dely_list_space, residual_squared_sum_hfls3)



    #im=ax2.pcolor(X, Y, Z, norm=colors.LogNorm(vmin=Z.min(), vmax=Z.max()))
    #im_aless=ax3.pcolor(X_aless, Y_aless, Z_aless, norm=colors.LogNorm(vmin=Z_aless.min(), vmax=Z_aless.max()))
    #im_hfls3=ax4.pcolor(X_hfls3, Y_hfls3, Z_hfls3, norm=colors.LogNorm(vmin=Z_hfls3.min(), vmax=Z_hfls3.max()))

    #cbar=fig.colorbar(im, ax=ax2, label='chi_squared (Arp220)')
    #cbar_aless=fig.colorbar(im_aless, ax=ax3, label='chi_squared (ALESS)')
    #cbar_hfls3=fig.colorbar(im_hfls3, ax=ax4, label='chi_squared (HFLS3)')



    #ax2.set_yscale('log',nonposy='clip')
    #ax2.set_ylim([0.1,10.0])  #==========================================================================================================================HERE: change parameter grids

    #ax2.tick_params(width=2, length=16, which='major')
    #ax2.tick_params(width=2, length=5, which='minor')
    #cbar.ax.tick_params(width=2, length=8, which='major')

    #ax2.legend(loc=4)
    #ax2.grid()
    #ax2.set_xlabel('Redshift')  
    ##ax2.set_title('(Arp220)') 
    #ax2.set_ylabel('Normalization factor')  




    #ax3.set_yscale('log',nonposy='clip')
    #ax3.set_ylim([0.1,10.0])  #==========================================================================================================================HERE: change parameter grids

    #ax3.tick_params(width=2, length=16, which='major')
    #ax3.tick_params(width=2, length=5, which='minor')
    #cbar_aless.ax.tick_params(width=2, length=8, which='major')

    #ax3.legend(loc=4)
    #ax3.grid()
    #ax3.set_xlabel('Redshift')  
    ##ax3.set_title('(ALESS)')
    #ax3.set_ylabel('Normalization factor')  





    #ax4.set_yscale('log',nonposy='clip')
    #ax4.set_ylim([0.1,10.0])  #==========================================================================================================================HERE: change parameter grids

    #ax4.tick_params(width=2, length=16, which='major')
    #ax4.tick_params(width=2, length=5, which='minor')
    #cbar_hfls3.ax.tick_params(width=2, length=8, which='major')

    #ax4.legend(loc=4)
    #ax4.grid()
    #ax4.set_xlabel('Redshift')  
    ##ax4.set_title('(HFLS3)')  
    #ax4.set_ylabel('Normalization factor')  







    ##2. SED fitting plot

    #ax1.plot(best_wavelength_micron, best_flux_mjy , color='blue',label='(Arp220) z_phot='+str(best_redshift)+', min_chi^2='+str(min(residual_squared_sum)) , linewidth=1.5)


    #ax1.plot(best_wavelength_micron_aless, best_flux_mjy_aless , color='green',label='(ALESS) z_phot='+str(best_redshift_aless)+', min_chi^2='+str(min(residual_squared_sum_aless)) , linewidth=1.5)
    #ax1.plot(best_wavelength_micron_hfls3, best_flux_mjy_hfls3 , color='purple',label='(HFLS3) z_phot='+str(best_redshift_hfls3)+',min_chi^2='+str(min(residual_squared_sum_hfls3)) , linewidth=1.5)


    #ax1.scatter([250.0, 350.0, 500.0, 850.0],[flux_250_observed[i], flux_350_observed[i], flux_500_observed[i], flux_850_observed[i]], color='red')
    #ax1.errorbar([250.0, 350.0, 500.0, 850.0],[flux_250_observed[i], flux_350_observed[i], flux_500_observed[i], flux_850_observed[i]], yerr=[error_250_observed[i], error_350_observed[i], error_500_observed[i], error_850_observed[i]], color='red',fmt='o', capsize=5, elinewidth=2,markeredgewidth=2)

    #ax1.tick_params(width=2, length=16, which='major')
    #ax1.tick_params(width=2, length=5, which='minor')

    #ax1.set_xscale('log',nonposy='clip')
    #ax1.set_yscale('log',nonposy='clip')

    #ax1.legend(loc=4)
    #ax1.grid()

    #ax1.set_xlabel('Wavelength (micron)')  
    #ax1.set_title(str(source_name[i]))
    #ax1.set_ylabel('flux (mJy)')  
    #plt.rc('font', size=15)


    #fig.set_size_inches(20,10)
    #fig.savefig('./NGP9/'+str(source_name[i])+'_all_templates.jpg')


  







    ##3. Marginalized fitting plot (+SED)



    #fig3 = plt.figure()

    #ax5 = fig3.add_subplot(221)
    #ax6 = fig3.add_subplot(222)

    #ax7 = fig3.add_subplot(223)
    #ax8 = fig3.add_subplot(224)



    #ax5.plot(best_wavelength_micron, best_flux_mjy , color='blue',label='(Arp220) z_phot='+str(best_redshift)+', min_chi^2='+str(min(residual_squared_sum)) , linewidth=1.5)

    #ax5.plot(best_wavelength_micron_aless, best_flux_mjy_aless , color='green',label='(ALESS) z_phot='+str(best_redshift_aless)+', min_chi^2='+str(min(residual_squared_sum_aless)) , linewidth=1.5)
    #ax5.plot(best_wavelength_micron_hfls3, best_flux_mjy_hfls3 , color='purple',label='(HFLS3) z_phot='+str(best_redshift_hfls3)+',min_chi^2='+str(min(residual_squared_sum_hfls3)) , linewidth=1.5)


    #ax5.scatter([250.0, 350.0, 500.0, 850.0],[flux_250_observed[i], flux_350_observed[i], flux_500_observed[i], flux_850_observed[i]], color='red')
    #ax5.errorbar([250.0, 350.0, 500.0, 850.0],[flux_250_observed[i], flux_350_observed[i], flux_500_observed[i], flux_850_observed[i]], yerr=[error_250_observed[i], error_350_observed[i], error_500_observed[i], error_850_observed[i]], color='red',fmt='o', capsize=5, elinewidth=2,markeredgewidth=2)

    #ax5.tick_params(width=2, length=16, which='major')
    #ax5.tick_params(width=2, length=5, which='minor')

    #ax5.set_xscale('log',nonposy='clip')
    #ax5.set_yscale('log',nonposy='clip')

    #ax5.legend(loc=4)
    #ax5.grid()
    #ax5.set_xlabel('Wavelength (micron)')  
    #ax5.set_title(str(source_name[i]))
    #ax5.set_ylabel('flux (mJy)')  
    #plt.rc('font', size=15)






    #X_marginalized, Y_marginalized, Z_marginalized = grid(z_list_space, dely_list_space, likelihood_sum_marginalize)


    #im_marginalized=ax6.pcolor(X_marginalized, Y_marginalized, Z_marginalized)

    #cbar_marginalized=fig.colorbar(im_marginalized, ax=ax6, label='likelihood (marginalized)')

    #ax6.set_yscale('log',nonposy='clip')
    #ax6.set_ylim([0.1,10.0])  #==========================================================================================================================HERE: change parameter grids

    #ax6.tick_params(width=2, length=16, which='major')
    #ax6.tick_params(width=2, length=5, which='minor')
    #cbar.ax.tick_params(width=2, length=8, which='major')

    #ax6.legend(loc=4)
    #ax6.grid()
    #ax6.set_xlabel('Redshift')  
    #ax6.set_title(str(source_name[i])+', marginalized_z='+str(best_redshift_marginalized)+', max_likelihood='+str(max(likelihood_sum_marginalize))) 
    #ax6.set_ylabel('Normalization factor')  





    #ax7.plot(z_list,likelihood_marginalize_over_z, color='blue', linewidth=1.5)
    #ax7.tick_params(width=2, length=16, which='major')
    #ax7.tick_params(width=2, length=5, which='minor')

    ##ax7.set_xscale('log',nonposy='clip')
    #ax7.set_yscale('log',nonposy='clip')
    #ax7.set_ylim(bottom=1.0e-10)

    ##ax7.legend(loc=4)
    #ax7.grid()
    #ax7.set_xlabel('z')  
    ##ax7.set_title(str(source_name[i]))
    #ax7.set_ylabel('P(z|F)')  
    #plt.rc('font', size=15)


    #ax8.plot(dely_list,likelihood_marginalize_over_a, color='blue', linewidth=1.5)
    #ax8.tick_params(width=2, length=16, which='major')
    #ax8.tick_params(width=2, length=5, which='minor')

    #ax8.set_xscale('log',nonposy='clip')
    #ax8.set_yscale('log',nonposy='clip')
    #ax8.set_ylim(bottom=1.0e-10)

    ##ax8.legend(loc=4)
    #ax8.grid()
    #ax8.set_xlabel('normalization factor')  
    ##ax8.set_title(str(source_name[i]))
    #ax8.set_ylabel('P(normalization factor|F)')  
    #plt.rc('font', size=15)







    #fig3.set_size_inches(20,10)
    #fig3.savefig('./NGP9/'+str(source_name[i])+'_marginalized.jpg')  #==========================================================================================================================HERE






print best_redshift_all
print best_dely_all
print best_redshift_all_aless
print best_dely_all_aless
print best_redshift_all_hfls3
print best_dely_all_hfls3
print best_redshift_all_marginalized
print best_dely_all_marginalized


#from astropy.table import Table


#check_Herschel_table = Table()    

#check_Herschel_table['source_name'] = source_name
#check_Herschel_table['RA'] = RA
#check_Herschel_table['Dec'] = Dec

#check_Herschel_table['best_redshift_arp220'] = best_redshift_all
#check_Herschel_table['best_dely_arp220'] = best_dely_all
#check_Herschel_table['min_chi2_arp220'] = min_chi2_all

#check_Herschel_table['best_redshift_aless'] = best_redshift_all_aless
#check_Herschel_table['best_dely_aless'] = best_dely_all_aless
#check_Herschel_table['min_chi2_aless'] = min_chi2_all_aless

#check_Herschel_table['best_redshift_hfls3'] = best_redshift_all_hfls3
#check_Herschel_table['best_dely_hfls3'] = best_dely_all_hfls3
#check_Herschel_table['min_chi2_hfls3'] = min_chi2_all_hfls3

#check_Herschel_table['best_redshift_marginalized'] = best_redshift_all_marginalized
#check_Herschel_table['best_dely_marginalized'] = best_dely_all_marginalized
#check_Herschel_table['max_likelihood_marginalized'] = max_likelihood_marginalized_all

#check_Herschel_table.write('./NGP9/photoz_estimate_cat_NGP9.fits', format='fits')   #=============================================================================================HERE







best_redshift_all_marginalized = np.array(best_redshift_all_marginalized)
best_redshift_all_marginalized=np.array([float(i) for i in best_redshift_all_marginalized])




#==================Plot z-distribution===========================================================================================================================
fig2 = plt.figure()

          
plt.hist(best_redshift_all_marginalized, label='Marginalized', alpha=0.4, bins=np.arange(min(best_redshift_all_marginalized), max(best_redshift_all_marginalized) + 0.5, 0.5))                     


plt.title("All clusters") #===============================================================================================================================HERE
plt.xlabel("Photometric Redshift")
plt.ylabel("Number of Sources")

plt.legend()
plt.rc('font', size=15)
fig2.set_size_inches(20,10)

fig2.savefig('./allcluster_hist.jpg')  #===============================================================================================================================HERE



mock_catalogue=np.column_stack((best_redshift_all_marginalized))

np.savetxt('./allcluster_hist.txt', mock_catalogue, delimiter=' ') #================================================================================HERE





#---------------------Plot multiplied marginalized likelohood----------------------------------------------------------------------------


mock_catalogue_marginalize_z=np.column_stack((z_list,multiply_value_z))
mock_catalogue_marginalize_a=np.column_stack((dely_list,multiply_value_a))

np.savetxt('./allcluster_marginalize_all_z.txt', mock_catalogue_marginalize_z, delimiter=' ') #================================================================================HERE
np.savetxt('./allcluster_marginalize_all_a.txt', mock_catalogue_marginalize_a, delimiter=' ') #================================================================================HERE

fig4=plt.figure()

ax9 = fig4.add_subplot(121)
ax10 = fig4.add_subplot(122)





ax9.plot(z_list,multiply_value_z, color='blue', linewidth=1.5)
ax9.tick_params(width=2, length=16, which='major')
ax9.tick_params(width=2, length=5, which='minor')

##ax9.set_xscale('log',nonposy='clip')
#ax9.set_yscale('log',nonposy='clip')
#ax9.set_ylim(bottom=1.0e-30)
ax9.set_ylim([min(multiply_value_z),max(multiply_value_z)])

#ax9.legend(loc=4)
ax9.grid()
ax9.set_xlabel('z')  
ax9.set_title('All clusters')  #=========================================================================================HERE
ax9.set_ylabel('P(z|F), all sources')  
plt.rc('font', size=15)



ax10.plot(dely_list,multiply_value_a, color='blue', linewidth=1.5)
ax10.tick_params(width=2, length=16, which='major')
ax10.tick_params(width=2, length=5, which='minor')

ax10.set_xscale('log',nonposy='clip')
#ax10.set_yscale('log',nonposy='clip')
#ax10.set_ylim(bottom=1.0e-30)
ax10.set_ylim([min(multiply_value_a),max(multiply_value_a)])

#ax10.legend(loc=4)
ax10.grid()
ax10.set_xlabel('normalization factor')  
ax10.set_title('All clusters')  #=========================================================================================HERE
ax10.set_ylabel('P(normalization factor|F), all sources')  
plt.rc('font', size=15)



fig4.set_size_inches(20,10)
fig4.savefig('./allcluster_marginalize_all.jpg')  #==========================================================================================================================HERE



