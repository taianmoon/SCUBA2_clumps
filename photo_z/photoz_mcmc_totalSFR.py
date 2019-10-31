#import matplotlib.pyplot as plt
#import numpy as np
#import math
#import matplotlib.cm as cm
#from astropy.io import fits
#import astropy.wcs as wcs
#from scipy.interpolate import interp1d
#import pyfits
#import os
#import subprocess


#=======================import the catalogues, get flux and err=============================================

if cluster_list[l]=='PCL1002' or cluster_list[l]=='S2CLS_random':
    names_wht = pyfits.open('./'+cluster_list[l]+'/'+cluster_list[l]+'_herschel_match_3p5sigma_cat_edgesourcedelete_allscuba2source_herschelflux.fits')    #===============================HERE
else:
    names_wht = pyfits.open('../'+cluster_list[l]+'/'+cluster_list[l]+'_herschel_match_3p5sigma_cat_edgesourcedelete_allscuba2source_herschelflux.fits')    #===============================HERE
    #names_wht = pyfits.open('../'+cluster_list[l]+'/test.fits')                            #================================HERE

names_wht_data = names_wht[1].data


if cluster_list[l]=='Bootes1':
    #error_250_observed_pre= names_wht_data.field('et_f250')    #========================================================HERE
    error_350_observed_pre= names_wht_data.field('et_f350')    #========================================================HERE
    error_500_observed_pre= names_wht_data.field('et_f500')    #========================================================HERE

if cluster_list[l]=='EGS' or cluster_list[l]=='Lockman' or cluster_list[l]=='PCL1002' or cluster_list[l]=='S2CLS_random':
    #error_250_observed_pre= names_wht_data.field('et250')    #========================================================HERE
    error_350_observed_pre= names_wht_data.field('et350')    #========================================================HERE
    error_500_observed_pre= names_wht_data.field('et500')    #========================================================HERE

if cluster_list[l]=='G12' or cluster_list[l]=='NGP1' or cluster_list[l]=='NGP2' or cluster_list[l]=='NGP3' or cluster_list[l]=='NGP4' or cluster_list[l]=='NGP5' or cluster_list[l]=='NGP6' or cluster_list[l]=='NGP7' or cluster_list[l]=='NGP8' or cluster_list[l]=='NGP9' or cluster_list[l]=='totalSFR':
    #error_250_observed_pre= names_wht_data.field('e250')    #========================================================HERE
    error_350_observed_pre= names_wht_data.field('e350')    #========================================================HERE
    error_500_observed_pre= names_wht_data.field('e500')    #========================================================HERE

#flux_250_observed_pre= names_wht_data.field('F250')    #========================================================HERE(Boote.s: F250, et_F250; EG.S/Lockma.n: f250, et250; G1.2/NG.P: F250, E250)
#error_250_observed_pre= names_wht_data.field('et_f250')    #========================================================HERE
flux_350_observed_pre= names_wht_data.field('F350')    #========================================================HERE
#error_350_observed_pre= names_wht_data.field('et_f350')    #========================================================HERE
flux_500_observed_pre= names_wht_data.field('F500')    #========================================================HERE
#error_500_observed_pre= names_wht_data.field('et_f500')    #========================================================HERE
flux_850_observed_pre= names_wht_data.field('flux_mJy_perbeam_deboost')

#error_850_observed_pre= names_wht_data.field('flux_error_mJy_deboost_confusion')
error_850_observed_pre= names_wht_data.field('flux_error_mJy_deboost_calibration')
#error_850_observed_pre= names_wht_data.field('test_col')

RA_pre= names_wht_data.field('RA_1')
Dec_pre= names_wht_data.field('Dec_1')
source_name_pre= names_wht_data.field('name')
SN_ratio_pre= names_wht_data.field('SN_ratio_1')
SN_ratio_confusion_pre= names_wht_data.field('SN_ratio_confusion')


#flux_250_observed_pre = np.array(flux_250_observed_pre)
#flux_250_observed_pre=np.array([float(i) for i in flux_250_observed_pre])
#error_250_observed_pre = np.array(error_250_observed_pre)
#error_250_observed_pre=np.array([float(i) for i in error_250_observed_pre])
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
SN_ratio_pre = np.array(SN_ratio_pre)
SN_ratio_pre=np.array([float(i) for i in SN_ratio_pre])
SN_ratio_confusion_pre = np.array(SN_ratio_confusion_pre)
SN_ratio_confusion_pre=np.array([float(i) for i in SN_ratio_confusion_pre])



#flux_250_observed=[]
#error_250_observed=[]
flux_350_observed=[]
error_350_observed=[]
flux_500_observed=[]
error_500_observed=[]
flux_850_observed=[]
error_850_observed=[]
RA=[]
Dec=[]
source_name=[]

for i in range(0, len(flux_350_observed_pre)):

    if  math.isnan(flux_350_observed_pre[i]) == False and SN_ratio_confusion_pre[i]>=3.5:

        #flux_250_observed.append(flux_250_observed_pre[i])
        #error_250_observed.append(error_250_observed_pre[i])
        flux_350_observed.append(flux_350_observed_pre[i])
        error_350_observed.append(error_350_observed_pre[i])
        flux_500_observed.append(flux_500_observed_pre[i])
        error_500_observed.append(error_500_observed_pre[i])
        flux_850_observed.append(flux_850_observed_pre[i])
        error_850_observed.append(error_850_observed_pre[i])
        RA.append(RA_pre[i])
        Dec.append(Dec_pre[i])
        source_name.append(source_name_pre[i])


#print source_name


#flux_250_observed = np.array(flux_250_observed)
#flux_250_observed=np.array([float(i) for i in flux_250_observed])
#error_250_observed = np.array(error_250_observed)
#error_250_observed=np.array([float(i) for i in error_250_observed])
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
if cluster_list[l]=='G12' or cluster_list[l]=='NGP1' or cluster_list[l]=='NGP2' or cluster_list[l]=='NGP3' or cluster_list[l]=='NGP4' or cluster_list[l]=='NGP5' or cluster_list[l]=='NGP6' or cluster_list[l]=='NGP7' or cluster_list[l]=='NGP8' or cluster_list[l]=='NGP9':
    #flux_250_observed=flux_250_observed*1000.0
    #error_250_observed=error_250_observed*1000.0
    flux_350_observed=flux_350_observed*1000.0
    error_350_observed=error_350_observed*1000.0
    flux_500_observed=flux_500_observed*1000.0
    error_500_observed=error_500_observed*1000.0
#----------------------------------------------------------------------------------------




##=====================import template(s) redshift it, interpolate the model values at the correspoinding bands, and create parameter space (z and L --> delx and dely)==================
#
#
#
##template_file = open("./arp220_template_donley2007_apj_660_167.txt", "r")    #======================================================================HERE: change template
#template_file = open("./aless_average_seds.dat", "r")    #===========================================================================================HERE: change template
##template_file = open("./HFLS3.txt", "r")    #=======================================================================================================HERE: change template
##template_file = open("./cosmic_eyelash_sed_lapi.txt", "r")    #======================================================================HERE: change template
#
#
#lines = template_file.readlines()
#template_file.close()
#
#response_curve=[]
#
#for i in range(0, len(lines)):
#    separated_lines=lines[i].split() 
#    response_curve.append(separated_lines)
#
#
#response_curve = np.array(response_curve)
#wavelength_micron=response_curve[:,0]
#flux_mjy=response_curve[:,1]
#
#
#wavelength_micron = np.array(wavelength_micron)
#flux_mjy = np.array(flux_mjy)
#wavelength_micron=np.array([float(i) for i in wavelength_micron])  #(wavelength in angstrom)
#flux_mjy=np.array([float(i) for i in flux_mjy])  #(wavelength in angstrom)
#
#
##+++++++++++++++++++++For HFLS3+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=====================================HERE: for HFLS3 template
##wavelength_micron=wavelength_micron/7.34  #redshift back to restframe (z=6.34)
##wavelength_micron=wavelength_micron*1.0e6  #convert from m to micron
##flux_mjy=flux_mjy*1000.0  #convert from Jy to mJy
#
#
##+++++++++++++++++++++For Cosmic Eyelash+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=====================================HERE: for Cosmic Eyelash template
#
##flux_mjy_loop=[]
##for j in range(0, len(wavelength_micron)):
#
##    flux_mjy_loop.append(1.0e-14*1000.0*(1.0e26)*flux_mjy[j]*(wavelength_micron[j]*1.0e-6)/2.99792458e8)
#
#
##flux_mjy=flux_mjy_loop
#
##flux_mjy = np.array(flux_mjy)
##flux_mjy=np.array([float(i) for i in flux_mjy])  #(wavelength in angstrom)
#
##plt.plot(wavelength_micron, flux_mjy)
##plt.xscale('log',nonposy='clip')
##plt.yscale('log',nonposy='clip')
##plt.show()
#
#============================================Doing MCMC==============================================================






#model data

model_x=wavelength_micron  
model_y=flux_mjy  


model_x = np.array(model_x)
model_x=np.array([float(i) for i in model_x])  #(wavelength in angstrom)
model_y = np.array(model_y)
model_y=np.array([float(i) for i in model_y])




def model_flux(z_1, a_1, wanted_wavelength):
    model_x_z = model_x*(1.0+z_1)
    model_y_a = model_y*pow(10.0,a_1)
    #f=np.array([np.interp(wanted_wavelength[0], model_x_z, model_y_a), np.interp(wanted_wavelength[1], model_x_z, model_y_a)])
    return np.interp(wanted_wavelength, model_x_z, model_y_a)
    #return f


#assuming prior=1


def lnlike(theta_2, x_2, y_2, yerr_2):
    z, a = theta_2
    #model = model_flux(z, a, 250.0)
    model = model_flux(z, a, x_2)
    inv_sigma2 = 1.0/(yerr_2**2)
    return -0.5*(np.sum((y_2-model)**2*inv_sigma2))




def lnprior(theta_1):
    z, a = theta_1
    if 0.0 < z < 15.0 and -2.0 < a < 2.0:  #==============================================================================HERE: parameter constraints
        return 0.0
    return -np.inf



def lnprob(theta, x_3, y_3, yerr_3):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x_3, y_3, yerr_3)






# observed data.


x=[350.0, 550.0, 850.0]  #======================================190730 new!!
x = np.array(x)
x=np.array([float(i) for i in x])  #(wavelength in angstrom)


#-----------------import min chi^2 values from grid method---------------------------------------(temporarily not used)

#names_chi2 = pyfits.open('./'+cluster_list[l]+'/photoz_estimate_cat_'+cluster_list[l]+'.fits')                            #========================================================HERE
#names_chi2_data = names_chi2[1].data


#best_z_chi2= names_chi2_data.field('best_redshift_aless')    #============================================================================================HERE: change template (no change for cosmic eyelash)
#best_a_chi2= names_chi2_data.field('best_dely_aless')    #================================================================================================HERE: change template (no change for cosmic eyelash)

#------------------------------------------------------------------------------------------------

#--------------------------------MCMC preparations----------------------------------------------

steps=5000   #===========================================================================================================================================HERE: Change steps
steps_to_be_cut=1000    #=================================================================================================================================HERE: Change steps to be cut
ndim, nwalkers = 2, 20  #==============================================================================================================================HERE: Change nwalkers

#import emcee
#import corner
#import random



#source_combine_z=[]
#source_combine_a=[]
#multiplier_z=1.0
#multiplier_a=1.0


lower_limit_z_table=[]
fitted_median_z_table=[]
higher_limit_z_table=[]
lower_limit_a_table=[]
fitted_median_a_table=[]
higher_limit_a_table=[]


for j in range(0, len(flux_350_observed)):  #looping different sources
#for j in range(0, 2):  #looping different sources

    
    y=[flux_350_observed[j], flux_500_observed[j], flux_850_observed[j]]
    yerr=[error_350_observed[j], error_500_observed[j], error_850_observed[j]]


    #print y
    #print yerr
    #print 'ccccccccccccccccccccccccccccccccccccccc'



    y = np.array(y)
    y=np.array([float(i) for i in y])  #(wavelength in angstrom)
    yerr = np.array(yerr)
    yerr=np.array([float(i) for i in yerr])  #(wavelength in angstrom)





    #-----------------import min chi^2 values from grid method (continued)---------------------------------------(temporarily not used)

    #best_chi2=[best_z_chi2[j], best_a_chi2[j]]

    #best_chi2 = np.array(best_chi2)
    #best_chi2=np.array([float(i) for i in best_chi2])  #(wavelength in angstrom)

    #print best_chi2
    #print 'uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu'

    #------------------------------------------------------------------------------------------------


    #import random
    #pos = [best_chi2 + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]    #============================================================================HERE: Change starting points
    pos=[[pow(10.0,np.random.uniform(-2,1)), random.uniform(-2.0, 2.0)] for i in range(nwalkers)]
    #pos=[[pow(10.0,np.random.uniform(-2,1)), pow(10.0,np.random.uniform(-2,1))] for i in range(nwalkers)]

    #print pos
    #print 'oooooooooooooooooooooooooooooooo'



    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))



    sampler.run_mcmc(pos, steps)



    samples = sampler.chain[:, steps_to_be_cut:, :].reshape((-1, ndim))

    np.save('./MCMC_results/'+cluster_list[l]+'/'+str(source_name[j])+'_mcmc_values_'+str(template_name[w])+'.npy', samples)    #===========================================================HERE (columns: z, a for each step/walker)
    #np.savetxt('./'+str(source_name[j])+'_mcmc_values.cat', samples, delimiter=' ') 

    print np.shape(samples)
    print '#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'

    #fig = corner.corner(samples, labels=["$z$", "$normalization factor$"], truths=[z_true, a_true])
    fig = corner.corner(samples, labels=["$z$", "$log(a)$"],
                        quantiles=[0.16, 0.5, 0.84],
                        show_titles=True, title_kwargs={"fontsize": 12})
    fig.savefig("./MCMC_results/"+cluster_list[l]+"/triangle_"+template_name[w]+"_"+str(source_name[j])+".png")  

    plt.close(fig)
    #fig.show()


    #for kk in range(0, len(samples[:,0])):  #looping different sources

        #source_combine_z.append(samples[kk,0])
        #source_combine_a.append(samples[kk,1])




    #-----------plot chain plot-----------------------------------------------------------

    plt.close('all')
    plt.figure()
    for k in range(0, nwalkers):

        plt.subplot(ndim, 1, 1)
        plt.plot(np.arange(steps)+1, sampler.chain[k,:,0], color='k', alpha = 0.4)
        plt.ylabel('z')

        plt.subplot(ndim, 1, 2)
        plt.plot(np.arange(steps)+1, sampler.chain[k,:,1], color='k', alpha = 0.4)
        plt.ylabel('normalization factor')

    plt.xlabel('step number')
    plt.savefig("./MCMC_results/"+cluster_list[l]+"/chain_"+template_name[w]+"_"+str(source_name[j])+".png") 
    plt.close(fig)


    #-----------------------Plot SED--------------------------------------------------------
    
    lower_limit_z, fitted_median_z, higher_limit_z = corner.quantile(samples[:,0], q=[0.16, 0.5, 0.84], weights=None)
    print source_name[j]
    print fitted_median_z
    print lower_limit_z
    print higher_limit_z


    lower_limit_a, fitted_median_a, higher_limit_a = corner.quantile(samples[:,1], q=[0.16, 0.5, 0.84], weights=None)
    print fitted_median_a
    print lower_limit_a
    print higher_limit_a


    lower_limit_z_table.append(lower_limit_z)
    fitted_median_z_table.append(fitted_median_z)
    higher_limit_z_table.append(higher_limit_z)

    lower_limit_a_table.append(lower_limit_a)
    fitted_median_a_table.append(fitted_median_a)
    higher_limit_a_table.append(higher_limit_a)

    model_x_plot_median = model_x*(1.0+fitted_median_z) 
    model_y_plot_median = model_y*pow(10.0,fitted_median_a)

    model_x_plot_low = model_x*(1.0+lower_limit_z) 
    model_y_plot_low = model_y*pow(10.0,lower_limit_a)

    model_x_plot_high = model_x*(1.0+higher_limit_z) 
    model_y_plot_high = model_y*pow(10.0,higher_limit_a)



    plt.close('all')
    plt.figure()
    plt.plot(model_x_plot_median, model_y_plot_median , color='blue' , linewidth=1.5)
    plt.plot(model_x_plot_low, model_y_plot_low , color='black', linewidth=1.0, linestyle='--')
    plt.plot(model_x_plot_high, model_y_plot_high , color='black', linewidth=1.0, linestyle='--')
    plt.plot(model_x, model_y , color='black' , linewidth=1.0, linestyle=':')


    plt.scatter(x, y, color='red')
    plt.errorbar(x,y, yerr=yerr, color='red',fmt='o', capsize=5, elinewidth=2,markeredgewidth=2)

    plt.tick_params(width=2, length=16, which='major')
    plt.tick_params(width=2, length=5, which='minor')

    plt.xscale('log',nonposy='clip')
    plt.yscale('log',nonposy='clip')

    #plt.legend(loc=4)
    plt.grid()
    plt.xlabel('Wavelength (micron)')  
    #plt.set_title(str(source_name[i]))
    plt.ylabel('flux (mJy)')  
    plt.rc('font', size=15)
    plt.savefig("./MCMC_results/"+cluster_list[l]+"/sed_"+template_name[w]+"_"+str(source_name[j])+".png")   
    plt.close(fig)




#source_combine_z = np.array(source_combine_z)
#source_combine_z=np.array([float(i) for i in source_combine_z])
#source_combine_a = np.array(source_combine_a)
#source_combine_a=np.array([float(i) for i in source_combine_a])

#mock_catalogue=np.column_stack((source_combine_z, source_combine_a))
#np.savetxt('./Bootes1_mcmc_values.cat', mock_catalogue, delimiter=' ') #================================================================================HERE (columns: z, a for each step/walker)


#print source_combine_z
#print source_combine_a
#print np.shape(mock_catalogue)
#print '#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'



#plt.close('all')
#plt.figure()
#fig = corner.corner(mock_catalogue, labels=["$z$", "$a$"],
#                    quantiles=[0.16, 0.5, 0.84],
#                    show_titles=True, title_kwargs={"fontsize": 12})
#fig.savefig("./allsource_triangle_arp220_Bootes1.png")  #==========================================================================================HERE: change template

#plt.close(fig)

print 'dddddddddddddddddddddddddddddd'
print lower_limit_z_table
print fitted_median_z_table
print higher_limit_z_table
print lower_limit_a_table
print fitted_median_a_table
print higher_limit_a_table




check_Herschel_table = Table.read('./'+cluster_list[l]+'/photoz_estimate_cat_'+cluster_list[l]+'.fits', format='fits')    

check_Herschel_table['z_'+str(template_name[w])+'_mcmc'] = fitted_median_z_table
check_Herschel_table['z_'+str(template_name[w])+'_mcmc_err_low'] = lower_limit_z_table
check_Herschel_table['z_'+str(template_name[w])+'_mcmc_err_high'] = higher_limit_z_table
check_Herschel_table['a_'+str(template_name[w])+'_mcmc'] = fitted_median_a_table
check_Herschel_table['a_'+str(template_name[w])+'_mcmc_err_low'] = lower_limit_a_table
check_Herschel_table['a_'+str(template_name[w])+'_mcmc_err_high'] = higher_limit_a_table

check_Herschel_table.write('./MCMC_results/'+cluster_list[l]+'/photoz_estimate_cat_'+cluster_list[l]+'_'+str(template_name[w])+'.fits', format='fits')   







