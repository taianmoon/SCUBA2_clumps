import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm
from astropy.io import fits
import pyfits

cluster_list=['Bootes1', 'EGS', 'Lockman', 'G12', 'NGP1', 'NGP2', 'NGP3', 'NGP4', 'NGP5', 'NGP6', 'NGP7', 'NGP8', 'NGP9']



fig = plt.figure(figsize=(16.0, 10.0))
fig2 = plt.figure(figsize=(16.0, 10.0))
one2one=np.arange(-0.5, 11.5, 1.0)

for l in range(0, len(cluster_list)):

    print cluster_list[l]

    names_wht_aless = pyfits.open('./'+cluster_list[l]+'/photoz_estimate_cat_'+cluster_list[l]+'_aless.fits')                            #========================================================HERE
    names_wht_data_aless = names_wht_aless[1].data


    z_gridding_aless= names_wht_data_aless.field('best_redshift_aless')
    a_gridding_aless= names_wht_data_aless.field('best_dely_aless')
    z_gridding_arp220= names_wht_data_aless.field('best_redshift_arp220')
    a_gridding_arp220= names_wht_data_aless.field('best_dely_arp220')
    z_gridding_hfls3= names_wht_data_aless.field('best_redshift_hfls3')
    a_gridding_hfls3= names_wht_data_aless.field('best_dely_hfls3')
    z_gridding_eyelash= names_wht_data_aless.field('best_redshift_eyelash')
    a_gridding_eyelash= names_wht_data_aless.field('best_dely_eyelash')

    z_mcmc_aless= names_wht_data_aless.field('z_aless_mcmc')
    z_mcmc_aless_err_low= names_wht_data_aless.field('z_aless_mcmc_err_low')
    z_mcmc_aless_err_high= names_wht_data_aless.field('z_aless_mcmc_err_high')
    a_mcmc_aless= names_wht_data_aless.field('a_aless_mcmc')
    a_mcmc_aless_err_low= names_wht_data_aless.field('a_aless_mcmc_err_low')
    a_mcmc_aless_err_high= names_wht_data_aless.field('a_aless_mcmc_err_high')

    names_wht_arp220 = pyfits.open('./'+cluster_list[l]+'/photoz_estimate_cat_'+cluster_list[l]+'_arp220.fits')                            #========================================================HERE
    names_wht_data_arp220 = names_wht_arp220[1].data

    names_wht_hfls3 = pyfits.open('./'+cluster_list[l]+'/photoz_estimate_cat_'+cluster_list[l]+'_hfls3.fits')                            #========================================================HERE
    names_wht_data_hfls3 = names_wht_hfls3[1].data

    names_wht_eyelash = pyfits.open('./'+cluster_list[l]+'/photoz_estimate_cat_'+cluster_list[l]+'_eyelash.fits')                            #========================================================HERE
    names_wht_data_eyelash = names_wht_eyelash[1].data

    z_mcmc_arp220= names_wht_data_arp220.field('z_arp220_mcmc')
    z_mcmc_arp220_err_low= names_wht_data_arp220.field('z_arp220_mcmc_err_low')
    z_mcmc_arp220_err_high= names_wht_data_arp220.field('z_arp220_mcmc_err_high')
    a_mcmc_arp220= names_wht_data_arp220.field('a_arp220_mcmc')
    a_mcmc_arp220_err_low= names_wht_data_arp220.field('a_arp220_mcmc_err_low')
    a_mcmc_arp220_err_high= names_wht_data_arp220.field('a_arp220_mcmc_err_high')

    z_mcmc_hfls3= names_wht_data_hfls3.field('z_hfls3_mcmc')
    z_mcmc_hfls3_err_low= names_wht_data_hfls3.field('z_hfls3_mcmc_err_low')
    z_mcmc_hfls3_err_high= names_wht_data_hfls3.field('z_hfls3_mcmc_err_high')
    a_mcmc_hfls3= names_wht_data_hfls3.field('a_hfls3_mcmc')
    a_mcmc_hfls3_err_low= names_wht_data_hfls3.field('a_hfls3_mcmc_err_low')
    a_mcmc_hfls3_err_high= names_wht_data_hfls3.field('a_hfls3_mcmc_err_high')

    z_mcmc_eyelash= names_wht_data_eyelash.field('z_eyelash_mcmc')
    z_mcmc_eyelash_err_low= names_wht_data_eyelash.field('z_eyelash_mcmc_err_low')
    z_mcmc_eyelash_err_high= names_wht_data_eyelash.field('z_eyelash_mcmc_err_high')
    a_mcmc_eyelash= names_wht_data_eyelash.field('a_eyelash_mcmc')
    a_mcmc_eyelash_err_low= names_wht_data_eyelash.field('a_eyelash_mcmc_err_low')
    a_mcmc_eyelash_err_high= names_wht_data_eyelash.field('a_eyelash_mcmc_err_high')


    ax1 = fig.add_subplot(3,5,l+1)

    ax1.plot(one2one, one2one, color='black', linestyle='--', linewidth=1.5)

    ax1.scatter(z_mcmc_aless, z_gridding_aless,  color='blue', linewidth=1.5, s=40.0, label=cluster_list[l]) #blue: ALESS
    #asymmetric_error_ax1_aless = [z_mcmc_aless_err_low, z_mcmc_aless_err_high]
    #ax1.errorbar(z_mcmc_aless, z_gridding_aless, xerr=asymmetric_error_ax1_aless, color='blue')
    ax1.scatter(z_mcmc_arp220, z_gridding_arp220,  color='red', linewidth=1.5, s=40.0) #red: arp220
    ax1.scatter(z_mcmc_hfls3, z_gridding_hfls3,  color='gray', linewidth=1.5, s=40.0) #gray: hfls3
    ax1.scatter(z_mcmc_eyelash, z_gridding_eyelash,  color='green', linewidth=1.5, s=40.0) #green: eyelash

    ax1.set_xlim(0.0, 10.0)
    ax1.set_ylim(0.0, 10.0)
    ax1.grid()
    ax1.set_xlabel(r'$z_{MCMC}$')  
    ax1.set_ylabel(r'$z_{gridding}$')  

    #plt.tick_params(width=2, length=16, which='major')
    #plt.tick_params(width=2, length=5, which='minor')



    #ax1.set_xscale('log',nonposy='clip')
    #ax1.set_yscale('log',nonposy='clip')


    #ax1.legend(handlelength=0, handletextpad=0, fancybox=True, loc='lower left')
    leg = ax1.legend(handlelength=0, handletextpad=0, fancybox=True, loc='upper left')
    #print leg
    for item in leg.legendHandles:
        item.set_visible(False)
    leg.get_frame().set_linewidth(0.0)
    leg.get_frame().set_alpha(0.75)

    plt.rc('font', size=15)

    #ax1.tick_params(axis='y', which='both', bottom=False, top=False)
    #ax1.set_yticks([])





    ax2 = fig2.add_subplot(3,5,l+1)

    #ax2.plot(one2one, one2one, color='black', linestyle='--', linewidth=1.5)

    ax2.scatter(a_mcmc_aless, a_gridding_aless,  color='blue', linewidth=1.5, s=40.0, label=cluster_list[l]) #blue: ALESS
    ax2.scatter(a_mcmc_arp220, a_gridding_arp220,  color='red', linewidth=1.5, s=40.0) #red: arp220
    ax2.scatter(a_mcmc_hfls3, a_gridding_hfls3,  color='gray', linewidth=1.5, s=40.0) #gray: hfls3
    ax2.scatter(a_mcmc_eyelash, a_gridding_eyelash,  color='green', linewidth=1.5, s=40.0) #green: eyelash

    ax2.set_xlim(1e-2, 1e2)
    ax2.set_ylim(1e-2, 1e2)
    ax2.grid()
    ax2.set_xlabel(r'$a_{MCMC}$')  
    ax2.set_ylabel(r'$a_{gridding}$')  

    ax2.set_xscale('log',nonposy='clip')
    ax2.set_yscale('log',nonposy='clip')

    #ax2.legend(handlelength=0, handletextpad=0, fancybox=True, loc='lower left')
    leg2 = ax2.legend(handlelength=0, handletextpad=0, fancybox=True, loc='upper left')
    #print leg2
    for item in leg2.legendHandles:
        item.set_visible(False)
    leg2.get_frame().set_linewidth(0.0)
    leg2.get_frame().set_alpha(0.75)

    plt.rc('font', size=15)

    #ax2.tick_params(axis='y', which='both', bottom=False, top=False)
    #ax2.set_yticks([])














plt.show()








