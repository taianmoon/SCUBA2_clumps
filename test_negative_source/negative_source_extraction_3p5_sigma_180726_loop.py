from astropy.io import fits
import astropy.wcs as wcs
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import scipy.interpolate



cluster_list=['Bootes1', 'EGS', 'Lockman', 'G12', 'NGP1', 'NGP2', 'NGP3', 'NGP4', 'NGP5', 'NGP6', 'NGP7', 'NGP8', 'NGP9', 'S2CLS']
#cluster_list=['Bootes1', 'EGS', 'Lockman', 'G12', 'NGP1', 'NGP2', 'NGP3', 'NGP4', 'NGP5', 'NGP6', 'NGP7', 'NGP8', 'NGP9']



sn_want=3.5  #========================================================================HERE: change threshold for sn ratio






for l in range(0, len(cluster_list)):

    print cluster_list[l]
    #execfile("./negative_source_extraction_3p5_sigma_180726.py")  #commented out after negative catalogues (without edgesourcedelete) are built





number_of_true_sources=[]  #==================================================================HERE: must match the cluster list
fig = plt.figure(figsize=(16.0, 10.0))


for kk in range(0, len(cluster_list)):

    print cluster_list[kk]

    template_file = open("./for_code/"+cluster_list[kk]+"_3p5sigma_negative_edgesourcedelete.cat", "r")
    lines = template_file.readlines()
    template_file.close()
    response_curve=[]
    for i in range(0, len(lines)):
        separated_lines=lines[i].split() 
        response_curve.append(separated_lines)
    response_curve = np.array(response_curve)
    sn_ratio=response_curve[:,4]
    flux_mjy_negative=response_curve[:,6]
    fluxerr_mjy=response_curve[:,7]

    counts_array=np.arange(1.0,len(sn_ratio)+1.0)

    sn_ratio = np.array(sn_ratio)
    flux_mjy_negative = np.array(flux_mjy_negative)
    fluxerr_mjy = np.array(fluxerr_mjy)
    counts_array = np.array(counts_array)
    sn_ratio=np.array([float(i) for i in sn_ratio])  #(wavelength in angstrom)
    flux_mjy_negative=np.array([float(i) for i in flux_mjy_negative])  #(wavelength in angstrom)
    fluxerr_mjy=np.array([float(i) for i in fluxerr_mjy])  #(wavelength in angstrom)
    counts_array=np.array([float(i) for i in counts_array])  #(wavelength in angstrom)

    
    sn_ratio=sn_ratio[::-1]
    counts_array=counts_array[::-1]
    #print sn_ratio
    #print counts_array





    if cluster_list[kk]=='S2CLS':
        #template_file = open("./catalogues/S2CLS_3p5sigma.cat", "r")    #==================================================================================================HERE
        template_file_true = open("../"+cluster_list[kk]+"/S2CLS_3p5sigma_deboost.cat", "r")    #==================================================================================================HERE
    else:
        #template_file = open("./catalogues/"+cluster_list[l]+"_3p5sigma_edgesourcedelete.cat", "r")    #======================================================================================HERE
        template_file_true = open("../"+cluster_list[kk]+"/"+cluster_list[kk]+"_3p5sigma_edgesourcedelete_deboost.cat", "r")    #=====================================================================HERE


    #template_file_true = open("../number_counts/catalogues/"+cluster_list[kk]+"_3p5sigma_edgesourcedelete.cat", "r")
    #lines_true = template_file_true.readlines()
    lines_true = template_file_true.readlines()[1:]
    template_file_true.close()
    response_curve_true=[]
    for i in range(0, len(lines_true)):
        separated_lines_true=lines_true[i].split() 
        response_curve_true.append(separated_lines_true)
    response_curve_true = np.array(response_curve_true)
    #sn_ratio_true=response_curve_true[:,4]
    #flux_mjy_true=response_curve_true[:,6]
    #fluxerr_mjy_true=response_curve_true[:,7]


    if cluster_list[kk]=='S2CLS':
        flux_mjy_true=response_curve_true[:,8]      #=========================================================HERE: 7 for Todd's cluster, 6 for our clusters
        fluxerr_mjy_true=response_curve_true[:,9]    #=========================================================HERE: 8 for Todd's cluster, 7 for our clusters
        sn_ratio_true=response_curve_true[:,4]     #=========================================================HERE: 5 for Todd's cluster, 4 for our clusters
    else:
        flux_mjy_true=response_curve_true[:,9]      #=========================================================HERE: 7 for Todd's cluster, 6 for our clusters
        fluxerr_mjy_true=response_curve_true[:,10]    #=========================================================HERE: 8 for Todd's cluster, 7 for our clusters
        sn_ratio_true=response_curve_true[:,5]     #=========================================================HERE: 5 for Todd's cluster, 4 for our clusters





    counts_array_true=np.arange(1.0,len(sn_ratio_true)+1.0)

    sn_ratio_true = np.array(sn_ratio_true)
    flux_mjy_true = np.array(flux_mjy_true)
    fluxerr_mjy_true = np.array(fluxerr_mjy_true)
    counts_array_true = np.array(counts_array_true)
    sn_ratio_true=np.array([float(i) for i in sn_ratio_true])  #(wavelength in angstrom)
    flux_mjy_true=np.array([float(i) for i in flux_mjy_true])  #(wavelength in angstrom)
    fluxerr_mjy_true=np.array([float(i) for i in fluxerr_mjy_true])  #(wavelength in angstrom)
    counts_array_true=np.array([float(i) for i in counts_array_true])  #(wavelength in angstrom)

    
    sn_ratio_true=sn_ratio_true[::-1]
    counts_array_true=counts_array_true[::-1]
    #print sn_ratio_true
    #print counts_array_true





    #===================Calculate percentage of negative sources===============================================   

    check_sn=[3.5, 4.0, 4.5, 5.0, 5.5]



    #print sn_ratio
    #print counts_array

    if len(sn_ratio)!=1:
        y_interp_neg = scipy.interpolate.interp1d(x=sn_ratio, y=counts_array, bounds_error=False, fill_value=(counts_array[0], 0.0))
    else:
        sn_ratio_new=[0.0, sn_ratio[0]]
        counts_array_new=[counts_array[0], counts_array[0]]
        y_interp_neg = scipy.interpolate.interp1d(x=sn_ratio_new, y=counts_array_new, bounds_error=False, fill_value=(counts_array_new[0], 0.0))

    #print y_interp_neg(3.5)



    if len(sn_ratio_true)!=1:
        y_interp_true = scipy.interpolate.interp1d(x=sn_ratio_true, y=counts_array_true, bounds_error=False, fill_value=(counts_array_true[0], 0.0))
    else:
        sn_ratio_true_new=[0.0, sn_ratio_true[0]]
        counts_array_true_new=[counts_array_true[0], counts_array_true[0]]
        y_interp_true = scipy.interpolate.interp1d(x=sn_ratio_true_new, y=counts_array_true_new, bounds_error=False, fill_value=(counts_array_true_new[0], 0.0))

    #print y_interp_true(3.5)


    fraction_true2neg=[]
    for k in range(0, len(check_sn)):
        fraction_true2neg.append(y_interp_neg(check_sn[k])/y_interp_true(check_sn[k]))   #===============================================Check here!
    
    print 'vvvvvvvvvvvvvvvvvvvvvvv'
    print fraction_true2neg





    #==========================================Plotting===========================================================


    ax = fig.add_subplot(3,5,kk+1)

    ax.plot(sn_ratio, counts_array,  color='blue', linewidth=1.5, label='Negative')
    ax.scatter(sn_ratio, counts_array,  color='blue', s=10.0)
    ax.plot(sn_ratio_true, counts_array_true,  color='red', linewidth=1.5, label='Positive')
    ax.scatter(sn_ratio_true, counts_array_true,  color='red', s=10.0)
    ax.tick_params(width=2, length=16, which='major')
    if kk==len(cluster_list)-1:
        ax.legend()
    ax.grid()
    ax.set_xlabel('S/N')  
    ax.set_title(cluster_list[kk])
    ax.set_ylabel('number of sources')
    ax.set_xlim([3,6])
    ax.set_xticks(np.arange(3, 6, step=0.5))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    ax2 = ax.twinx()
    ax2.plot(check_sn, fraction_true2neg, color='black', linewidth=1.5, linestyle='--')
    ax2.scatter(check_sn, fraction_true2neg, color='black', s=10.0)
    ax2.set_ylabel('percentage of negative sources')
    ax2.set_xlim([3,6])
    ax2.set_xticks(np.arange(3, 6, step=0.5))
    ax2.set_ylim([0.0,1.0])
    ax2.set_yticks(np.arange(0.0, 1.0, step=0.1))


    plt.rc('font', size=12)


plt.show()



