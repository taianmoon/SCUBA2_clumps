from astropy.io import fits
import astropy.wcs as wcs
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import scipy.interpolate



#cluster_list=['Planck18p735', 'PLCK_DU_G045.7-41.2', 'PLCK_DU_G059.1-67.1', 'PLCK_DU_G073.4-57.5', 'PLCK_G006.1+61.8', 'PLCK_G009.8+72.6', 'PLCK_G056.7+62.6', 'PLCK_G068.3+31.9', 'PLCK_G075.1+33.2', 'PLCK_G077.7+32.6', 'PLCK_G082.5+38.4']

#cluster_list=[ 'PLCK_G093.6+55.9', 'PLCK_G144.1+81.0', 'PLCK_G160.7+41.0', 'PLCK_G162.1-59.3', 'PLCK_G165.8+45.3' , 'PLCK_G177.0+35.9', 'PLCK_G179.3+50.7', 'PLCK_G186.3-72.7', 'PLCK_G186.6+66.7', 'PLCK_G188.6-68.9', 'PLCK_G191.3+62.0', 'PLCK_G191.8-83.4']

cluster_list=['PLCK_G213.0+65.9', 'PLCK_G223.9+41.2', 'PLCK_G328.9+71.4', 'PLCK_G84.0-71.5', 'PLCK_HZ_G067.2-63.8', 'PLCK_HZ_G103.1-73.6', 'PLCK_HZ_G106.8-83.3', 'PLCK_HZ_G119.4-76.6', 'PLCK_HZ_G132.6-81.1', 'PLCK_HZ_G171.1-78.7', 'PLCK_HZ_G214.1+48.3']

#cluster_list=[]

#clisters that have no source >=3.5 (SN_confusion): 'Planck18p194', 'PLCK_G078.9+48.2', 'PLCK_G083.3+51.0', 'PLCK_G091.9+43.0', 'PLCK_G173.8+59.3', 'PLCK_G201.1+50.7',  'PLCK_HZ_G038.0-51.5',  'PLCK_HZ_G173.9+57.0',  'PLCK_HZ_G176.6+59.0'


sn_want=3.5  #========================================================================HERE: change threshold for sn ratio






for l in range(0, len(cluster_list)):

    print cluster_list[l]
    #execfile("./negative_source_extraction_3p5_sigma_180726.py")  #commented out after negative catalogues (without edgesourcedelete) are built



print '-----------------------Good, now comment out the last line and update the catalogue (edge remove) and run this code again----------------------'

number_of_true_sources=[]  #==================================================================HERE: must match the cluster list
fig = plt.figure(figsize=(16.0, 10.0))


for kk in range(0, len(cluster_list)):

    print cluster_list[kk]

    template_file = open("./for_code/"+cluster_list[kk]+"_3p5sigma_negative.cat", "r")
    lines = template_file.readlines()[1:]
    template_file.close()
    response_curve=[]
    for i in range(0, len(lines)):
        separated_lines=lines[i].split() 
        response_curve.append(separated_lines)
    response_curve = np.array(response_curve)
    sn_ratio_pre=response_curve[:,5]
    flux_mjy_negative_pre=response_curve[:,9]  #debooosted
    fluxerr_mjy_pre=response_curve[:,14]  #deboosted and calibration
    sn_ratio_confusion_pre=response_curve[:,12]


    sn_ratio_pre = np.array(sn_ratio_pre)
    sn_ratio_confusion_pre = np.array(sn_ratio_confusion_pre)
    flux_mjy_negative_pre = np.array(flux_mjy_negative_pre)
    fluxerr_mjy_pre = np.array(fluxerr_mjy_pre)
    sn_ratio_pre=np.array([float(i) for i in sn_ratio_pre])  #(wavelength in angstrom)
    sn_ratio_confusion_pre=np.array([float(i) for i in sn_ratio_confusion_pre])  #(wavelength in angstrom)
    flux_mjy_negative_pre=np.array([float(i) for i in flux_mjy_negative_pre])  #(wavelength in angstrom)
    fluxerr_mjy_pre=np.array([float(i) for i in fluxerr_mjy_pre])  #(wavelength in angstrom)


    sn_ratio=[]
    flux_mjy_negative=[]
    fluxerr_mjy=[]
    sn_ratio_confusion=[]
    for rc in range(0, len(flux_mjy_negative_pre)):
        if abs(sn_ratio_confusion_pre[rc]) >= 3.5:
            sn_ratio.append(abs(sn_ratio_pre[rc]))
            flux_mjy_negative.append(abs(flux_mjy_negative_pre[rc]))
            fluxerr_mjy.append(abs(fluxerr_mjy_pre[rc]))
            sn_ratio_confusion.append(abs(sn_ratio_confusion_pre[rc]))


    counts_array=np.arange(1.0,len(sn_ratio)+1.0)

    counts_array = np.array(counts_array)
    counts_array=np.array([float(i) for i in counts_array])  #(wavelength in angstrom)

    print counts_array
    print sn_ratio
    print 'ssssssssssssss'




    sn_ratio=sn_ratio[::-1]
    counts_array=counts_array[::-1]
    #print sn_ratio
    #print counts_array





    if cluster_list[kk]=='S2CLS':
        #template_file = open("./catalogues/S2CLS_3p5sigma.cat", "r")    #==================================================================================================HERE
        template_file_true = open("../"+cluster_list[kk]+"/S2CLS_3p5sigma_deboost.cat", "r")    #==================================================================================================HERE
    else:
        #template_file = open("./catalogues/"+cluster_list[l]+"_3p5sigma_edgesourcedelete.cat", "r")    #======================================================================================HERE
        template_file_true = open("../catalogues/edgesourcedelete_deboost_190812/"+cluster_list[kk]+"_3p5sigma_deboost.cat", "r")    #=====================================================================HERE


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
        flux_mjy_true_pre=response_curve_true[:,9]      #=========================================================HERE: 7 for Todd's cluster, 6 for our clusters
        fluxerr_mjy_true_pre=response_curve_true[:,14]    #=========================================================HERE: 8 for Todd's cluster, 7 for our clusters
        sn_ratio_true_pre=response_curve_true[:,5]     #=========================================================HERE: 5 for Todd's cluster, 4 for our clusters
        sn_ratio_confusion_true_pre=response_curve_true[:,12]     #=========================================================HERE: 5 for Todd's cluster, 4 for our clusters

    sn_ratio_true_pre = np.array(sn_ratio_true_pre)
    sn_ratio_confusion_true_pre = np.array(sn_ratio_confusion_true_pre)
    flux_mjy_true_pre = np.array(flux_mjy_true_pre)
    fluxerr_mjy_true_pre = np.array(fluxerr_mjy_true_pre)
    sn_ratio_true_pre=np.array([float(i) for i in sn_ratio_true_pre])  #(wavelength in angstrom)
    sn_ratio_confusion_true_pre=np.array([float(i) for i in sn_ratio_confusion_true_pre])  #(wavelength in angstrom)
    flux_mjy_true_pre=np.array([float(i) for i in flux_mjy_true_pre])  #(wavelength in angstrom)
    fluxerr_mjy_true_pre=np.array([float(i) for i in fluxerr_mjy_true_pre])  #(wavelength in angstrom)

    sn_ratio_true=[]
    sn_ratio_confusion_true=[]
    flux_mjy_true=[]
    fluxerr_mjy_true=[]
    for hd in range(0, len(sn_ratio_true_pre)):
        if abs(sn_ratio_confusion_true_pre[hd]) >= 3.5:
            sn_ratio_true.append(sn_ratio_true_pre[hd])
            sn_ratio_confusion_true.append(sn_ratio_confusion_true_pre[hd])
            flux_mjy_true.append(flux_mjy_true_pre[hd])
            fluxerr_mjy_true.append(fluxerr_mjy_true_pre[hd])


    counts_array_true=np.arange(1.0,len(sn_ratio_true)+1.0)

    counts_array_true = np.array(counts_array_true)
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



