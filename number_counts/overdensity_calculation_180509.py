#import matplotlib.pyplot as plt
#import numpy as np
#import math
#import matplotlib.cm as cm
#from astropy.io import fits
#import astropy.wcs as wcs
#from scipy.interpolate import interp1d



#===============================Cumulative number counts=============================================================


#geach_flux_cum=[3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5]  #mJy
#geach_cum_number_cum=[1012.3, 508.0, 271.9, 151.8, 85.3, 47.1, 26.4, 14.5, 8.7, 5.5, 3.2, 2.4, 1.8]  #/deg^2
#geach_cum_number_err_up_cum=[19.6, 12.3, 8.5, 6.2, 4.7, 3.6, 2.8, 2.2, 1.8, 1.5, 1.2, 1.1, 1.0] #Poisson error
#geach_cum_number_err_low_cum=[19.2, 12.0, 8.2, 6.0, 4.4, 3.3, 2.5, 1.9, 1.5, 1.2, 0.9, 0.8, 0.7]


#------------------------------Load in cum number counts for our clusters---------------------------------------------



#template_file_cum = open("./cum_counts_catalogues/cum_counts_"+cluster_list[l]+".cat", "r")    #========================================================HERE
template_file_cum = open("./180524_cum_sensitivity_area_deboost_bin2/cum_counts_"+cluster_list[l]+"_bin2.cat", "r")    #========================================================HERE

lines_cum = template_file_cum.readlines()
template_file_cum.close()

response_curve_cum=[]

for i in range(0, len(lines_cum)):
    separated_lines_cum=lines_cum[i].split() 
    response_curve_cum.append(separated_lines_cum)


response_curve_cum = np.array(response_curve_cum)
flux_cum=response_curve_cum[:,0]      #=========================================================HERE: 7 for Todd's cluster, 6 for our clusters
cum_number=response_curve_cum[:,1]    #=========================================================HERE: 8 for Todd's cluster, 7 for our clusters
cum_number_err=response_curve_cum[:,2]     #=========================================================HERE: 5 for Todd's cluster, 4 for our clusters


flux_cum = np.array(flux_cum)
flux_cum=np.array([float(i) for i in flux_cum])  #(wavelength in angstrom)
cum_number = np.array(cum_number)
cum_number=np.array([float(i) for i in cum_number])  #(wavelength in angstrom)
cum_number_err = np.array(cum_number_err)
cum_number_err=np.array([float(i) for i in cum_number_err])  #(wavelength in angstrom)


#print flux_cum
#print cum_number
#print cum_number_err


#---------------------------Calculate overdensities in cum number counts-----------------------------------------------------

overdensity_cum=[]
flux_final_cum=[]
n_cum_final=[]
n_cum_err_final=[]
ratio_cum=[]

for i in range(0, len(flux_cum)):

    for j in range(0, len(geach_flux_cum)):

        if flux_cum[i] == geach_flux_cum[j]:

            overdensity_cum.append((cum_number[i]-geach_cum_number_cum[j])/geach_cum_number_err_up_cum[j])
            flux_final_cum.append(flux_cum[i])
            n_cum_final.append(cum_number[i])
            n_cum_err_final.append(cum_number_err[i])
            ratio_cum.append(cum_number[i]/geach_cum_number_cum[j])



overdensity_cum = np.array(overdensity_cum)
overdensity_cum=np.array([float(i) for i in overdensity_cum])  #(wavelength in angstrom)
flux_final_cum = np.array(flux_final_cum)
flux_final_cum=np.array([float(i) for i in flux_final_cum])  #(wavelength in angstrom)
n_cum_final = np.array(n_cum_final)
n_cum_final=np.array([float(i) for i in n_cum_final])  #(wavelength in angstrom)
n_cum_err_final = np.array(n_cum_err_final)
n_cum_err_final=np.array([float(i) for i in n_cum_err_final])  #(wavelength in angstrom)
ratio_cum = np.array(ratio_cum)
ratio_cum=np.array([float(i) for i in ratio_cum])  #(wavelength in angstrom)


print flux_final_cum
print overdensity_cum
print ratio_cum


mock_catalogue_cum=np.column_stack((flux_final_cum, n_cum_final, n_cum_err_final, overdensity_cum, ratio_cum))

#np.savetxt('./cum_counts_'+cluster_list[l]+'_overdensity.cat', mock_catalogue_cum, delimiter=' ') #================================================================================HERE
np.savetxt('./counts_overdensity_catalogues_180607/cum_counts_'+cluster_list[l]+'_bin2_overdensity.cat', mock_catalogue_cum, delimiter=' ') #===================================================HERE


#plt.plot(flux_final_cum, overdensity_cum , color='blue' , linewidth=1.5)


#===============================Diff number counts=============================================================


#geach_flux_diff=[3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5]  #mJy
#geach_diff_number_diff=[451.0, 204.4, 102.6, 56.1, 32.5, 18.0, 9.8, 5.8, 3.4, 2.1, 0.8, 0.5, 0.3]  #/deg^2 /mJy
#geach_diff_number_err_up_diff=[17.1, 9.3, 6.0, 4.3, 3.2, 2.5, 1.9, 1.5, 1.2, 1.1, 0.8, 0.7, 0.6] #Poisson error
#geach_diff_number_err_low_diff=[16.4, 8.9, 5.7, 4.0, 2.9, 2.2, 1.6, 1.2, 0.9, 0.7, 0.4, 0.3, 0.2]


#------------------------------Load in diff number counts for our clusters---------------------------------------------




#template_file_diff = open("./diff_counts_catalogues_180416/diff_counts_"+cluster_list[l]+".cat", "r")    #========================================================HERE
template_file_diff = open("./180524_diff_sensitivity_area_deboost_bin2/diff_counts_"+cluster_list[l]+"_bin2.cat", "r")    #========================================================HERE

lines_diff = template_file_diff.readlines()
template_file_diff.close()

response_curve_diff=[]

for i in range(0, len(lines_diff)):
    separated_lines_diff=lines_diff[i].split() 
    response_curve_diff.append(separated_lines_diff)


response_curve_diff = np.array(response_curve_diff)
flux_diff=response_curve_diff[:,0]      #=========================================================HERE: 7 for Todd's cluster, 6 for our clusters
diff_number=response_curve_diff[:,1]    #=========================================================HERE: 8 for Todd's cluster, 7 for our clusters
diff_number_err=response_curve_diff[:,2]     #=========================================================HERE: 5 for Todd's cluster, 4 for our clusters


flux_diff = np.array(flux_diff)
flux_diff=np.array([float(i) for i in flux_diff])  #(wavelength in angstrom)
diff_number = np.array(diff_number)
diff_number=np.array([float(i) for i in diff_number])  #(wavelength in angstrom)
diff_number_err = np.array(diff_number_err)
diff_number_err=np.array([float(i) for i in diff_number_err])  #(wavelength in angstrom)




#---------------------------Calculate overdensities in diff number counts-----------------------------------------------------

overdensity_diff=[]
flux_final_diff=[]
n_diff_final=[]
n_diff_err_final=[]
ratio_diff=[]

for i in range(0, len(flux_diff)):

    for j in range(0, len(geach_flux_diff)):

        if flux_diff[i] == geach_flux_diff[j]:

            overdensity_diff.append((diff_number[i]-geach_diff_number_diff[j])/geach_diff_number_err_up_diff[j])
            flux_final_diff.append(flux_diff[i])
            n_diff_final.append(diff_number[i])
            n_diff_err_final.append(diff_number_err[i])
            ratio_diff.append(diff_number[i]/geach_diff_number_diff[j])

overdensity_diff = np.array(overdensity_diff)
overdensity_diff=np.array([float(i) for i in overdensity_diff])  #(wavelength in angstrom)
flux_final_diff = np.array(flux_final_diff)
flux_final_diff=np.array([float(i) for i in flux_final_diff])  #(wavelength in angstrom)
n_diff_final = np.array(n_diff_final)
n_diff_final=np.array([float(i) for i in n_diff_final])  #(wavelength in angstrom)
n_diff_err_final = np.array(n_diff_err_final)
n_diff_err_final=np.array([float(i) for i in n_diff_err_final])  #(wavelength in angstrom)
ratio_diff = np.array(ratio_diff)
ratio_diff=np.array([float(i) for i in ratio_diff])  #(wavelength in angstrom)


print flux_final_diff
print overdensity_diff
print ratio_diff


mock_catalogue_diff=np.column_stack((flux_final_diff, n_diff_final, n_diff_err_final, overdensity_diff, ratio_diff))


#np.savetxt('./diff_counts_'+cluster_list[l]+'_overdensity.cat', mock_catalogue_diff, delimiter=' ') #================================================================================HERE
np.savetxt('./counts_overdensity_catalogues_180607/diff_counts_'+cluster_list[l]+'_bin2_overdensity.cat', mock_catalogue_diff, delimiter=' ') #==========================================================HERE

#plt.plot(flux_final_diff, overdensity_diff , color='blue' , linewidth=1.5)







