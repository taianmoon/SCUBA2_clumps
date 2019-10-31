import numpy as np
import corner
import matplotlib.pyplot as plt
import pyfits
import math


#cluster_list=['PCL1002']
cluster_list=['S2CLS_random']
#cluster_list=['Bootes1', 'EGS', 'G12', 'Lockman', 'NGP1', 'NGP2', 'NGP3', 'NGP4', 'NGP5', 'NGP6', 'NGP7', 'NGP8', 'NGP9']
#cluster_list=[ 'Bootes1', 'EGS','NGP4', 'NGP5',  'NGP7', 'NGP8']

#source_name=['Bootes1_01', 'Bootes1_02', 'Bootes1_03', 'Bootes1_04']

cat1_switch = False  #===============================HERE
plot_corner_switch = True                                                                                                                    #===============================HERE

z_distribution_all=[]
z_distribution_all_uperr=[]
z_distribution_all_lowerr=[]

a_distribution_all=[]
a_distribution_all_uperr=[]
a_distribution_all_lowerr=[]

L_distribution_all=[]
L_distribution_all_uperr=[]
L_distribution_all_lowerr=[]

SFR_distribution_all=[]
SFR_distribution_all_uperr=[]
SFR_distribution_all_lowerr=[]

for l in range(0, len(cluster_list)):    #=================================Loop over clusters


    print cluster_list[l]

    if cluster_list[0]=='PCL1002' or cluster_list[0]=='S2CLS_random':
        names_wht = pyfits.open('../'+cluster_list[l]+'/'+cluster_list[l]+'_herschel_match_3p5sigma_cat_edgesourcedelete_allscuba2source_herschelflux.fits')    #===============================HERE
    else:
        names_wht = pyfits.open('../../'+cluster_list[l]+'/'+cluster_list[l]+'_herschel_match_3p5sigma_cat_edgesourcedelete_allscuba2source_herschelflux.fits')    #===============================HERE
    names_wht_data = names_wht[1].data
    flux_250_observed_pre= names_wht_data.field('F250')
    source_name_pre= names_wht_data.field('name')
    SN_ratio_pre= names_wht_data.field('SN_ratio_1')
    SN_ratio_confusion_pre= names_wht_data.field('SN_ratio_confusion')

    flux_250_observed_pre = np.array(flux_250_observed_pre)
    flux_250_observed_pre=np.array([float(i) for i in flux_250_observed_pre])
    SN_ratio_pre = np.array(SN_ratio_pre)
    SN_ratio_pre=np.array([float(i) for i in SN_ratio_pre])
    SN_ratio_confusion_pre = np.array(SN_ratio_confusion_pre)
    SN_ratio_confusion_pre=np.array([float(i) for i in SN_ratio_confusion_pre])

    #flux_250_observed=[]
    source_name=[]

    for i in range(0, len(flux_250_observed_pre)):     

        if  math.isnan(flux_250_observed_pre[i]) == False and SN_ratio_confusion_pre[i]>=3.5:

            #flux_250_observed.append(flux_250_observed_pre[i])
            source_name.append(source_name_pre[i])
    #print source_name

    for x in range(0, len(source_name)):      #=================================Loop over sources

        print source_name[x]
        #print 'samples_values_all_'+str(cluster_list[l])+'_'+str(source_name[x])+'.npy'
        marginalized_array = np.load('./marginalized/samples_values_all_'+str(cluster_list[l])+'_'+str(source_name[x])+'.npy')

        print np.shape(marginalized_array)

        if plot_corner_switch==True:
            fig = corner.corner(marginalized_array[:,0:3], labels=["$z$", "$log(a)$", r"$log(L(L_{\odot}))$"],
                                quantiles=[0.16, 0.5, 0.84],
                                show_titles=True, title_kwargs={"fontsize": 12})
      
            #plt.show()

            fig.savefig("./marginalized/triangle_marginalized_"+str(cluster_list[l])+"_"+str(source_name[x])+".png")  
            plt.close(fig)
    



        lower_limit_z, fitted_median_z, higher_limit_z = corner.quantile(marginalized_array[:,0], q=[0.16, 0.5, 0.84], weights=None)
        z_distribution_all.append(fitted_median_z)
        z_distribution_all_uperr.append(higher_limit_z-fitted_median_z)
        z_distribution_all_lowerr.append(fitted_median_z-lower_limit_z)

        #print '2222222222222222222222222222'
        #print source_name[x]
        #print fitted_median_z
        #print lower_limit_z
        #print higher_limit_z


        lower_limit_a, fitted_median_a, higher_limit_a = corner.quantile(marginalized_array[:,1], q=[0.16, 0.5, 0.84], weights=None)

        a_distribution_all.append(fitted_median_a)
        a_distribution_all_uperr.append(higher_limit_a-fitted_median_a)
        a_distribution_all_lowerr.append(fitted_median_a-lower_limit_a)

        #print fitted_median_a
        #print lower_limit_a
        #print higher_limit_a


        #==================================II. Calculate Luminosity and SFR=================================================================
        #marginalized_array_random=[]

        #np.save('./marginalized/samples_values_all_'+str(cluster_list[l])+'_'+str(source_name[x])+'.npy', marginalized_array)
        #execfile("./luminosity_etc_181113.py")

        lower_limit_L, fitted_median_L, higher_limit_L = corner.quantile(marginalized_array[:,2], q=[0.16, 0.5, 0.84], weights=None)
        L_distribution_all.append(fitted_median_L)
        L_distribution_all_uperr.append(higher_limit_L-fitted_median_L)
        L_distribution_all_lowerr.append(fitted_median_L-lower_limit_L)
        #print fitted_median_L
        #print lower_limit_L
        #print higher_limit_L

        lower_limit_SFR, fitted_median_SFR, higher_limit_SFR = corner.quantile(marginalized_array[:,3], q=[0.16, 0.5, 0.84], weights=None)
        SFR_distribution_all.append(fitted_median_SFR)
        SFR_distribution_all_uperr.append(higher_limit_SFR-fitted_median_SFR)
        SFR_distribution_all_lowerr.append(fitted_median_SFR-lower_limit_SFR)
        #print fitted_median_SFR
        #print lower_limit_SFR
        #print higher_limit_SFR


if cluster_list[0]=='PCL1002' or cluster_list[0]=='S2CLS_random':

    mock_catalogue_z=np.column_stack((z_distribution_all, z_distribution_all_uperr, z_distribution_all_lowerr))
    np.savetxt('./marginalized/allsource_z_'+cluster_list[0]+'.cat', mock_catalogue_z, delimiter=' ') #================================================================================HERE

    mock_catalogue_a=np.column_stack((a_distribution_all, a_distribution_all_uperr, a_distribution_all_lowerr))
    np.savetxt('./marginalized/allsource_a_'+cluster_list[0]+'.cat', mock_catalogue_a, delimiter=' ') #========================================================================HERE:actually it's log(a)


    mock_catalogue_L=np.column_stack((L_distribution_all, L_distribution_all_uperr, L_distribution_all_lowerr))
    np.savetxt('./marginalized/allsource_L_'+cluster_list[0]+'.cat', mock_catalogue_L, delimiter=' ') #================================================================================HERE

    mock_catalogue_SFR=np.column_stack((SFR_distribution_all, SFR_distribution_all_uperr, SFR_distribution_all_lowerr))
    np.savetxt('./marginalized/allsource_SFR_'+cluster_list[0]+'.cat', mock_catalogue_SFR, delimiter=' ') #================================================================================HERE


else:

    if cat1_switch==True:
        mock_catalogue_z=np.column_stack((z_distribution_all, z_distribution_all_uperr, z_distribution_all_lowerr))
        np.savetxt('./marginalized/allsource_z_cat1.cat', mock_catalogue_z, delimiter=' ') #================================================================================HERE

        mock_catalogue_a=np.column_stack((a_distribution_all, a_distribution_all_uperr, a_distribution_all_lowerr))
        np.savetxt('./marginalized/allsource_a_cat1.cat', mock_catalogue_a, delimiter=' ') #================================================================================HERE:actually it's log(a)


        mock_catalogue_L=np.column_stack((L_distribution_all, L_distribution_all_uperr, L_distribution_all_lowerr))
        np.savetxt('./marginalized/allsource_L_cat1.cat', mock_catalogue_L, delimiter=' ') #================================================================================HERE

        mock_catalogue_SFR=np.column_stack((SFR_distribution_all, SFR_distribution_all_uperr, SFR_distribution_all_lowerr))
        np.savetxt('./marginalized/allsource_SFR_cat1.cat', mock_catalogue_SFR, delimiter=' ') #================================================================================HERE


    else:
        mock_catalogue_z=np.column_stack((z_distribution_all, z_distribution_all_uperr, z_distribution_all_lowerr))
        np.savetxt('./marginalized/allsource_z.cat', mock_catalogue_z, delimiter=' ') #================================================================================HERE

        mock_catalogue_a=np.column_stack((a_distribution_all, a_distribution_all_uperr, a_distribution_all_lowerr))
        np.savetxt('./marginalized/allsource_a.cat', mock_catalogue_a, delimiter=' ') #================================================================================HERE:actually it's log(a)


        mock_catalogue_L=np.column_stack((L_distribution_all, L_distribution_all_uperr, L_distribution_all_lowerr))
        np.savetxt('./marginalized/allsource_L.cat', mock_catalogue_L, delimiter=' ') #================================================================================HERE

        mock_catalogue_SFR=np.column_stack((SFR_distribution_all, SFR_distribution_all_uperr, SFR_distribution_all_lowerr))
        np.savetxt('./marginalized/allsource_SFR.cat', mock_catalogue_SFR, delimiter=' ') #================================================================================HERE







