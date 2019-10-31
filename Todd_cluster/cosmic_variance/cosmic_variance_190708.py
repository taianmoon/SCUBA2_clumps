import pyfits
import matplotlib.pyplot as plt
import numpy as np
import math
import random


random_number_amount=10000 #======================================================================HERE

#=====================Do COSMOS field

names_wht = pyfits.open('./COSMOS_cat_myown.fits')
names_wht_data = names_wht[1].data
RA_DEG= names_wht_data.field('RA')   #=====================================================HERE: change column names when confirm whether to use confusion
Dec_DEG= names_wht_data.field('Dec')
detection_SNR= names_wht_data.field('SN_ratio')
S_850_deboost= names_wht_data.field('flux_mJy_perbeam_deboost')
delta_S_850_deboost= names_wht_data.field('flux_error_mJy_deboost')


RA_DEG = np.array(RA_DEG)
RA_DEG=np.array([float(i) for i in RA_DEG])
Dec_DEG = np.array(Dec_DEG)
Dec_DEG=np.array([float(i) for i in Dec_DEG])
detection_SNR = np.array(detection_SNR)
detection_SNR=np.array([float(i) for i in detection_SNR])
S_850_deboost = np.array(S_850_deboost)
S_850_deboost=np.array([float(i) for i in S_850_deboost])
delta_S_850_deboost = np.array(delta_S_850_deboost)
delta_S_850_deboost=np.array([float(i) for i in delta_S_850_deboost])


RA_cat=[]
Dec_cat=[]
SNR_cat=[]
S850_deboost_cat=[]
S850_deboost_err_cat=[]
for i in range(0, len(RA_DEG)):
    if S_850_deboost[i] >= 6.0 and detection_SNR[i] >= 3.5:  #=================================================HERE: change flux density and S/N threshold
        RA_cat.append(RA_DEG[i])
        Dec_cat.append(Dec_DEG[i])
        SNR_cat.append(detection_SNR[i])
        S850_deboost_cat.append(S_850_deboost[i])
        S850_deboost_err_cat.append(delta_S_850_deboost[i])


print len(SNR_cat)
#print SNR_cat


random_centre_RA=[]
random_centre_Dec=[]
for j in range(0, random_number_amount): 
    random_centre_RA.append(random.uniform(149.70814, 150.52018))
    random_centre_Dec.append(random.uniform(1.8223541, 2.6296625))


#print random_centre_RA
#print random_centre_Dec



cosmic_variance_dist=[]
for k in range(0, len(random_centre_RA)):
    counter=0
    for cx in range(0, len(RA_cat)):
        if math.sqrt(pow(random_centre_RA[k]-RA_cat[cx], 2)+pow(random_centre_Dec[k]-Dec_cat[cx], 2)) <= 350.0/3600.0:
            counter=counter+1
    cosmic_variance_dist.append(counter)


print cosmic_variance_dist


print np.mean(cosmic_variance_dist)
print np.std(cosmic_variance_dist)



#=====================Check cosmic variance probability for each field=================================

#print 'cosmic variance probability'
#print 'Bootes1:'
#print np.count_nonzero(np.array(cosmic_variance_dist) >= 6.0)  #=========================change when table is updated after confusion noise is removed
#print np.count_nonzero(np.array(cosmic_variance_dist) >= 6.0)/10000.0  #=========================change when table is updated after confusion noise is removed

#print 'EGS:'
#print np.count_nonzero(np.array(cosmic_variance_dist) >= 14.0)  #=========================change when table is updated after confusion noise is removed
#print np.count_nonzero(np.array(cosmic_variance_dist) >= 14.0)/10000.0  #=========================change when table is updated after confusion noise is removed

#print 'G12:'
#print np.count_nonzero(np.array(cosmic_variance_dist) >= 8.7)  #=========================change when table is updated after confusion noise is removed
#print np.count_nonzero(np.array(cosmic_variance_dist) >= 8.7)/10000.0  #=========================change when table is updated after confusion noise is removed

#print 'Lockman:'
#print np.count_nonzero(np.array(cosmic_variance_dist) >= 1.7)  #=========================change when table is updated after confusion noise is removed
#print np.count_nonzero(np.array(cosmic_variance_dist) >= 1.7)/10000.0  #=========================change when table is updated after confusion noise is removed

#print 'NGP1:'
#print np.count_nonzero(np.array(cosmic_variance_dist) >= 3.8)  #=========================change when table is updated after confusion noise is removed
#print np.count_nonzero(np.array(cosmic_variance_dist) >= 3.8)/10000.0  #=========================change when table is updated after confusion noise is removed

#print 'NGP2:'
#print np.count_nonzero(np.array(cosmic_variance_dist) >= 6.5)  #=========================change when table is updated after confusion noise is removed
#print np.count_nonzero(np.array(cosmic_variance_dist) >= 6.5)/10000.0  #=========================change when table is updated after confusion noise is removed

#print 'NGP3:'
#print np.count_nonzero(np.array(cosmic_variance_dist) >= 5.1)  #=========================change when table is updated after confusion noise is removed
#print np.count_nonzero(np.array(cosmic_variance_dist) >= 5.1)/10000.0  #=========================change when table is updated after confusion noise is removed

#print 'NGP4:'
#print np.count_nonzero(np.array(cosmic_variance_dist) >= 17.0) #=========================change when table is updated after confusion noise is removed
#print np.count_nonzero(np.array(cosmic_variance_dist) >= 17.0)/10000.0  #=========================change when table is updated after confusion noise is removed

#print 'NGP5:'
#print np.count_nonzero(np.array(cosmic_variance_dist) >= 5.7)  #=========================change when table is updated after confusion noise is removed
#print np.count_nonzero(np.array(cosmic_variance_dist) >= 5.7)/10000.0  #=========================change when table is updated after confusion noise is removed

#print 'NGP6:'
#print np.count_nonzero(np.array(cosmic_variance_dist) >= 1.0)  #=========================change when table is updated after confusion noise is removed
#print np.count_nonzero(np.array(cosmic_variance_dist) >= 1.0)/10000.0  #=========================change when table is updated after confusion noise is removed

#print 'NGP7:'
#print np.count_nonzero(np.array(cosmic_variance_dist) >= 13.0)  #=========================change when table is updated after confusion noise is removed
#print np.count_nonzero(np.array(cosmic_variance_dist) >= 13.0)/10000.0  #=========================change when table is updated after confusion noise is removed

#print 'NGP8:'
#print np.count_nonzero(np.array(cosmic_variance_dist) >= 16.0)  #=========================change when table is updated after confusion noise is removed
#print np.count_nonzero(np.array(cosmic_variance_dist) >= 16.0)/10000.0  #=========================change when table is updated after confusion noise is removed

#print 'NGP9:'
#print np.count_nonzero(np.array(cosmic_variance_dist) >= 3.5)  #=========================change when table is updated after confusion noise is removed
#print np.count_nonzero(np.array(cosmic_variance_dist) >= 3.5)/10000.0  #=========================change when table is updated after confusion noise is removed


S_8mJy = [7.1, 1.1, 6.3, 8.4, 9.8, 2.7, 12.0, 13.0, 13.0, 21.0, 9.0, 9.6, 4.7, 9.1, 12.0, 17.0, 5.6, 2.2, 6.9, 29.0, 3.5, 20.0, 4.6, 7.4, 9.4, 8.8, 18.0, 21.0, 9.9, 20.0, 7.9, 16.0, 8.5, 37.0, 3.7, 6.6, 12.0, 9.7, 7.4, 11.0, 7.9, 4.4, 20.0, 7.7, 12.0, 14.0]


print len(S_8mJy)

N_overdensity = []
P_overdensity = []

for i in range(0, len(S_8mJy)):
    N_overdensity.append(np.count_nonzero(np.array(cosmic_variance_dist) >= S_8mJy[i]))
    P_overdensity.append(np.count_nonzero(np.array(cosmic_variance_dist) >= S_8mJy[i])/10000.0)


mock_catalogue=np.column_stack((N_overdensity, P_overdensity))



np.savetxt('./cosmic_variance_result_6mJy.cat', mock_catalogue, delimiter=' ', header='N_overdendity P_overdensity')









