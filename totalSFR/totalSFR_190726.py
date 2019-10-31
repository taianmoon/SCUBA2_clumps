import numpy as np
import corner
import pyfits
import math
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord


#note: G12 and NGP1 use de-magnified flux for the lens!!


names_wht = pyfits.open('./source_cat_190722.fits')    #===============================HERE
names_wht_data = names_wht[1].data
RA_pre= names_wht_data.field('RA')
Dec_pre= names_wht_data.field('Dec')
S850_pre= names_wht_data.field('S850')
S250_pre= names_wht_data.field('S250')
S350_pre= names_wht_data.field('S350')
S500_pre= names_wht_data.field('S500')
z_pre= names_wht_data.field('z')
z_uperr_pre= names_wht_data.field('z_uperr')
z_lowerr_pre= names_wht_data.field('z_lowerr')
L_pre= names_wht_data.field('L')
L_uperr_pre= names_wht_data.field('L_uperr')
L_lowerr_pre= names_wht_data.field('L_lowerr')
SFR_pre= names_wht_data.field('SFR')
SFR_uperr_pre= names_wht_data.field('SFR_uperr')
SFR_lowerr_pre= names_wht_data.field('SFR_lowerr')


cluster_list=['Bootes1', 'EGS', 'G12', 'Lockman', 'NGP1', 'NGP2', 'NGP3', 'NGP4', 'NGP5', 'NGP6', 'NGP7', 'NGP8', 'NGP9']  #must be in this sequence
#cluster_list=['Bootes1', 'EGS']
clusrer_list_sourceindexnumber_start=[0, 12, 16, 31, 33, 40, 45, 48, 57, 64, 72, 85, 90]
clusrer_list_sourceindexnumber_end=[11, 15, 30, 32, 39, 44, 47, 56, 63, 71, 84, 89, 93]


no_members=[]
SFR_total=[]
SFR_total_uperr=[]
SFR_total_lowerr=[]
L_total=[]
L_total_uperr=[]
L_total_lowerr=[]
size_arcmin=[]
cluster_name=[]
z_mean=[]
z_mean_err=[]
for l in range(0, len(cluster_list)):    #=================================Loop over clusters
    print cluster_list[l]
    cluster_name.append(cluster_list[l])

    totalSFR=0.0
    totalSFR_uperr=0.0
    totalSFR_lowerr=0.0
    totalL=0.0
    totalL_uperr=0.0
    totalL_lowerr=0.0
    counter=0
    z_mean_nominator=0.0
    z_mean_denominator=0.0
    z_mean_err_loop=0.0
    RA_members=[]
    Dec_members=[]
    for i in range(clusrer_list_sourceindexnumber_start[l], clusrer_list_sourceindexnumber_end[l]+1):   #=================================Loop over sources

        if cluster_list[l]=='NGP1' or cluster_list[l]=='NGP6':
            #if  z_pre[i] <= 6.0 and S250_pre[i]/S350_pre[i] < 1.7 and S350_pre[i]/S850_pre[i] < 10.0 and not (S250_pre[i]/S350_pre[i] <= 1.3 and S350_pre[i]/S850_pre[i] <= 6.0): 
            #if S250_pre[i]/S350_pre[i] < 1.7 and S350_pre[i]/S850_pre[i] < 10.0 and not (S250_pre[i]/S350_pre[i] <= 1.3 and S350_pre[i]/S850_pre[i] <= 6.0): 
            if  z_pre[i] <= 6.0:

                counter=counter+1
        
                #print pow(10, SFR_pre[i])
                totalSFR=totalSFR+pow(10, SFR_pre[i])
                #print pow(abs(pow(10, SFR_pre[i]+SFR_uperr_pre[i])-pow(10, SFR_pre[i])), 2)
                #totalSFR_uperr=totalSFR_uperr+pow(abs(pow(10, SFR_pre[i]+SFR_uperr_pre[i])-pow(10, SFR_pre[i])), 2)
                totalSFR_uperr=totalSFR_uperr+pow(2.303*pow(10, SFR_pre[i])*SFR_uperr_pre[i], 2)
                #print pow(abs(pow(10, SFR_pre[i]-SFR_lowerr_pre[i])-pow(10, SFR_pre[i])), 2)
                #totalSFR_lowerr=totalSFR_lowerr+pow(abs(pow(10, SFR_pre[i]-SFR_lowerr_pre[i])-pow(10, SFR_pre[i])), 2)
                totalSFR_lowerr=totalSFR_lowerr+pow(2.303*pow(10, SFR_pre[i])*SFR_lowerr_pre[i], 2)
        
                #print pow(10, L_pre[i])
                totalL=totalL+pow(10, L_pre[i])
                #print pow(abs(pow(10, L_pre[i]+L_uperr_pre[i])-pow(10, L_pre[i])), 2)
                #totalL_uperr=totalL_uperr+pow(abs(pow(10, L_pre[i]+L_uperr_pre[i])-pow(10, L_pre[i])), 2)
                totalL_uperr=totalL_uperr+pow(2.303*pow(10, L_pre[i])*L_uperr_pre[i], 2)
                #print pow(abs(pow(10, L_pre[i]-L_lowerr_pre[i])-pow(10, L_pre[i])), 2)
                #totalL_lowerr=totalL_lowerr+pow(abs(pow(10, L_pre[i]-L_lowerr_pre[i])-pow(10, L_pre[i])), 2)
                totalL_lowerr=totalL_lowerr+pow(2.303*pow(10, L_pre[i])*L_lowerr_pre[i], 2)

                z_mean_nominator=z_mean_nominator+(z_pre[i]/((z_uperr_pre[i]+z_uperr_pre[i])/2.0))
                z_mean_denominator=z_mean_denominator+(1.0/((z_uperr_pre[i]+z_uperr_pre[i])/2.0))  
                z_mean_err_loop=z_mean_err_loop+(1.0/((z_uperr_pre[i]+z_uperr_pre[i])/2.0))              

                RA_members.append(RA_pre[i])
                Dec_members.append(Dec_pre[i])

        else:
            #if  z_pre[i] <= 6.0 and S250_pre[i]/S350_pre[i] <= 1.3 and S350_pre[i]/S850_pre[i] <= 6.0:
            #if S250_pre[i]/S350_pre[i] <= 1.3 and S350_pre[i]/S850_pre[i] <= 6.0:
            if  z_pre[i] <= 6.0:

                counter=counter+1
        
                #print pow(10, SFR_pre[i])
                totalSFR=totalSFR+pow(10, SFR_pre[i])
                #print pow(abs(pow(10, SFR_pre[i]+SFR_uperr_pre[i])-pow(10, SFR_pre[i])), 2)
                #totalSFR_uperr=totalSFR_uperr+pow(abs(pow(10, SFR_pre[i]+SFR_uperr_pre[i])-pow(10, SFR_pre[i])), 2)
                totalSFR_uperr=totalSFR_uperr+pow(2.303*pow(10, SFR_pre[i])*SFR_uperr_pre[i], 2)
                #print pow(abs(pow(10, SFR_pre[i]-SFR_lowerr_pre[i])-pow(10, SFR_pre[i])), 2)
                #totalSFR_lowerr=totalSFR_lowerr+pow(abs(pow(10, SFR_pre[i]-SFR_lowerr_pre[i])-pow(10, SFR_pre[i])), 2)
                totalSFR_lowerr=totalSFR_lowerr+pow(2.303*pow(10, SFR_pre[i])*SFR_lowerr_pre[i], 2)

        
                #print pow(10, L_pre[i])
                totalL=totalL+pow(10, L_pre[i])
                #print pow(abs(pow(10, L_pre[i]+L_uperr_pre[i])-pow(10, L_pre[i])), 2)
                #totalL_uperr=totalL_uperr+pow(abs(pow(10, L_pre[i]+L_uperr_pre[i])-pow(10, L_pre[i])), 2)
                totalL_uperr=totalL_uperr+pow(2.303*pow(10, L_pre[i])*L_uperr_pre[i], 2)
                #print pow(abs(pow(10, L_pre[i]-L_lowerr_pre[i])-pow(10, L_pre[i])), 2)
                #totalL_lowerr=totalL_lowerr+pow(abs(pow(10, L_pre[i]-L_lowerr_pre[i])-pow(10, L_pre[i])), 2)
                totalL_lowerr=totalL_lowerr+pow(2.303*pow(10, L_pre[i])*L_lowerr_pre[i], 2)

                z_mean_nominator=z_mean_nominator+(z_pre[i]/((z_uperr_pre[i]+z_uperr_pre[i])/2.0))
                z_mean_denominator=z_mean_denominator+(1.0/((z_uperr_pre[i]+z_uperr_pre[i])/2.0))  
                z_mean_err_loop=z_mean_err_loop+(1.0/((z_uperr_pre[i]+z_uperr_pre[i])/2.0))                    

                RA_members.append(RA_pre[i])
                Dec_members.append(Dec_pre[i])


    print counter
    print totalSFR
    print math.sqrt(totalSFR_uperr)
    print math.sqrt(totalSFR_lowerr)
    print np.array(totalL/1e13)
    print np.array(math.sqrt(totalL_uperr)/1e13)
    print np.array(math.sqrt(totalL_lowerr)/1e13)
    #print RA_members
    #print Dec_members
    no_members.append(counter)
    SFR_total.append(totalSFR)
    SFR_total_uperr.append(math.sqrt(totalSFR_uperr))
    SFR_total_lowerr.append(math.sqrt(totalSFR_lowerr))
    L_total.append(np.array(totalL/1e13))
    L_total_uperr.append(np.array(math.sqrt(totalL_uperr)/1e13))
    L_total_lowerr.append(np.array(math.sqrt(totalL_lowerr)/1e13))

    z_mean.append(z_mean_nominator/z_mean_denominator)
    z_mean_err.append(1.0/(math.sqrt(z_mean_err_loop)))

    size_loop=[]
    for pi in range(0, len(RA_members)):
        for pj in range(0, len(Dec_members)):
            #size_loop.append(math.sqrt(pow(RA_members[pi]-RA_members[pj], 2)+pow(Dec_members[pi]-Dec_members[pj], 2)))
            c1 = SkyCoord(RA_members[pi], Dec_members[pi], unit='deg')
            c2 = SkyCoord(RA_members[pj], Dec_members[pj], unit='deg')
            size_loop.append(c1.separation(c2).arcminute)

    size_loop=np.array(size_loop)
    #print size_loop
    print max(size_loop)  #in arcmin
    size_arcmin.append(max(size_loop))

print 'sdddddddddddddddddddddddddddddddd'
print z_mean
print z_mean_err


mock_catalogue_z=np.column_stack((no_members, z_mean, z_mean_err, SFR_total, SFR_total_uperr, SFR_total_lowerr, L_total, L_total_uperr, L_total_lowerr, size_arcmin))
np.savetxt('./cluster_totalSFR.cat', mock_catalogue_z, delimiter=' ') #================================================================================HERE






