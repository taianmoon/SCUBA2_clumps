import numpy as np
import corner
import matplotlib.pyplot as plt
import pyfits
import math
 






template_file = open("./marginalized/allsource_z_S2CLS_random.cat", "r")    #==================================================================================================HERE
lines = template_file.readlines()
template_file.close()
response_curve=[]
for i in range(0, len(lines)):
    separated_lines=lines[i].split() 
    response_curve.append(separated_lines)
response_curve = np.array(response_curve)
distribution_all=response_curve[:,0]
distribution_all_errup=response_curve[:,1]      
distribution_all_errlow=response_curve[:,2]            
distribution_all = np.array(distribution_all)
distribution_all=np.array([float(i) for i in distribution_all])  #(wavelength in angstrom)
#print distribution_all

template_file_SFR = open("./marginalized/allsource_SFR_S2CLS_random.cat", "r")    #==================================================================================================HERE
lines_SFR = template_file_SFR.readlines()
template_file_SFR.close()
response_curve_SFR=[]
for i in range(0, len(lines_SFR)):
    separated_lines_SFR=lines_SFR[i].split() 
    response_curve_SFR.append(separated_lines_SFR)
response_curve_SFR = np.array(response_curve_SFR)
distribution_all_SFR=response_curve_SFR[:,0]
distribution_all_SFR_errup=response_curve_SFR[:,1]      
distribution_all_SFR_errlow=response_curve_SFR[:,2]                  
distribution_all_SFR = np.array(distribution_all_SFR)
distribution_all_SFR=np.array([float(i) for i in distribution_all_SFR])  #(wavelength in angstrom)
#print distribution_all_SFR

#template_file_a = open("./marginalized/allsource_a.cat", "r")    #==================================================================================================HERE
#lines_a = template_file_a.readlines()
#template_file_a.close()
#response_curve_a=[]
#for i in range(0, len(lines_a)):
#    separated_lines_a=lines_a[i].split() 
#    response_curve_a.append(separated_lines_a)
#response_curve_a = np.array(response_curve_a)
#distribution_all_a=response_curve_a[:,0]      
#distribution_all_a_errup=response_curve_a[:,1]      
#distribution_all_a_errlow=response_curve_a[:,2]                  
#distribution_all_a = np.array(distribution_all_a)
#distribution_all_a=np.array([float(i) for i in distribution_all_a])  #(wavelength in angstrom)
##print distribution_all_a

template_file_L = open("./marginalized/allsource_L_S2CLS_random.cat", "r")    #==================================================================================================HERE
lines_L = template_file_L.readlines()
template_file_L.close()
response_curve_L=[]
for i in range(0, len(lines_L)):
    separated_lines_L=lines_L[i].split() 
    response_curve_L.append(separated_lines_L)
response_curve_L = np.array(response_curve_L)
distribution_all_L=response_curve_L[:,0]      
distribution_all_L_errup=response_curve_L[:,1]      
distribution_all_L_errlow=response_curve_L[:,2]                  
distribution_all_L = np.array(distribution_all_L)
distribution_all_L=np.array([float(i) for i in distribution_all_L])  #(wavelength in angstrom)
#print distribution_all_L


distribution_all_z_cut=[]
distribution_all_z_errup_cut=[]
distribution_all_z_errlow_cut=[]
#distribution_all_a_cut=[]
#distribution_all_a_errup_cut=[]
#distribution_all_a_errlow_cut=[]
distribution_all_L_cut=[]
distribution_all_L_errup_cut=[]
distribution_all_L_errlow_cut=[]
distribution_all_SFR_cut=[]
distribution_all_SFR_errup_cut=[]
distribution_all_SFR_errlow_cut=[]

for i in range(0, len(distribution_all)):

    if distribution_all[i] <= 6.0:

        distribution_all_z_cut.append(distribution_all[i])
        distribution_all_z_errup_cut.append(distribution_all_errup[i])
        distribution_all_z_errlow_cut.append(distribution_all_errlow[i])
        #distribution_all_a_cut.append(distribution_all_a[i])
        #distribution_all_a_errup_cut.append(distribution_all_a_errup[i])
        #distribution_all_a_errlow_cut.append(distribution_all_a_errlow[i])
        distribution_all_L_cut.append(distribution_all_L[i])
        distribution_all_L_errup_cut.append(distribution_all_L_errup[i])
        distribution_all_L_errlow_cut.append(distribution_all_L_errlow[i])
        distribution_all_SFR_cut.append(distribution_all_SFR[i])
        distribution_all_SFR_errup_cut.append(distribution_all_SFR_errup[i])
        distribution_all_SFR_errlow_cut.append(distribution_all_SFR_errlow[i])

distribution_all_SFR_cut = np.array(distribution_all_SFR_cut)
distribution_all_SFR_errup_cut = np.array(distribution_all_SFR_errup_cut)
distribution_all_SFR_errlow_cut = np.array(distribution_all_SFR_errlow_cut)
distribution_all_z_cut = np.array(distribution_all_z_cut)
distribution_all_z_errup_cut = np.array(distribution_all_z_errup_cut)
distribution_all_z_errlow_cut = np.array(distribution_all_z_errlow_cut)
#distribution_all_a_cut = np.array(distribution_all_a_cut)
#distribution_all_a_errup_cut = np.array(distribution_all_a_errup_cut)
#distribution_all_a_errlow_cut = np.array(distribution_all_a_errlow_cut)
distribution_all_L_cut = np.array(distribution_all_L_cut)
distribution_all_L_errup_cut = np.array(distribution_all_L_errup_cut)
distribution_all_L_errlow_cut = np.array(distribution_all_L_errlow_cut)

distribution_all_SFR_cut=np.array([float(i) for i in distribution_all_SFR_cut])  #(wavelength in angstrom)
distribution_all_SFR_errup_cut=np.array([float(i) for i in distribution_all_SFR_errup_cut])  #(wavelength in angstrom)
distribution_all_SFR_errlow_cut=np.array([float(i) for i in distribution_all_SFR_errlow_cut])  #(wavelength in angstrom)
distribution_all_z_cut=np.array([float(i) for i in distribution_all_z_cut])  #(wavelength in angstrom)
distribution_all_z_errup_cut=np.array([float(i) for i in distribution_all_z_errup_cut])  #(wavelength in angstrom)
distribution_all_z_errlow_cut=np.array([float(i) for i in distribution_all_z_errlow_cut])  #(wavelength in angstrom)
#distribution_all_a_cut=np.array([float(i) for i in distribution_all_a_cut])  #(wavelength in angstrom)
#distribution_all_a_errup_cut=np.array([float(i) for i in distribution_all_a_errup_cut])  #(wavelength in angstrom)
#distribution_all_a_errlow_cut=np.array([float(i) for i in distribution_all_a_errlow_cut])  #(wavelength in angstrom)
distribution_all_L_cut=np.array([float(i) for i in distribution_all_L_cut])  #(wavelength in angstrom)
distribution_all_L_errup_cut=np.array([float(i) for i in distribution_all_L_errup_cut])  #(wavelength in angstrom)
distribution_all_L_errlow_cut=np.array([float(i) for i in distribution_all_L_errlow_cut])  #(wavelength in angstrom)



#print distribution_all_z_cut
#print distribution_all_a_cut
#print distribution_all_L_cut
#print distribution_all_SFR_cut

print len(distribution_all_z_cut)
print len(distribution_all_z_errup_cut)
print len(distribution_all_z_errlow_cut)

mock_catalogue_z=np.column_stack((distribution_all_z_cut, distribution_all_z_errup_cut, distribution_all_z_errlow_cut))
np.savetxt('./marginalized/allsource_z_S2CLS_random_cutz6.cat', mock_catalogue_z, delimiter=' ') #================================================================================HERE

mock_catalogue_SFR=np.column_stack((distribution_all_SFR_cut, distribution_all_SFR_errup_cut, distribution_all_SFR_errlow_cut))
np.savetxt('./marginalized/allsource_SFR_S2CLS_random_cutz6.cat', mock_catalogue_SFR, delimiter=' ') #================================================================================HERE

#mock_catalogue_a=np.column_stack((distribution_all_a_cut, distribution_all_a_errup_cut, distribution_all_a_errlow_cut))
#np.savetxt('./marginalized/allsource_a_cutz6.cat', mock_catalogue_a, delimiter=' ') #================================================================================HERE

mock_catalogue_L=np.column_stack((distribution_all_L_cut, distribution_all_L_errup_cut, distribution_all_L_errlow_cut))
np.savetxt('./marginalized/allsource_L_S2CLS_random_cutz6.cat', mock_catalogue_L, delimiter=' ') #================================================================================HERE



