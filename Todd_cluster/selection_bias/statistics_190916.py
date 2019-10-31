import numpy as np
import corner
import matplotlib.pyplot as plt
import pyfits
import math



template_file_z = open("./C19_z_cutz6.cat", "r")    #==================================================================================================HERE
lines_z = template_file_z.readlines()
template_file_z.close()
response_curve_z=[]
for i in range(0, len(lines_z)):
    separated_lines_z=lines_z[i].split() 
    response_curve_z.append(separated_lines_z)
response_curve_z = np.array(response_curve_z)
distribution_z=response_curve_z[:,0]      
distribution_z = np.array(distribution_z)
distribution_z=np.array([float(i) for i in distribution_z])  #(wavelength in angstrom)
#print distribution_z


template_file_L = open("./C19_L_cutz6.cat", "r")    #==================================================================================================HERE
lines_L = template_file_L.readlines()
template_file_L.close()
response_curve_L=[]
for i in range(0, len(lines_L)):
    separated_lines_L=lines_L[i].split() 
    response_curve_L.append(separated_lines_L)
response_curve_L = np.array(response_curve_L)
distribution_L=response_curve_L[:,0]      
distribution_L = np.array(distribution_L)
distribution_L=np.array([float(i) for i in distribution_L])  #(wavelength in angstrom)
#print distribution_L




template_file_z_Todd = open("./z_distribution_Todd.cat", "r")    #==================================================================================================HERE
lines_z_Todd = template_file_z_Todd.readlines()
template_file_z_Todd.close()
response_curve_z_Todd=[]
for i in range(0, len(lines_z_Todd)):
    separated_lines_z_Todd=lines_z_Todd[i].split() 
    response_curve_z_Todd.append(separated_lines_z_Todd)
response_curve_z_Todd = np.array(response_curve_z_Todd)
distribution_z_Todd=response_curve_z_Todd[:,0]      
distribution_z_Todd = np.array(distribution_z_Todd)
distribution_z_Todd=np.array([float(i) for i in distribution_z_Todd])  #(wavelength in angstrom)
#print distribution_z_Todd


template_file_L_Todd = open("./L_distribution_Todd.cat", "r")    #==================================================================================================HERE
lines_L_Todd = template_file_L_Todd.readlines()
template_file_L_Todd.close()
response_curve_L_Todd=[]
for i in range(0, len(lines_L_Todd)):
    separated_lines_L_Todd=lines_L_Todd[i].split() 
    response_curve_L_Todd.append(separated_lines_L_Todd)
response_curve_L_Todd = np.array(response_curve_L_Todd)
distribution_L_Todd=response_curve_L_Todd[:,0]      
distribution_L_Todd = np.array(distribution_L_Todd)
distribution_L_Todd=np.array([float(i) for i in distribution_L_Todd])  #(wavelength in angstrom)
#print distribution_L_Todd



z_C19_mean = np.mean(distribution_z)
z_C19_std = np.std(distribution_z)

z_M17_mean = np.mean(distribution_z_Todd)
z_M17_std = np.std(distribution_z_Todd)

L_C19_mean = np.mean(distribution_L)
L_C19_std =  np.std(distribution_L)

L_M17_mean =  np.mean(distribution_L_Todd)
L_M17_std =  np.std(distribution_L_Todd)



print 'C19_z: mean='+str(z_C19_mean)+', std='+str(z_C19_std)
print 'C19_L: mean='+str(L_C19_mean)+', std='+str(L_C19_std)

print 'M17_z: mean='+str(z_M17_mean)+', std='+str(z_M17_std)
print 'M17_L: mean='+str(L_M17_mean)+', std='+str(L_M17_std)








