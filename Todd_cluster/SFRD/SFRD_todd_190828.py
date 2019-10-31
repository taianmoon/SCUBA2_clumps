import numpy as np
#import corner
import matplotlib.pyplot as plt
#import pyfits
import math


template_file_Todd = open("./z_distribution_Todd.cat", "r")    #==================================================================================================HERE
lines_Todd = template_file_Todd.readlines()
template_file_Todd.close()
response_curve_Todd=[]
for i in range(0, len(lines_Todd)):
    separated_lines_Todd=lines_Todd[i].split() 
    response_curve_Todd.append(separated_lines_Todd)
response_curve_Todd = np.array(response_curve_Todd)
z_all_Todd=response_curve_Todd[:,0]      
z_all_Todd = np.array(z_all_Todd)
z_all_Todd=np.array([float(i) for i in z_all_Todd])  #(wavelength in angstrom)
print len(z_all_Todd)


L_file_Todd = open("./L_distribution_Todd.cat", "r")    #==================================================================================================HERE
L_lines_Todd = L_file_Todd.readlines()
L_file_Todd.close()
L_curve_Todd=[]
for i in range(0, len(L_lines_Todd)):
    separatedL_lines_Todd=L_lines_Todd[i].split() 
    L_curve_Todd.append(separatedL_lines_Todd)
L_curve_Todd = np.array(L_curve_Todd)
L_all_Todd=L_curve_Todd[:,0]      
L_all_Todd = np.array(L_all_Todd)
L_all_Todd=np.array([float(i) for i in L_all_Todd])  #(wavelength in angstrom)
print len(L_all_Todd)
