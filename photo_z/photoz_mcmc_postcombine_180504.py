import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm
from astropy.io import fits
import astropy.wcs as wcs
from scipy.interpolate import interp1d
import pyfits


#======================================extract source names===============================================================================================================

names_wht = pyfits.open('../Bootes1/Bootes1_herschel_match_3p5sigma_cat_edgesourcedelete_allscuba2source.fits')                            #========================================================HERE
names_wht_data = names_wht[1].data


flux_250_observed_pre= names_wht_data.field('F250')    #========================================================HERE(Boote.s: F250, et_F250; EG.S/Lockma.n: f250, et250; G1.2/NG.P: F250, E250)
error_250_observed_pre= names_wht_data.field('et_f250')    #========================================================HERE
source_name_pre= names_wht_data.field('name')


flux_250_observed_pre = np.array(flux_250_observed_pre)
flux_250_observed_pre=np.array([float(i) for i in flux_250_observed_pre])


source_name=[]

for i in range(0, len(flux_250_observed_pre)):

    if  math.isnan(flux_250_observed_pre[i]) == False:

        source_name.append(source_name_pre[i])


print source_name






#======================================extract source names===============================================================================================================

#for l in range(0, len(cluster_list)):
for l in range(0, 2):


    template_file = open("./MCMC_results/180503_combine/"+source_name[l]+"_mcmc_values.cat", "r")    #======================================================================HERE: change template


    lines = template_file.readlines()
    template_file.close()

    response_curve=[]

    for i in range(0, len(lines)):
        separated_lines=lines[i].split() 
        response_curve.append(separated_lines)


    response_curve = np.array(response_curve)
    walker_z=response_curve[:,0]
    walker_a=response_curve[:,1]


    walker_z = np.array(walker_z)
    walker_z=np.array([float(i) for i in walker_z])  #(wavelength in angstrom)
    walker_a = np.array(walker_a)
    walker_a=np.array([float(i) for i in walker_a])  #(wavelength in angstrom)



    hist, bin_edges =np.histogram(walker_z, bins=20, range=(0.0, 10.0), normed=True)

    print hist
    print bin_edges

    plt.hist(walker_z, bins=20, range=(0.0, 10.0), normed=1)
    plt.show()










