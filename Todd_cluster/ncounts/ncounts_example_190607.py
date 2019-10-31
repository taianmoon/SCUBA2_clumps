import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm
from astropy.io import fits
import astropy.wcs as wcs
from scipy.interpolate import interp1d


#cluster_list=['PLCK_G49.6-42.9', 'PLCK_G082.5+38.4', 'PLCK_G165.8+45.3']

cluster_list=['Planck18p194', 'Planck18p735', 'Planck24p194', 'PLCK_DU_G045.7-41.2', 'PLCK_DU_G059.1-67.1', 'PLCK_DU_G073.4-57.5', 'PLCK_G006.1+61.8', 'PLCK_G009.8+72.6', 'PLCK_G056.7+62.6', 'PLCK_G068.3+31.9', 'PLCK_G075.1+33.2', 'PLCK_G077.7+32.6', 'PLCK_G078.9+48.2', 'PLCK_G082.5+38.4', 'PLCK_G083.3+51.0', 'PLCK_G091.9+43.0', 'PLCK_G093.6+55.9', 'PLCK_G132.9-76.0', 'PLCK_G144.1+81.0', 'PLCK_G160.7+41.0', 'PLCK_G162.1-59.3', 'PLCK_G165.8+45.3', 'PLCK_G173.8+59.3', 'PLCK_G177.0+35.9', 'PLCK_G179.3+50.7', 'PLCK_G186.3-72.7', 'PLCK_G186.6+66.7', 'PLCK_G188.6-68.9', 'PLCK_G191.3+62.0', 'PLCK_G191.8-83.4', 'PLCK_G201.1+50.7', 'PLCK_G213.0+65.9', 'PLCK_G223.9+41.2', 'PLCK_G328.9+71.4', 'PLCK_G49.6-42.9', 'PLCK_G84.0-71.5', 'PLCK_HZ_G038.0-51.5', 'PLCK_HZ_G067.2-63.8', 'PLCK_HZ_G103.1-73.6', 'PLCK_HZ_G106.8-83.3', 'PLCK_HZ_G119.4-76.6', 'PLCK_HZ_G132.6-81.1', 'PLCK_HZ_G171.1-78.7', 'PLCK_HZ_G173.9+57.0', 'PLCK_HZ_G176.6+59.0', 'PLCK_HZ_G214.1+48.3']


template_file = open("./cum_counts_geach_bin2.cat", "r")   #===================================================HERE

lines = template_file.readlines()[1:]
#lines = template_file.readlines()
template_file.close()

response_curve=[]

for i in range(0, len(lines)):
    separated_lines=lines[i].split() 
    response_curve.append(separated_lines)


response_curve = np.array(response_curve)
geach_flux = response_curve[:,0]      #=========================================================HERE: 7 for Todd's cluster, 6 for our clusters
geach_cum_number_plot = response_curve[:,1]    #=========================================================HERE: 8 for Todd's cluster, 7 for our clusters
geach_cum_number_err_low_plot = response_curve[:,2]     #=========================================================HERE: 5 for Todd's cluster, 4 for our clusters
geach_cum_number_err_up_plot = response_curve[:,3]     #=========================================================HERE: 5 for Todd's cluster, 4 for our clusters


geach_flux = np.array(geach_flux)
geach_cum_number_plot = np.array(geach_cum_number_plot)
geach_cum_number_err_low_plot = np.array(geach_cum_number_err_low_plot)
geach_cum_number_err_up_plot = np.array(geach_cum_number_err_up_plot)
geach_flux = np.array([float(i) for i in geach_flux])  #(wavelength in angstrom)
geach_cum_number_plot = np.array([float(i) for i in geach_cum_number_plot])  #(wavelength in angstrom)
geach_cum_number_err_low_plot = np.array([float(i) for i in geach_cum_number_err_low_plot])  #(wavelength in angstrom)
geach_cum_number_err_up_plot = np.array([float(i) for i in geach_cum_number_err_up_plot])  #(wavelength in angstrom)

asymmetric_error_plot = [geach_cum_number_err_low_plot, geach_cum_number_err_up_plot]

print geach_flux
print geach_cum_number_plot
print geach_cum_number_err_low_plot
print geach_cum_number_err_up_plot



for l in range(0, len(cluster_list)):
    print cluster_list[l]


    template_file = open("./cum_190813/cum_counts_"+cluster_list[l]+"_bin2.cat", "r")   #===================================================HERE

    lines = template_file.readlines()[1:]
    #lines = template_file.readlines()
    template_file.close()

    response_curve=[]

    for i in range(0, len(lines)):
        separated_lines=lines[i].split() 
        response_curve.append(separated_lines)


    response_curve = np.array(response_curve)

    if response_curve.shape == (0,):
        bin_flux=[]
        n_cum=[]
        completeness_err=[]
    else:
        bin_flux = response_curve[:,0]      #=========================================================HERE: 7 for Todd's cluster, 6 for our clusters
        n_cum = response_curve[:,1]    #=========================================================HERE: 8 for Todd's cluster, 7 for our clusters
        completeness_err = response_curve[:,2]     #=========================================================HERE: 5 for Todd's cluster, 4 for our clusters


    bin_flux = np.array(bin_flux)
    n_cum = np.array(n_cum)
    completeness_err = np.array(completeness_err)
    bin_flux = np.array([float(i) for i in bin_flux])  #(wavelength in angstrom)
    n_cum = np.array([float(i) for i in n_cum])  #(wavelength in angstrom)
    completeness_err = np.array([float(i) for i in completeness_err])  #(wavelength in angstrom)

    print bin_flux
    print n_cum
    print completeness_err

    if cluster_list[l]== 'PLCK_G49.6-42.9':
        color='blue'
        plt.plot(bin_flux, n_cum, label=cluster_list[l], color=color, linewidth=3.0)
        plt.scatter(bin_flux, n_cum, s=80.0, marker='s', color=color)
        plt.errorbar(bin_flux, n_cum, yerr=completeness_err, fmt='o', capsize=5, elinewidth=3,markeredgewidth=3, color=color)
    elif cluster_list[l]=='PLCK_G201.1+50.7':
        color='green'
        plt.plot(bin_flux, n_cum, label=cluster_list[l], color=color, linewidth=3.0)
        plt.scatter(bin_flux, n_cum, s=80.0, marker='s', color=color)
        plt.errorbar(bin_flux, n_cum, yerr=completeness_err, fmt='o', capsize=5, elinewidth=3,markeredgewidth=3, color=color)
    elif cluster_list[l]=='PLCK_G165.8+45.3':
        color='purple'
        plt.plot(bin_flux, n_cum, label=cluster_list[l], color=color, linewidth=3.0)
        plt.scatter(bin_flux, n_cum, s=80.0, marker='s', color=color)
        plt.errorbar(bin_flux, n_cum, yerr=completeness_err, fmt='o', capsize=5, elinewidth=3,markeredgewidth=3, color=color)

    plt.plot(bin_flux, n_cum, color='lightgray', zorder=1)


    plt.tick_params(width=2, length=16, which='major')
    plt.tick_params(width=2, length=5, which='minor')


    plt.xscale('log',nonposy='clip')
    plt.yscale('log',nonposy='clip')

    plt.xlim((3.0, 30.0))
    plt.ylim((1e-2, 1e2))


    #leg = plt.legend(handlelength=0, handletextpad=0, fancybox=True, loc='upper right')
    ##print leg
    ##print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    #for item in leg.legendHandles:
    #    item.set_visible(False)
    #leg.get_frame().set_linewidth(0.0)
    #leg.get_frame().set_alpha(0.0)



    plt.grid()

    plt.xlabel('Flux Density (mJy)')  
    plt.ylabel(r'Cumulative Count' +'\n'+'(per map size)')  
    plt.rc('font', size=18)


plt.plot(geach_flux, geach_cum_number_plot, color='red', label='Geach 17 (field)', linewidth=3.0)
plt.scatter(geach_flux, geach_cum_number_plot, color='red')
plt.errorbar(geach_flux, geach_cum_number_plot, yerr=asymmetric_error_plot, color='red',fmt='o', capsize=5, elinewidth=3,markeredgewidth=3)

plt.legend(loc=3)
plt.show()




