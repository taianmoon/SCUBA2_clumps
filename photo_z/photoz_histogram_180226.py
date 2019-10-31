import matplotlib.pyplot as plt
import numpy as np
import math
#from astropy.io import fits
#import pyfits


#===============Load in data===============================================

template_file = open("./S2COSMOS/S2COSMOS_hist.txt", "r")

lines = template_file.readlines()
template_file.close()

response_curve=[]

for i in range(0, len(lines)):
    separated_lines=lines[i].split() 
    response_curve.append(separated_lines)


response_curve = np.array(response_curve)
redshift_distribution=response_curve[0,:]



redshift_distribution = np.array(redshift_distribution)
redshift_distribution=np.array([float(i) for i in redshift_distribution])  #(wavelength in angstrom)


print redshift_distribution


#==================Put in models===========================================

template_file_model = open("./S2COSMOS/Chapman_smg_specz.txt", "r")

lines_model = template_file_model.readlines()
template_file_model.close()

response_curve_model=[]

for i in range(0, len(lines_model)):
    separated_lines_model=lines_model[i].split() 
    response_curve_model.append(separated_lines_model)


response_curve_model = np.array(response_curve_model)
redshift_distribution_model=response_curve_model[:,0]



redshift_distribution_model = np.array(redshift_distribution_model)
redshift_distribution_model=np.array([float(i) for i in redshift_distribution_model])  #(wavelength in angstrom)



print redshift_distribution_model

#=================Plotting=================================================


n_data, bins_data, patches_data = plt.hist(redshift_distribution, label='S2CLS', alpha=0.4, bins=np.arange(min(redshift_distribution), max(redshift_distribution) + 0.5, 0.5), histtype='stepfilled')                     
#n_model, bins_model, patches_model = plt.hist(redshift_distribution_model, label='Chapman 2005', alpha=0.4, bins=np.arange(min(redshift_distribution_model), max(redshift_distribution_model) + 0.5, 0.5), histtype='stepfilled', normed=1) 
n_model, bins_model, patches_model = plt.hist(redshift_distribution_model, label='Chapman 2005', alpha=0.4, bins=9, range=[0.0, 3.87], histtype='stepfilled')                                         

print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
#print n
#print bins
#print patches
#print max(n_data)
#print max(n_model)

#ratio=max([max(n_data),max(n_model)])/min([max(n_data),max(n_model)])

#print ratio

n_data_normalized=n_data/max(n_data)
n_model_normalized=n_model/max(n_model)

print n_data_normalized
print bins_data
print n_model_normalized
print bins_model


#plt.plot(n_data_normalized,bins_data)


#plt.title("S2COSMOS") #===============================================================================================================================HERE
plt.xlabel("Photometric Redshift")
plt.ylabel("Number of Sources")

plt.legend()
plt.rc('font', size=15)

#plt.show()

#fig2.set_size_inches(20,10)

#fig2.savefig('./S2COSMOS/S2COSMOS_hist.jpg')  #===============================================================================================================================HERE

plt.close()



n_data_normalized_plot=[]
for i in range(0, len(n_data_normalized)):

    n_data_normalized_plot.append(n_data_normalized[i])
    n_data_normalized_plot.append(n_data_normalized[i])


print n_data_normalized_plot


bins_data_plot=[]
for i in range(0, len(bins_data)):

    bins_data_plot.append(bins_data[i])

    if i != 0 and i != len(bins_data)-1:
        bins_data_plot.append(bins_data[i])


print bins_data_plot



n_model_normalized_plot=[]
for i in range(0, len(n_model_normalized)):

    n_model_normalized_plot.append(n_model_normalized[i])
    n_model_normalized_plot.append(n_model_normalized[i])


print n_model_normalized_plot


bins_model_plot=[]
for i in range(0, len(bins_model)):

    bins_model_plot.append(bins_model[i])

    if i != 0 and i != len(bins_model)-1:
        bins_model_plot.append(bins_model[i])


print bins_model_plot




plt.plot(bins_data_plot, n_data_normalized_plot, label='S2CLS', linewidth=3.0)
plt.plot(bins_model_plot, n_model_normalized_plot, label='Chapman 2005', linewidth=3.0)


plt.ylim(top=1.1)
plt.xlabel("Photometric Redshift")
plt.ylabel("Number of Sources (Normalized)")

plt.legend()
plt.rc('font', size=30)

plt.show()






