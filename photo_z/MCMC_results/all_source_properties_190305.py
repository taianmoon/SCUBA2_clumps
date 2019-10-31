import numpy as np
import corner
import matplotlib.pyplot as plt
import pyfits
import math


want_property='L'     #===============================================================================================HERE:'z','a','L',('SFR')
cat1_true= False     #===============================================================================================HERE
Todd_switch=True     #===============================================================================================HERE


if cat1_true==True:
    map_size_all= 6.0*math.pi*pow(350.0/3600.0, 2.0)    #in deg^2======================================================================HERE: need to calculate if not from our SCUBA-2 maps!!
else:
    map_size_all= 13.0*math.pi*pow(350.0/3600.0, 2.0)    #in deg^2======================================================================HERE: need to calculate if not from our SCUBA-2 maps!!

print map_size_all

map_size_all_S2CLS_random= 1802.0*1802.0*5.5555555555555E-4*0.00055555555555555*(20.0/507.0)  #==================================================S2CLS_random
map_size_all_S2CLS=   1802.0*1802.0*5.5555555555555E-4*0.00055555555555555 #==================================================S2CLS
map_size_all_PCL1002= math.pi*pow(350.0/3600.0, 2.0) #===========================================================================================PCL1002
map_size_all_Todd=2.31495185185

if cat1_true==True:
    template_file = open("./marginalized/allsource_"+want_property+"_cutz6_cat1.cat", "r")    #==================================================================================================HERE
else:
    template_file = open("./marginalized/allsource_"+want_property+"_cutz6.cat", "r")    #==================================================================================================HERE
lines = template_file.readlines()
template_file.close()
response_curve=[]
for i in range(0, len(lines)):
    separated_lines=lines[i].split() 
    response_curve.append(separated_lines)
response_curve = np.array(response_curve)
distribution_all=response_curve[:,0]      
distribution_all = np.array(distribution_all)
distribution_all=np.array([float(i) for i in distribution_all])  #(wavelength in angstrom)
print distribution_all


template_file_Todd = open("./marginalized/"+want_property+"_distribution_Todd.cat", "r")    #==================================================================================================HERE
lines_Todd = template_file_Todd.readlines()
template_file_Todd.close()
response_curve_Todd=[]
for i in range(0, len(lines_Todd)):
    separated_lines_Todd=lines_Todd[i].split() 
    response_curve_Todd.append(separated_lines_Todd)
response_curve_Todd = np.array(response_curve_Todd)
distribution_all_Todd=response_curve_Todd[:,0]      
distribution_all_Todd = np.array(distribution_all_Todd)
distribution_all_Todd=np.array([float(i) for i in distribution_all_Todd])  #(wavelength in angstrom)
print distribution_all_Todd


template_file_PCL1002 = open("./marginalized/allsource_"+want_property+"_PCL1002_cutz6.cat", "r")    #==================================================================================================HERE
lines_PCL1002 = template_file_PCL1002.readlines()
template_file_PCL1002.close()
response_curve_PCL1002=[]
for i in range(0, len(lines_PCL1002)):
    separated_lines_PCL1002=lines_PCL1002[i].split() 
    response_curve_PCL1002.append(separated_lines_PCL1002)
response_curve_PCL1002 = np.array(response_curve_PCL1002)
distribution_all_PCL1002=response_curve_PCL1002[:,0]      
distribution_all_PCL1002 = np.array(distribution_all_PCL1002)
distribution_all_PCL1002=np.array([float(i) for i in distribution_all_PCL1002])  #(wavelength in angstrom)
print distribution_all_PCL1002


if want_property=='z':
    #template_file_S2CLS_random = open("./marginalized/allsource_"+want_property+"_S2CLS_cutz6.cat", "r")    #==============================================================================HERE
    #lines_S2CLS_random = template_file_S2CLS_random.readlines()
    #template_file_S2CLS_random.close()
    #response_curve_S2CLS_random=[]
    #for i in range(0, len(lines_S2CLS_random)):
    #    separated_lines_S2CLS_random=lines_S2CLS_random[i].split() 
    #response_curve_S2CLS_random.append(separated_lines_S2CLS_random)
    #response_curve_S2CLS_random = np.array(response_curve_S2CLS_random)
    #distribution_all_S2CLS_random=response_curve_S2CLS_random[:,0]      
    #distribution_all_S2CLS_random = np.array(distribution_all_S2CLS_random)
    #distribution_all_S2CLS_random=np.array([float(i) for i in distribution_all_S2CLS_random])  #(wavelength in angstrom)
    #print distribution_all_S2CLS_random
    names_wht_S2CLS = pyfits.open('../S2CLS/photoz_estimate_cat_S2CLS.fits') 
    names_wht_data_S2CLS = names_wht_S2CLS[1].data  
    distribution_all_S2CLS_tmp= names_wht_data_S2CLS.field('best_redshift_marginalized')
    distribution_all_S2CLS_tmp = np.array(distribution_all_S2CLS_tmp)
    distribution_all_S2CLS_tmp=np.array([float(i) for i in distribution_all_S2CLS_tmp])
    distribution_all_S2CLS=[]
    for i in range(0, len(distribution_all_S2CLS_tmp)):
        if distribution_all_S2CLS_tmp[i] < 6.0:
            distribution_all_S2CLS.append(distribution_all_S2CLS_tmp[i])
    distribution_all_S2CLS = np.array(distribution_all_S2CLS)
    distribution_all_S2CLS=np.array([float(i) for i in distribution_all_S2CLS])

#----------------------Plot 1: plt.hist-------------------------------------------------------
#fig = plt.figure()

#axs = fig.gca()

#plt.hist(distribution_all, bins=10, range=[1.5, 4.0])  # arguments are passed to np.histogram

#plt.hist(distribution_all, bins=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])  # arguments are passed to np.histogram


#axs.set_yticklabels(plt.yticks()[0]/(map_size_all*3600.0))




#plt.title("Histogram with 'auto' bins")
#plt.show()




#----------------------Plot 2: np.histogram-------------------------------------------------------
if want_property=='z':
    hist, bin_edges = np.histogram(distribution_all, bins=6, range=[0.0, 6.0])  #==================================================================================================HERE:z
    hist_PCL1002, bin_edges_PCL1002 = np.histogram(distribution_all_PCL1002, bins=6, range=[0.0, 6.0])  #==========================================================================HERE:z
    hist_S2CLS_random, bin_edges_S2CLS_random = np.histogram(distribution_all_S2CLS, bins=6, range=[0.0, 6.0])  #===========================================HERE:a
    hist_Todd, bin_edges_Todd = np.histogram(distribution_all_Todd, bins=6, range=[0.0, 6.0])  #===========================================HERE:a
if want_property=='a':
    hist, bin_edges = np.histogram(distribution_all, bins=10, range=[-2.0, 2.0])  #==================================================================================================HERE:a
    hist_PCL1002, bin_edges_PCL1002 = np.histogram(distribution_all_PCL1002, bins=10, range=[-2.0, 2.0])  #=========================================================================HERE:a
    #hist_S2CLS_random, bin_edges_S2CLS_random = np.histogram(distribution_all_S2CLS, bins=10, range=[-2.0, 2.0])  #================================================HERE:a
if want_property=='L':
    hist, bin_edges = np.histogram(distribution_all, bins=8, range=[12.0, 14.0])  #==================================================================================================HERE:L
    hist_PCL1002, bin_edges_PCL1002 = np.histogram(distribution_all_PCL1002, bins=8, range=[12.0, 14.0])  #=========================================================================HERE:L
    #hist_S2CLS_random, bin_edges_S2CLS_random = np.histogram(distribution_all_S2CLS, bins=8, range=[12.0, 14.0])  #==========================================================HERE:L
    hist_Todd, bin_edges_Todd = np.histogram(distribution_all_Todd, bins=8, range=[12.0, 14.0])  #=========================================================================HERE:L
if want_property=='SFR':
    hist, bin_edges = np.histogram(distribution_all, bins=10, range=[1.5, 4.0])  #==================================================================================================HERE:SFR
    hist_PCL1002, bin_edges_PCL1002 = np.histogram(distribution_all_PCL1002, bins=10, range=[1.5, 4.0])  #=========================================================================HERE:SFR
    #hist_S2CLS_random, bin_edges_S2CLS_random = np.histogram(distribution_all_S2CLS, bins=10, range=[1.5, 4.0])  #========================================================HERE:SFR
    hist_Todd, bin_edges_Todd = np.histogram(distribution_all_Todd, bins=10, range=[1.5, 4.0])  #=========================================================================HERE:SFR


#hist = np.array(hist)
#hist=np.array([float(i) for i in hist])

#print hist

bin_edges = bin_edges[1:]
#bin_edges_PCL1002 =bin_edges_PCL1002[1:]
#bin_edges_S2CLS_random=bin_edges_S2CLS_random[1:]


bin_edges_true=[]
hist_true=[]
#bin_edges_true_PCL1002=[]
hist_true_PCL1002=[]
#bin_edges_true_S2CLS_random=[]
hist_true_S2CLS_random=[]
hist_true_Todd=[]
for i in range(0, len(bin_edges)):
    #bin_edges_true.append(bin_edges[i]+0.5)
    bin_edges_true.append(bin_edges[i])
    #hist_true.append(hist[i])   #divide by map size to get number per arcmin^2
    hist_true.append(hist[i]/(map_size_all*3600.0))   #divide by map size to get number per arcmin^2  ===================================HERE: becareful of other fields!
    hist_true_PCL1002.append(hist_PCL1002[i]/(map_size_all_PCL1002*3600.0))   #divide by map size to get number per arcmin^2  ===================================HERE: becareful of other fields!
    #hist_true_PCL1002.append(hist_PCL1002[i])   #divide by map size to get number per arcmin^2  ===================================HERE: becareful of other fields!
    if want_property=='z':
        hist_true_S2CLS_random.append(hist_S2CLS_random[i]/(map_size_all_S2CLS*3600.0))   #divide by map size to get number per arcmin^2  =============================HERE: becareful of other fields!
        #hist_true_S2CLS_random.append(hist_S2CLS_random[i])   #divide by map size to get number per arcmin^2  =============================HERE: becareful of other fields!
        hist_true_Todd.append(hist_Todd[i]/(map_size_all_Todd*3600.0))   #divide by map size to get number per arcmin^2  =============================HERE: becareful of other fields!
    if want_property=='L' or want_property=='SFR':
        hist_true_Todd.append(hist_Todd[i]/(map_size_all_Todd*3600.0))   #divide by map size to get number per arcmin^2  =============================HERE: becareful of other fields!


#print bin_edges_true
#print map_size_all*3600.0


fig = plt.figure(figsize=(15,10))
#plt.plot(bin_edges_true, hist_true, linewidth=1.5, label='This work')
ax1 = fig.add_subplot(111)

ax1.step(bin_edges_true, hist_true_PCL1002, color='cyan', linewidth=4, label='PCL1002')
ax1.step(bin_edges_true, hist_true, color='black', linewidth=4, label='This work')
if want_property=='z':
    ax1.step(bin_edges_true, hist_true_S2CLS_random, color='red', linewidth=4, label='S2CLS')
    ax1.step(bin_edges_true, hist_true_Todd, color='magenta', linewidth=4, label='MacKenzie 2017')
if want_property=='L' or want_property=='SFR':


    ax1.step(bin_edges_true, hist_true_Todd, color='magenta', linewidth=4, label='MacKenzie 2017')
    ax2 = ax1.twiny()

    #new_tick_locations = np.array([11.9698700396, 12.4698700396, 12.9698700396, 13.4698700396, 13.9698700396])  #Use log10(SFR) -log10(2.8e-44*3.828e26*1.0e7)) to convert into log10(L(L_solar))
    #new_tick_locations = np.array([11.8270784696, 12.3270784696, 12.8270784696, 13.3270784696, 13.8270784696])  #Use log10(SFR) -log10(3.89e-44*3.828e26*1.0e7)) to convert into log10(L(L_solar))
    new_tick_locations = np.array([12.3270784696, 12.8270784696, 13.3270784696, 13.8270784696])  #Use log10(SFR) -log10(3.89e-44*3.828e26*1.0e7)) to convert into log10(L(L_solar))
    def tick_function(X):
        #V = X+ math.log10(2.8e-44*3.828e26*1.0e7)
        V = X+ math.log10(3.89e-44*3.828e26*1.0e7)
        return ["%.1f" % z for z in V]

    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(tick_function(new_tick_locations))
    ax2.set_xlabel(r'$log(SFR(M_{\odot}/yr))$')











#plt.scatter(bin_edges_true, hist_true, s=80.0)

#plt.grid()
ax1.legend(loc=1)

if want_property=='z':
    ax1.set_xlabel('Redshift')  #==================================================================================================HERE:z
if want_property=='L':
    ax1.set_xlabel(r'$log(L(L_{\odot}))$')  #==================================================================================================HERE:L
if want_property=='a':
    ax1.set_xlabel(r'$log(a)$')  #==================================================================================================HERE:a
if want_property=='SFR':
    ax1.set_xlabel(r'$log(SFR(M_{\odot}/yr))$')  #==================================================================================================HERE:SFR

ax1.set_ylabel(r'Source density ($\mathrm{arcmin}^{-2}$)')
#plt.ylabel(r'Number of Sources')
plt.rc('font', size=30)
plt.show()



