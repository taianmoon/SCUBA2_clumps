#++++++++++++++++++++++++++++++++++ALESS+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template_file_aless = open("./aless_average_seds.dat", "r")    #========================================================HERE

lines_aless = template_file_aless.readlines()
template_file_aless.close()

response_curve_aless=[]

for i in range(0, len(lines_aless)):
    separated_lines_aless=lines_aless[i].split() 
    response_curve_aless.append(separated_lines_aless)


response_curve_aless = np.array(response_curve_aless)
wavelength_micron_aless=response_curve_aless[:,0]
flux_mjy_aless=response_curve_aless[:,1]


wavelength_micron_aless = np.array(wavelength_micron_aless)
flux_mjy_aless = np.array(flux_mjy_aless)
wavelength_micron_aless=np.array([float(i) for i in wavelength_micron_aless])  #(wavelength in angstrom)
flux_mjy_aless=np.array([float(i) for i in flux_mjy_aless])  #(wavelength in angstrom)



#flux_250_z_aless=[]
flux_350_z_aless=[]
flux_500_z_aless=[]
flux_850_z_aless=[]

#z_list=np.arange(0.0, 5.0, 0.1)  #========================================================HERE: change parameter grids
#dely_list=np.arange(0.0, 4.0, 0.05)  #========================================================HERE: change parameter grids

#print z_list
#print dely_list


#z_list_space=[]
#dely_list_space=[]
for i in range(0, len(z_list)):
    for j in range(0, len(dely_list)):

        wavelength_micron_z_aless = wavelength_micron_aless*(1.0+z_list[i])
        flux_mjy_z_aless = flux_mjy_aless*dely_list[j]
        #flux_250_z_aless.append(np.interp(250.0, wavelength_micron_z_aless, flux_mjy_z_aless))
        flux_350_z_aless.append(np.interp(350.0, wavelength_micron_z_aless, flux_mjy_z_aless))
        flux_500_z_aless.append(np.interp(550.0, wavelength_micron_z_aless, flux_mjy_z_aless))  #==============================190730 new!
        flux_850_z_aless.append(np.interp(850.0, wavelength_micron_z_aless, flux_mjy_z_aless))
        #z_list_space.append(z_list[i])
        #dely_list_space.append(dely_list[j])

#z_list_space = np.array(z_list_space)
#z_list_space=np.array([float(i) for i in z_list_space])  #(wavelength in angstrom)
#dely_list_space = np.array(dely_list_space)
#dely_list_space=np.array([float(i) for i in dely_list_space])  #(wavelength in angstrom)





#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++HFLS3++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



template_file_hfls3 = open("./HFLS3.txt", "r")    #========================================================HERE

lines_hfls3 = template_file_hfls3.readlines()
template_file_hfls3.close()

response_curve_hfls3=[]

for i in range(0, len(lines_hfls3)):
    separated_lines_hfls3=lines_hfls3[i].split() 
    response_curve_hfls3.append(separated_lines_hfls3)


response_curve_hfls3 = np.array(response_curve_hfls3)
wavelength_micron_hfls3=response_curve_hfls3[:,0]
flux_mjy_hfls3=response_curve_hfls3[:,1]


wavelength_micron_hfls3 = np.array(wavelength_micron_hfls3)
flux_mjy_hfls3 = np.array(flux_mjy_hfls3)
wavelength_micron_hfls3=np.array([float(i) for i in wavelength_micron_hfls3])  #(wavelength in angstrom)
flux_mjy_hfls3=np.array([float(i) for i in flux_mjy_hfls3])  #(wavelength in angstrom)


wavelength_micron_hfls3=wavelength_micron_hfls3/7.34  #redshift back to restframe (z=6.34)
wavelength_micron_hfls3=wavelength_micron_hfls3*1.0e6  #convert from m to micron
flux_mjy_hfls3=flux_mjy_hfls3*1000.0  #convert from Jy to mJy


#flux_250_z_hfls3=[]
flux_350_z_hfls3=[]
flux_500_z_hfls3=[]
flux_850_z_hfls3=[]

#z_list=np.arange(0.0, 5.0, 0.1)  #========================================================HERE: change parameter grids
#dely_list=np.arange(0.0, 4.0, 0.05)  #========================================================HERE: change parameter grids

#print z_list
#print dely_list


#z_list_space=[]
#dely_list_space=[]
for i in range(0, len(z_list)):
    for j in range(0, len(dely_list)):

        wavelength_micron_z_hfls3 = wavelength_micron_hfls3*(1.0+z_list[i])
        flux_mjy_z_hfls3 = flux_mjy_hfls3*dely_list[j]
        #flux_250_z_hfls3.append(np.interp(250.0, wavelength_micron_z_hfls3, flux_mjy_z_hfls3))
        flux_350_z_hfls3.append(np.interp(350.0, wavelength_micron_z_hfls3, flux_mjy_z_hfls3))
        flux_500_z_hfls3.append(np.interp(550.0, wavelength_micron_z_hfls3, flux_mjy_z_hfls3))  #==============================190730 new!
        flux_850_z_hfls3.append(np.interp(850.0, wavelength_micron_z_hfls3, flux_mjy_z_hfls3))
        #z_list_space.append(z_list[i])
        #dely_list_space.append(dely_list[j])

#z_list_space = np.array(z_list_space)
#z_list_space=np.array([float(i) for i in z_list_space])  #(wavelength in angstrom)
#dely_list_space = np.array(dely_list_space)
#dely_list_space=np.array([float(i) for i in dely_list_space])  #(wavelength in angstrom)










#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++Cosmic eyelash++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



#template_file_eyelash = open("./cosmic_eyelash_sed_lapi.txt", "r")    #========================================================HERE
template_file_eyelash = open("./cosmic_eyelash_sed_josh.txt", "r")    #========================================================HERE

lines_eyelash = template_file_eyelash.readlines()
template_file_eyelash.close()

response_curve_eyelash=[]

for i in range(0, len(lines_eyelash)):
    separated_lines_eyelash=lines_eyelash[i].split() 
    response_curve_eyelash.append(separated_lines_eyelash)


response_curve_eyelash = np.array(response_curve_eyelash)
wavelength_micron_eyelash=response_curve_eyelash[:,0]
flux_mjy_eyelash=response_curve_eyelash[:,1]


wavelength_micron_eyelash = np.array(wavelength_micron_eyelash)
flux_mjy_eyelash = np.array(flux_mjy_eyelash)
wavelength_micron_eyelash=np.array([float(i) for i in wavelength_micron_eyelash])  #(wavelength in angstrom)
flux_mjy_eyelash=np.array([float(i) for i in flux_mjy_eyelash])  #(wavelength in angstrom)



#flux_mjy_loop=[]
#for k in range(0, len(wavelength_micron_eyelash)):

#    flux_mjy_loop.append(1.0e-14*1000.0*(1.0e26)*flux_mjy_eyelash[k]*(wavelength_micron_eyelash[k]*1.0e-6)/2.99792458e8)


#flux_mjy_eyelash=flux_mjy_loop

#flux_mjy_eyelash = np.array(flux_mjy_eyelash)
#flux_mjy_eyelash=np.array([float(i) for i in flux_mjy_eyelash])  #(wavelength in angstrom)





#flux_250_z_eyelash=[]
flux_350_z_eyelash=[]
flux_500_z_eyelash=[]
flux_850_z_eyelash=[]

#z_list=np.arange(0.0, 5.0, 0.1)  #========================================================HERE: change parameter grids
#dely_list=np.arange(0.0, 4.0, 0.05)  #========================================================HERE: change parameter grids

#print z_list
#print dely_list


#z_list_space=[]
#dely_list_space=[]
for i in range(0, len(z_list)):
    for j in range(0, len(dely_list)):

        wavelength_micron_z_eyelash = wavelength_micron_eyelash*(1.0+z_list[i])
        flux_mjy_z_eyelash = flux_mjy_eyelash*dely_list[j]
        #flux_250_z_eyelash.append(np.interp(250.0, wavelength_micron_z_eyelash, flux_mjy_z_eyelash))
        flux_350_z_eyelash.append(np.interp(350.0, wavelength_micron_z_eyelash, flux_mjy_z_eyelash))
        flux_500_z_eyelash.append(np.interp(550.0, wavelength_micron_z_eyelash, flux_mjy_z_eyelash))   #==============================190730 new!
        flux_850_z_eyelash.append(np.interp(850.0, wavelength_micron_z_eyelash, flux_mjy_z_eyelash))
        #z_list_space.append(z_list[i])
        #dely_list_space.append(dely_list[j])










