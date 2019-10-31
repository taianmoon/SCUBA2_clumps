import numpy as np
import math


#cluster_list=['Planck18p194', 'Planck18p735', 'Planck24p194', 'PLCK_DU_G045.7-41.2', 'PLCK_DU_G059.1-67.1', 'PLCK_DU_G073.4-57.5', 'PLCK_G006.1+61.8', 'PLCK_G009.8+72.6', 'PLCK_G056.7+62.6', 'PLCK_G068.3+31.9', 'PLCK_G075.1+33.2', 'PLCK_G077.7+32.6', 'PLCK_G078.9+48.2', 'PLCK_G082.5+38.4', 'PLCK_G083.3+51.0', 'PLCK_G091.9+43.0', 'PLCK_G093.6+55.9', 'PLCK_G132.9-76.0', 'PLCK_G144.1+81.0', 'PLCK_G160.7+41.0', 'PLCK_G162.1-59.3', 'PLCK_G165.8+45.3', 'PLCK_G173.8+59.3', 'PLCK_G177.0+35.9', 'PLCK_G179.3+50.7', 'PLCK_G186.3-72.7', 'PLCK_G186.6+66.7', 'PLCK_G188.6-68.9', 'PLCK_G191.3+62.0', 'PLCK_G191.8-83.4', 'PLCK_G201.1+50.7', 'PLCK_G213.0+65.9', 'PLCK_G223.9+41.2', 'PLCK_G328.9+71.4', 'PLCK_G49.6-42.9', 'PLCK_G84.0-71.5', 'PLCK_HZ_G038.0-51.5', 'PLCK_HZ_G067.2-63.8', 'PLCK_HZ_G103.1-73.6', 'PLCK_HZ_G106.8-83.3', 'PLCK_HZ_G119.4-76.6', 'PLCK_HZ_G132.6-81.1', 'PLCK_HZ_G171.1-78.7', 'PLCK_HZ_G173.9+57.0', 'PLCK_HZ_G176.6+59.0', 'PLCK_HZ_G214.1+48.3']


cluster_list=['Planck18p194', 'Planck18p735', 'PLCK_DU_G045.7-41.2', 'PLCK_DU_G059.1-67.1', 'PLCK_DU_G073.4-57.5', 'PLCK_G006.1+61.8', 'PLCK_G009.8+72.6', 'PLCK_G056.7+62.6', 'PLCK_G068.3+31.9', 'PLCK_G075.1+33.2', 'PLCK_G077.7+32.6', 'PLCK_G078.9+48.2', 'PLCK_G082.5+38.4', 'PLCK_G083.3+51.0', 'PLCK_G091.9+43.0', 'PLCK_G093.6+55.9', 'PLCK_G144.1+81.0', 'PLCK_G160.7+41.0', 'PLCK_G162.1-59.3', 'PLCK_G165.8+45.3', 'PLCK_G173.8+59.3', 'PLCK_G177.0+35.9', 'PLCK_G179.3+50.7', 'PLCK_G186.3-72.7', 'PLCK_G186.6+66.7', 'PLCK_G188.6-68.9', 'PLCK_G191.3+62.0', 'PLCK_G191.8-83.4', 'PLCK_G201.1+50.7', 'PLCK_G213.0+65.9', 'PLCK_G223.9+41.2', 'PLCK_G328.9+71.4', 'PLCK_G84.0-71.5', 'PLCK_HZ_G038.0-51.5', 'PLCK_HZ_G067.2-63.8', 'PLCK_HZ_G103.1-73.6', 'PLCK_HZ_G106.8-83.3', 'PLCK_HZ_G119.4-76.6', 'PLCK_HZ_G132.6-81.1', 'PLCK_HZ_G171.1-78.7', 'PLCK_HZ_G173.9+57.0', 'PLCK_HZ_G176.6+59.0', 'PLCK_HZ_G214.1+48.3']


for l in range(0, len(cluster_list)):

    print l
    print cluster_list[l]


    template_file = open("./for_code/"+cluster_list[l]+"_3p5sigma_negative.cat", "r")    #=========================================HERE

    lines = template_file.readlines()[1:]
    #lines = template_file.readlines()
    template_file.close()

    response_curve=[]

    for i in range(0, len(lines)):
        separated_lines=lines[i].split() 
        response_curve.append(separated_lines)


    response_curve = np.array(response_curve)


    #name=response_curve[:,0]
    pix_x=response_curve[:,0]
    pix_y=response_curve[:,1]
    RA=response_curve[:,2]
    Dec=response_curve[:,3]
    SN_ratio=response_curve[:,4]
    flux_pW=response_curve[:,5]  #this is already in unit of mJy!
    #flux_mJy_perbeam=response_curve[:,6]
    flux_error_mJy_wrong=response_curve[:,7]
    #flux_mJy_perbeam_deboost=response_curve[:,9]
    #flux_error_mJy_deboost=response_curve[:,10]


    #name = np.array(name)
    #name=np.array([float(i) for i in name])  #(wavelength in angstrom)
    pix_x = np.array(pix_x)
    pix_x=np.array([float(i) for i in pix_x])  #(wavelength in angstrom)
    pix_y = np.array(pix_y)
    pix_y=np.array([float(i) for i in pix_y])  #(wavelength in angstrom)
    RA = np.array(RA)
    RA=np.array([float(i) for i in RA])  #(wavelength in angstrom)
    Dec = np.array(Dec)
    Dec=np.array([float(i) for i in Dec])  #(wavelength in angstrom)
    SN_ratio = np.array(SN_ratio)
    SN_ratio=np.array([float(i) for i in SN_ratio])  #(wavelength in angstrom)
    flux_pW = np.array(flux_pW)
    flux_pW=np.array([float(i) for i in flux_pW])  #(wavelength in angstrom)
    #flux_mJy_perbeam = np.array(flux_mJy_perbeam)
    #flux_mJy_perbeam=np.array([float(i) for i in flux_mJy_perbeam])  #(wavelength in angstrom)
    flux_error_mJy_wrong = np.array(flux_error_mJy_wrong)
    flux_error_mJy_wrong=np.array([float(i) for i in flux_error_mJy_wrong])  #(wavelength in angstrom)
    #flux_mJy_perbeam_deboost = np.array(flux_mJy_perbeam_deboost)
    #flux_mJy_perbeam_deboost=np.array([float(i) for i in flux_mJy_perbeam_deboost])  #(wavelength in angstrom)
    #flux_error_mJy_deboost = np.array(flux_error_mJy_deboost)
    #flux_error_mJy_deboost=np.array([float(i) for i in flux_error_mJy_deboost])  #(wavelength in angstrom)


    name=[]
    flux_mJy_perbeam=[]
    flux_error_mJy=[]
    flux_mJy_perbeam_deboost=[]
    flux_error_mJy_deboost=[]    

    flux_error_mJy_deboost_confusion=[]
    SN_ratio_confusion=[]
    flux_error_mJy_deboost_confusion_calibration=[]
    flux_error_mJy_deboost_calibration=[]


    for i in range(0, len(pix_x)):
        name.append(i)
        flux_mJy_perbeam.append(flux_pW[i])
        flux_error_mJy.append(flux_error_mJy_wrong[i]/(537.0*1000.0))  #convert back from the wrong conversion from pW to mJy
        flux_error_mJy_loop = flux_error_mJy_wrong[i]/(537.0*1000.0)
        flux_mJy_perbeam_deboost.append(flux_pW[i]/(1.0 + 0.2*pow(SN_ratio[i]/5.0, -2.3)))
        flux_error_mJy_deboost.append(flux_error_mJy_loop/(1.0 + 0.2*pow(SN_ratio[i]/5.0, -2.3)))
        flux_mJy_perbeam_deboost_loop = flux_pW[i]/(1.0 + 0.2*pow(SN_ratio[i]/5.0, -2.3))
        flux_error_mJy_deboost_loop = flux_error_mJy_loop/(1.0 + 0.2*pow(SN_ratio[i]/5.0, -2.3))


        flux_error_mJy_deboost_confusion.append(math.sqrt(pow(0.7, 2)+pow(flux_error_mJy_deboost_loop, 2)))
        SN_ratio_confusion.append(flux_mJy_perbeam_deboost_loop/(math.sqrt(pow(0.7, 2)+pow(flux_error_mJy_deboost_loop, 2))))
        flux_error_mJy_deboost_confusion_calibration.append(math.sqrt(pow(math.sqrt(pow(0.7, 2)+pow(flux_error_mJy_deboost_loop, 2)), 2)+pow(flux_mJy_perbeam_deboost_loop*0.05, 2)))
        flux_error_mJy_deboost_calibration.append(math.sqrt(pow(flux_error_mJy_deboost_loop, 2)+pow(flux_mJy_perbeam_deboost_loop*0.05, 2)))


    name = np.array(name)
    name=np.array([float(i) for i in name])  #(wavelength in angstrom)
    flux_mJy_perbeam = np.array(flux_mJy_perbeam)
    flux_mJy_perbeam=np.array([float(i) for i in flux_mJy_perbeam])  #(wavelength in angstrom)
    flux_error_mJy = np.array(flux_error_mJy)
    flux_error_mJy=np.array([float(i) for i in flux_error_mJy])  #(wavelength in angstrom)
    flux_mJy_perbeam_deboost = np.array(flux_mJy_perbeam_deboost)
    flux_mJy_perbeam_deboost=np.array([float(i) for i in flux_mJy_perbeam_deboost])  #(wavelength in angstrom)
    flux_error_mJy_deboost = np.array(flux_error_mJy_deboost)
    flux_error_mJy_deboost=np.array([float(i) for i in flux_error_mJy_deboost])  #(wavelength in angstrom)

    #print flux_error_mJy_deboost
    #print flux_error_mJy_deboost_confusion
    #print flux_error_mJy_deboost_confusion_calibration
    #print SN_ratio
    #print SN_ratio_confusion


    mock_catalogue=np.column_stack((name, pix_x, pix_y, RA, Dec, SN_ratio, flux_pW, flux_mJy_perbeam, flux_error_mJy, flux_mJy_perbeam_deboost, flux_error_mJy_deboost, flux_error_mJy_deboost_confusion, SN_ratio_confusion, flux_error_mJy_deboost_confusion_calibration, flux_error_mJy_deboost_calibration))

    np.savetxt('./for_code/updated/'+cluster_list[l]+'_3p5sigma_negative.cat', mock_catalogue, delimiter=' ', header="name pix_x pix_y RA Dec SN_ratio flux_pW flux_mJy_perbeam flux_error_mJy flux_mJy_perbeam_deboost flux_error_mJy_deboost flux_error_mJy_deboost_confusion SN_ratio_confusion flux_error_mJy_deboost_confusion_calibration flux_error_mJy_deboost_calibration")



    mock_catalogue_ds9=np.column_stack((RA, Dec))

    np.savetxt('./for_code/updated/'+cluster_list[l]+'_3p5sigma_negative.reg', mock_catalogue_ds9, delimiter=' ') #=================================================HERE





