import matplotlib.pyplot as plt
import numpy as np
import math




cluster_list=['Planck18p194', 'Planck18p735', 'Planck24p194', 'PLCK_DU_G045.7-41.2', 'PLCK_DU_G059.1-67.1', 'PLCK_DU_G073.4-57.5', 'PLCK_G006.1+61.8', 'PLCK_G009.8+72.6', 'PLCK_G056.7+62.6', 'PLCK_G068.3+31.9', 'PLCK_G075.1+33.2', 'PLCK_G077.7+32.6', 'PLCK_G078.9+48.2', 'PLCK_G082.5+38.4', 'PLCK_G083.3+51.0', 'PLCK_G091.9+43.0', 'PLCK_G093.6+55.9', 'PLCK_G132.9-76.0', 'PLCK_G144.1+81.0', 'PLCK_G160.7+41.0', 'PLCK_G162.1-59.3', 'PLCK_G165.8+45.3', 'PLCK_G173.8+59.3', 'PLCK_G177.0+35.9', 'PLCK_G179.3+50.7', 'PLCK_G186.3-72.7', 'PLCK_G186.6+66.7', 'PLCK_G188.6-68.9', 'PLCK_G191.3+62.0', 'PLCK_G191.8-83.4', 'PLCK_G201.1+50.7', 'PLCK_G213.0+65.9', 'PLCK_G223.9+41.2', 'PLCK_G328.9+71.4', 'PLCK_G49.6-42.9', 'PLCK_G84.0-71.5', 'PLCK_HZ_G038.0-51.5', 'PLCK_HZ_G067.2-63.8', 'PLCK_HZ_G103.1-73.6', 'PLCK_HZ_G106.8-83.3', 'PLCK_HZ_G119.4-76.6', 'PLCK_HZ_G132.6-81.1', 'PLCK_HZ_G171.1-78.7', 'PLCK_HZ_G173.9+57.0', 'PLCK_HZ_G176.6+59.0', 'PLCK_HZ_G214.1+48.3']




for l in range(0, len(cluster_list)):

    print cluster_list[l]



    template_file = open("./catalogues/edgesourcedelete/"+cluster_list[l]+"_3p5sigma.cat", "r")

    lines = template_file.readlines()[1:]
    template_file.close()

    response_curve=[]

    for i in range(0, len(lines)):
        separated_lines=lines[i].split() 
        response_curve.append(separated_lines)


    response_curve = np.array(response_curve)
    index_name=response_curve[:,0]
    pix_coor_x_3sigma_sorted_extract=response_curve[:,1]
    pix_coor_y_3sigma_sorted_extract=response_curve[:,2]
    ra_3sigma_sorted_extract=response_curve[:,3]
    dec_3sigma_sorted_extract=response_curve[:,4]
    pix_value_array_3sigma_sorted_extract=response_curve[:,5]
    pix_value_flux_array_3sigma_sorted_extract=response_curve[:,6]
    pix_value_flux_array_3sigma_sorted_extract_mJy=response_curve[:,7]
    pix_value_flux_error_array_3igma_sorted_extract_mJy=response_curve[:,8]



    index_name = np.array(index_name)
    index_name=np.array([float(i) for i in index_name])  #(wavelength in angstrom)
    pix_coor_x_3sigma_sorted_extract = np.array(pix_coor_x_3sigma_sorted_extract)
    pix_coor_x_3sigma_sorted_extract=np.array([float(i) for i in pix_coor_x_3sigma_sorted_extract])  #(wavelength in angstrom)
    pix_coor_y_3sigma_sorted_extract = np.array(pix_coor_y_3sigma_sorted_extract)
    pix_coor_y_3sigma_sorted_extract=np.array([float(i) for i in pix_coor_y_3sigma_sorted_extract])  #(wavelength in angstrom)
    ra_3sigma_sorted_extract = np.array(ra_3sigma_sorted_extract)
    ra_3sigma_sorted_extract=np.array([float(i) for i in ra_3sigma_sorted_extract])  #(wavelength in angstrom)
    dec_3sigma_sorted_extract = np.array(dec_3sigma_sorted_extract)
    dec_3sigma_sorted_extract=np.array([float(i) for i in dec_3sigma_sorted_extract])  #(wavelength in angstrom)
    pix_value_array_3sigma_sorted_extract = np.array(pix_value_array_3sigma_sorted_extract)
    pix_value_array_3sigma_sorted_extract=np.array([float(i) for i in pix_value_array_3sigma_sorted_extract])  #(wavelength in angstrom)
    pix_value_flux_array_3sigma_sorted_extract = np.array(pix_value_flux_array_3sigma_sorted_extract)
    pix_value_flux_array_3sigma_sorted_extract=np.array([float(i) for i in pix_value_flux_array_3sigma_sorted_extract])  #(wavelength in angstrom)
    pix_value_flux_array_3sigma_sorted_extract_mJy = np.array(pix_value_flux_array_3sigma_sorted_extract_mJy)
    pix_value_flux_array_3sigma_sorted_extract_mJy=np.array([float(i) for i in pix_value_flux_array_3sigma_sorted_extract_mJy])  #(wavelength in angstrom)
    pix_value_flux_error_array_3igma_sorted_extract_mJy = np.array(pix_value_flux_error_array_3igma_sorted_extract_mJy)
    pix_value_flux_error_array_3igma_sorted_extract_mJy=np.array([float(i) for i in pix_value_flux_error_array_3igma_sorted_extract_mJy])  #(wavelength in angstrom)


    flux_mjy_deboost=[]
    flux_mjy_deboost_err=[]

    for k in range(0, len(pix_value_flux_array_3sigma_sorted_extract_mJy)):
        flux_mjy_deboost.append(pix_value_flux_array_3sigma_sorted_extract_mJy[k]/(1.0+0.2*pow(pix_value_array_3sigma_sorted_extract[k]/5.0,-2.3)))
        flux_mjy_deboost_err.append(pix_value_flux_error_array_3igma_sorted_extract_mJy[k]/(1.0+0.2*pow(pix_value_array_3sigma_sorted_extract[k]/5.0,-2.3)))


    flux_mjy_deboost = np.array(flux_mjy_deboost)
    flux_mjy_deboost=np.array([float(i) for i in flux_mjy_deboost])  #(wavelength in angstrom)
    flux_mjy_deboost_err = np.array(flux_mjy_deboost_err)
    flux_mjy_deboost_err=np.array([float(i) for i in flux_mjy_deboost_err])  #(wavelength in angstrom)




    mock_catalogue=np.column_stack((index_name, pix_coor_x_3sigma_sorted_extract, pix_coor_y_3sigma_sorted_extract, ra_3sigma_sorted_extract, dec_3sigma_sorted_extract, pix_value_array_3sigma_sorted_extract, pix_value_flux_array_3sigma_sorted_extract, pix_value_flux_array_3sigma_sorted_extract_mJy, pix_value_flux_error_array_3igma_sorted_extract_mJy, flux_mjy_deboost, flux_mjy_deboost_err))

    np.savetxt('./catalogues/edgesourcedelete_deboost_181008/'+cluster_list[l]+'_3p5sigma_deboost.cat', mock_catalogue, delimiter=' ', header="name pix_x pix_y RA Dec SN_ratio flux_pW flux_mJy_perbeam flux_error_mJy flux_mJy_perbeam_deboost flux_error_mJy_deboost")


    mock_catalogue_ds9=np.column_stack((ra_3sigma_sorted_extract, dec_3sigma_sorted_extract))

    np.savetxt('./catalogues/edgesourcedelete_deboost_181008/'+cluster_list[l]+'_3p5sigma_deboost.reg', mock_catalogue_ds9, delimiter=' ') #================================================================================HERE


