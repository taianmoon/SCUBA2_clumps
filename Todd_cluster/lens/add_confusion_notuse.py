import numpy as np
import math








template_file = open("./Planck18p194/Planck18p194_3p5sigma_deboost.cat", "r")    #========================================================HERE

lines = template_file.readlines()[1:]
template_file.close()

response_curve=[]

for i in range(0, len(lines)):
    separated_lines=lines[i].split() 
    response_curve.append(separated_lines)


response_curve = np.array(response_curve)
name=response_curve[:,0]
pix_x=response_curve[:,1]
pix_y=response_curve[:,2]
RA=response_curve[:,3]
Dec=response_curve[:,4]
SN_ratio=response_curve[:,5]
flux_pW=response_curve[:,6]
flux_mJy_perbeam=response_curve[:,7]
flux_error_mJy=response_curve[:,8]
flux_mJy_perbeam_deboost=response_curve[:,9]
flux_error_mJy_deboost=response_curve[:,10]
flux_error_mJy_deboost_confusion=response_curve[:,11]
SN_ratio_confusion=response_curve[:,12]
flux_error_mJy_deboost_confusion_calibration=response_curve[:,13]
Herschel_250_mJy=response_curve[:,14]
Herschel_250err_mJy=response_curve[:,15]
Herschel_350_mJy=response_curve[:,16]
Herschel_350err_mJy=response_curve[:,17]
Herschel_500_mJy=response_curve[:,18]
Herschel_500err_mJy=response_curve[:,19]



name = np.array(name)
name=np.array([float(i) for i in name])
pix_x = np.array(pix_x)
pix_x=np.array([float(i) for i in pix_x])
pix_y = np.array(pix_y)
pix_y=np.array([float(i) for i in pix_y])
RA = np.array(RA)
RA=np.array([float(i) for i in RA])
Dec = np.array(Dec)
Dec=np.array([float(i) for i in Dec])
SN_ratio = np.array(SN_ratio)
SN_ratio=np.array([float(i) for i in SN_ratio])
flux_pW = np.array(flux_pW)
flux_pW=np.array([float(i) for i in flux_pW])
flux_mJy_perbeam = np.array(flux_mJy_perbeam)
flux_mJy_perbeam=np.array([float(i) for i in flux_mJy_perbeam])
flux_error_mJy = np.array(flux_error_mJy)
flux_error_mJy=np.array([float(i) for i in flux_error_mJy])
flux_mJy_perbeam_deboost = np.array(flux_mJy_perbeam_deboost)
flux_mJy_perbeam_deboost=np.array([float(i) for i in flux_mJy_perbeam_deboost])
flux_error_mJy_deboost = np.array(flux_error_mJy_deboost)
flux_error_mJy_deboost=np.array([float(i) for i in flux_error_mJy_deboost])
flux_error_mJy_deboost_confusion = np.array(flux_error_mJy_deboost_confusion)
flux_error_mJy_deboost_confusion=np.array([float(i) for i in flux_error_mJy_deboost_confusion])
SN_ratio_confusion = np.array(SN_ratio_confusion)
SN_ratio_confusion=np.array([float(i) for i in SN_ratio_confusion])
flux_error_mJy_deboost_confusion_calibration = np.array(flux_error_mJy_deboost_confusion_calibration)
flux_error_mJy_deboost_confusion_calibration=np.array([float(i) for i in flux_error_mJy_deboost_confusion_calibration])
Herschel_250_mJy = np.array(Herschel_250_mJy)
Herschel_250_mJy=np.array([float(i) for i in Herschel_250_mJy])
Herschel_250err_mJy = np.array(Herschel_250err_mJy)
Herschel_250err_mJy=np.array([float(i) for i in Herschel_250err_mJy])
Herschel_350_mJy = np.array(Herschel_350_mJy)
Herschel_350_mJy=np.array([float(i) for i in Herschel_350_mJy])
Herschel_350err_mJy = np.array(Herschel_350err_mJy)
Herschel_350err_mJy=np.array([float(i) for i in Herschel_350err_mJy])
Herschel_500_mJy = np.array(Herschel_500_mJy)
Herschel_500_mJy=np.array([float(i) for i in Herschel_500_mJy])
Herschel_500err_mJy = np.array(Herschel_500err_mJy)
Herschel_500err_mJy=np.array([float(i) for i in Herschel_500err_mJy])


print Herschel_250err_mJy
print Herschel_350err_mJy
print Herschel_500err_mJy

Herschel_250err_mJy_confusion=[]
Herschel_350err_mJy_confusion=[]
Herschel_500err_mJy_confusion=[]

for i in range(0, len(name)):
    Herschel_250err_mJy_confusion.append(math.sqrt(pow(Herschel_250err_mJy[i], 2.0)+pow(5.8, 2.0)))
    Herschel_350err_mJy_confusion.append(math.sqrt(pow(Herschel_350err_mJy[i], 2.0)+pow(6.3, 2.0)))
    Herschel_500err_mJy_confusion.append(math.sqrt(pow(Herschel_500err_mJy[i], 2.0)+pow(6.8, 2.0)))


print Herschel_250err_mJy_confusion
print Herschel_350err_mJy_confusion
print Herschel_500err_mJy_confusion








