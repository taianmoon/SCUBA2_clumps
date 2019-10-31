import numpy as np
import corner
import pyfits
import math
import matplotlib.pyplot as plt




#z_want_arp220 = np.linspace(0.0, 6.0, num=50)  #====================================================Redshift for the SED


sed_names=['cosmic_eyelash', 'HFLS3', 'Arp220', 'ALESS']
#sed_names=['Arp220', 'ALESS']

for i in range(0, len(sed_names)):
    print sed_names[i]

    template_file = open("./"+str(sed_names[i])+"_LIR.cat", "r") 

    lines = template_file.readlines()
    template_file.close()
    
    response_curve=[]

    for i in range(0, len(lines)):
        separated_lines=lines[i].split() 
        response_curve.append(separated_lines)


    response_curve = np.array(response_curve)
    z_want_arp220=response_curve[:,0]
    L_IR=response_curve[:,1]

    #L_IR = np.array(L_IR)
    #L_IR=np.array([float(i) for i in L_IR])  #(wavelength in angstrom)

    print len(z_want_arp220)
    print len(L_IR)
    print 'dddddddddddddddddddd'
    plt.plot(z_want_arp220, L_IR, linewidth=1.5, color='b')



plt.scatter([2.86, 3.35], [12.85, 13.09], color='r', s=20.0)
plt.errorbar([2.86, 3.35], [12.85, 13.09], xerr=[0.96, 1.09], yerr=[0.22, 0.23], color='r')
plt.xlabel('z')
plt.ylabel('log10(L_IR)')
plt.legend()
plt.grid()
plt.rc('font', size=20)
plt.show()
