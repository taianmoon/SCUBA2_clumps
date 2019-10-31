import numpy as np
import matplotlib.pyplot as plt


template_file = open("./all_source_3p5sigma_edgesourcedelete_deboost_match_herschel_separation.cat", "r")    #=====================================================HERE


lines = template_file.readlines()[1:]
template_file.close()

response_curve=[]

for i in range(0, len(lines)):
    separated_lines=lines[i].split() 
    response_curve.append(separated_lines)


response_curve = np.array(response_curve)

separation=response_curve[:,0]



separation = np.array(separation)
separation=np.array([float(i) for i in separation])  #(wavelength in angstrom)


#print separation



template_file_random = open("./all_source_3p5sigma_edgesourcedelete_deboost_random_separation.cat", "r")    #=============================================================HERE


lines_random = template_file_random.readlines()[1:]
template_file_random.close()

response_curve_random=[]

for i in range(0, len(lines_random)):
    separated_lines_random=lines_random[i].split() 
    response_curve_random.append(separated_lines_random)


response_curve_random = np.array(response_curve_random)

separation_random=response_curve_random[:,0]



separation_random = np.array(separation_random)
separation_random=np.array([float(i) for i in separation_random])  #(wavelength in angstrom)


#print separation_random



plt.hist(separation, bins=np.arange(1.0, 26.0, 1.0), normed=True, alpha = 0.5, color='green', label='original sources')
plt.hist(separation_random, bins=np.arange(1.0, 26.0, 1.0), normed=True, alpha = 0.5, color='blue', label='random offset')
plt.xlabel('Separation (arcminute)')
plt.ylabel('Normalized counts')
plt.legend()
#plt.grid()
plt.rc('font', size=18)
plt.tick_params(width=2, length=16, which='major')
plt.tick_params(width=2, length=5, which='minor')

plt.show()




