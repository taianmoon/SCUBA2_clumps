import math
import pyfits
import numpy as np
from decimal import Decimal







names_wht = pyfits.open('./ncounts/ncounts.fits')                            #===========================================HERE

names_wht_data = names_wht[1].data


observe_n= names_wht_data.field('col5')    #========================================================HERE(Boote.s: F250, et_F250; EG.S/Lockma.n: f250, et250; G1.2/NG.P: F250, E250)
predict_n= names_wht_data.field('col7')    #========================================================HER



#observe_n = np.array(observe_n)
#observe_n=np.array([float(i) for i in observe_n])
#predict_n = np.array(predict_n)
#predict_n=np.array([float(i) for i in predict_n])

prob=[]
for i in range(0, len(observe_n)):

    prob_loop=math.exp(-1.0*predict_n[i])*pow(predict_n[i], observe_n[i])/math.factorial(int(observe_n[i]))
    print prob_loop
    prob.append(prob_loop)
    
