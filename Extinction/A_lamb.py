'''A module for computing the A_lambda given teh Ia_A_poly.pickle
file'''

import pickle
from numpy import power

f = open('Ia_A_poly.pickle')
d = pickle.load(f)
f.close()

def A_lamb(f, EBV, Rv, redlaw='ccm'):
   global d
   
   order = d[redlaw]['order'][f]
   coefs = d[redlaw][f]
   Al = EBV*0
   id = 0
   for j in range(order+1):
      for i in range(order+1):
         if i + j <= order:
            Al = Al + coefs[id]*power(Rv,i)*power(EBV,j+1)
            id += 1
   return Al

