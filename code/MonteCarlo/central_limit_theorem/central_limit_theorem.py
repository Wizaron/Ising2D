# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import random
import math
import numpy as np
import matplotlib.pyplot as plt
from array import *
import matplotlib.mlab as mlab

mean = 2 # for uniform distribution in between [0,4]  E[X]=2
var = 2 # for uniform distribution in between [0,4] Var[X]=2

# <codecell>

def clt(k,n):
    a=[]
    for i in xrange(0,4*n+1):
        a.append(0);
    
    for j in range (1,k+1):
        w_n=0
        for i in range (1,n+1):
            w_n += random.randrange(0,5)
        a[w_n] += 1.0*math.sqrt(n*var)/k
    
    return a
       
    

# <codecell>

k=10000
for n in range(1,5):
  x=(np.linspace(0,4*n,4*n+1)-n*mean)/math.sqrt(n*var)
  y=clt(k,n)
  plt.subplot(2,2,n)
  plt.plot(x,mlab.normpdf(x,0,1))
  plt.plot(x,y,'ro')

# <codecell>

k=100000
for n in range(5,9):
   x=(np.linspace(0,4*n,4*n+1)-n*mean)/math.sqrt(n*var)
  #x=linspace(0,4*n,4*n+1)
   #print x
   y=clt(k,n)
   plt.subplot(2,2,n-4)
   plt.plot(x,mlab.normpdf(x,0,1))
   plt.plot(x,y,'ro')
   

# <codecell>


