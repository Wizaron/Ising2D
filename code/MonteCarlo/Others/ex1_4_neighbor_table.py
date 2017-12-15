import random
import math
import numpy as np
import matplotlib.pyplot as plt
from array import *
import pylab as pl

def nbr(n,k):
    if k == 1: 
        if n == 1:
            return 2
        elif n == 2:
            return 4
        else: 
            return 0
    if k == 2: 
        if n == 1:
            return 3
        elif n == 2:
            return 5
        elif n == 3: 
            return 1
        else: 
            return 0
    if k == 3: 
        if n == 2:
            return 6
        elif n == 3:
            return 2
        else: 
            return 0
    if k == 4: 
        if n == 1:
            return 5
        elif n == 2:
            return 7
        elif n == 4:
            return 1
        else: 
            return 0
    if k == 5: 
        if n == 1:
            return 6
        elif n == 2:
            return 8
        elif n == 3: 
            return 4
        elif n == 4:
            return 2
    if k == 6: 
        if n == 2:
            return 9
        elif n == 3:
            return 5
        elif n == 4:
            return 3
        else:
            return 0
    if k == 7:
        if n == 1:
            return 8
        elif n == 4:
            return 4
        else: 
            return 0
    if k == 8:
        if n == 1:
            return 9
        elif n == 3:
            return 7
        elif n == 4:
            return 5
        else: 
            return 0
    if k == 9: 
        if n == 3:
            return 8
        elif n == 4:
            return 6
        else: 
            return 0
            
def markov_discrete_pebble():
    n=1000
    konum = array('i',[0 for k in range(9)] )
    print "Enter the position of pebble:"
    position = int(raw_input()) 
    print "First position=" ,position
    konum[position-1] += 1

    for i in range(1, n):
        move = random.randrange(1, 5)
        #print "random=",move
        if nbr(move,position) != 0:
            position = nbr(move,position)    
        
        print "New Position=", position
        konum[position-1] += 1
    
    print "konum=",konum
    
    #x=array('i',[1,2,3,4,5,6,7,8,9])
    #plt.plot(yer,konum,'ro')
    
    #x=range(1,10)
    x = np.linspace(1,9,9)
    pl.clf()
    pl.plot(x,konum,'ro')
    
    p = np.polyfit(x,konum,0)
    #print " p=" ,p
    pl.plot( x, p[0] + 0*x, 'r')
    
    plt.ylabel('Bulunulma sayisi')
    plt.xlabel('Konum')
    plt.title('Bir konumda bulunma sayisi vs konum')  
    plt.show()
    
markov_discrete_pebble()