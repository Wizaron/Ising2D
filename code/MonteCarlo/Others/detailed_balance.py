import random
import math
import numpy as np
import matplotlib.pyplot as plt
from array import *

            
def transfer_matrix():
    n=8
    #Matrix = [[0 for x in xrange(5)] for x in xrange(5)] 
    #konum = [[]*9 for x in xrange(9)]
    konum= [ [ 0 for i in range(9) ] for j in range(9) ]
    #print konum
    prob = array('f',[]) 
    print "Enter elements of transfer matrix p(a->b):"
    
    for i in range(0, 9):
        for k in range(0, 9):
            konum[i][k] = float(raw_input()) 
            
     
    print "Enter probabilities of being in positions:"
    
    for j in range(0,9):
        prob.append(float(raw_input()))
        
    
    newprob = array('f',[0 for k in range(9)]) 
    
    for a in range(0,9):
        for b in range(0,9):
            newprob[a] = newprob[a] + konum[b][a] * prob[b]
            
    for num in newprob:
        print "new probabilities=",newprob
        
        
        
transfer_matrix()
    