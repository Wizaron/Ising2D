import random
import math
import numpy as np
import matplotlib.pyplot as plt
from array import *

            
def two_side_metropolis():
    print "Enter the position of pebble: 0 or 1"
    k = int(raw_input()) 
    if k == 0:
        l=1
    if k == 1:
        l=0
    print "Enter prob. of being at initial position:"
    probk = float(raw_input())
    
    print "Enter prob. of being at other position:"
    probl = float(raw_input())
    
    gamma = probl/probk
    rand = random.uniform(0, 1) 
    
    if gamma > rand:
        k = l 
        
    print "final position:" ,k

two_side_metropolis()