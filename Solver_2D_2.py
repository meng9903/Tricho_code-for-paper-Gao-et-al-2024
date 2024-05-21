'''
Created on Sep 29, 2016

@author: Keisuke
'''

from pylab import *
from numpy import *

class Solver_2D:
    def __init__(self,X):
        self.X=X

def solver_2D(a,b,c):
    

    X0= -2*c/(b+sqrt(b**2-4*a*c))
    
    X=real(X0)
    X1=Solver_2D(X)
    return(X1)