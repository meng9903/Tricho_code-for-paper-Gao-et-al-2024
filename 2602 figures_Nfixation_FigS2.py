'''
Created on Mar 21, 2022

@author: 31417
'''
from pandas import*
from pylab import *
from numpy import *
from matplotlib.pyplot import *
## set figures size
## This file is for Figure S2
rcParams.update({'font.size': 18})
rcParams['figure.figsize']=8,6.5

## read data
n1 = read_csv('Nfix_1min_H1.csv',header=None)
n2 = read_csv('Nfix_1min_H2.csv',header=None)
n1_1 = read_csv('Nfix_6h_H1.csv',header=None)
n2_1 = read_csv('Nfix_6h_H2.csv',header=None)

dt = 1/86400
tt = arange(0,0.5+dt,dt)

# Nitrogen Fixation H1 Comparison Figure for 12 hours Figure
figure(1)
plot(tt*24,n1,label = 'H1 1min:1min',color='red')
plot(tt*24,n1_1,label = 'H1 6h:6h',color='black',linestyle = 'dashed')
xlabel('hour')
ylabel('${N_2}$ fixation rate (molN molC$^{-1}$ d$^{-1}$)')
xlim(0,12);xticks(arange(0,12+1e-6,2))
ylim(0,0.3);yticks(arange(0,0.3,0.1))
title('${N_2}$ fixation')
legend()
savefig('N2 fixation comparison0614_H1_T.png',dpi=300)
print(1)

# Nitrogen Fixation H2 Comparison Figure for 12 hours
figure(2)
plot(tt*24,n2,label = 'H2 1min:1min',color='blue')
plot(tt*24,n2_1,label = 'H2 6h:6h',color='black',linestyle = 'dashed')
xlabel('hour')
ylabel('${N_2}$ fixation rate (molN molC$^{-1}$ d$^{-1}$)')
xlim(0,12);xticks(arange(0,12+1e-6,2))
ylim(0,0.3);yticks(arange(0,0.3,0.1))
title('${N_2}$ fixation')
legend()
savefig('N2 fixation comparison0614_H2_T.png',dpi=300)
print(2)


show()

