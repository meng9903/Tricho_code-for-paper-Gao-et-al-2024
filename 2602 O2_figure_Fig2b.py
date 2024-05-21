'''
Created on Mar 1, 2022

@author: 31417
'''
from pandas import*
from pylab import *
from numpy import *
from matplotlib.pyplot import *

rcParams.update({'font.size': 15})
O2 = read_csv('O2_0524.csv',header=None)
O1 = read_csv('AO3_0524.csv',header=None)

dt = 0.000001/86400
tt = arange(0,120000001*dt,dt)
dt1 =1/86400
tt1 = arange(0,121*dt1,dt1)
print(O2)

figure('O2 Comparison')
plot(tt*86400,O2,label = 'real-state')
plot(tt1*86400,O1,label = 'pseudo-state')
xlabel('time (s)')
ylabel('${O_2}$(mol m$^{-3}$)')
legend()
savefig('figure_O2_0524.png',dpi=300)
show()


