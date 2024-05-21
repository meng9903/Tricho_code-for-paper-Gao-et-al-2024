'''
Created on Mar 21, 2022

@author: 31417
'''
from pandas import*
from pylab import *
from numpy import *
from matplotlib.pyplot import *

rcParams.update({'font.size': 15})

m1 = read_csv('mu_1min_H2_0207.csv',header=None)
m2 = read_csv('mu_6h_H2_0207.csv',header=None)
m1_1 = read_csv('mu_1min_H1_0524.csv',header=None)
m2_1 = read_csv('mu_6h_H1_0207.csv',header=None)
#nc1_1 = read_csv('n_c_1min_H1_0207.csv',header=None)
#nc1_2 = read_csv('n_c_6h_H1_0127.csv',header=None)
df = read_csv('SCT_0207.csv')
dfn = read_csv('SNT_0207.csv')
dt = 1/86400
tt = arange(0,0.5+dt,dt)
print(df)
print(dfn)
print(size(tt))

## figure 2c
figure('growth rate Comparison')
plot(tt*24,m1_1,label = 'H1 1min:1min',color='red')
plot(tt*24,m2_1,label = 'H1 6h:6h',color='red',linestyle = 'dashed')
xlabel('hour')
ylabel('growth rate (d$^{-1}$)')
xlim(0,12);xticks(arange(0,12+1e-6,2))
ylim(0,0.5);yticks(arange(0,0.5,0.1))
title('growth rate comparison')
legend()
savefig('growth rate comparison0614_H1.png',dpi=300)

## figure 2d
figure()
plot(tt*24,m1,label = 'H2 1min:1min',color='blue')
plot(tt*24,m2,label = 'H2 6h:6h',color='blue',linestyle = 'dashed')
xlabel('hour')
ylabel('growth rate (d$^{-1}$)')
xlim(0,12+1e-6);xticks(arange(0,12+1e-6,2))
ylim(0,0.5);yticks(arange(0,0.5,0.1))
title('growth rate comparison')
legend()
savefig('growth rate comparison0614_H2.png',dpi=300)


rcParams.update({'font.size': 20})
rcParams['axes.prop_cycle'] = cycler(color=['c','orange','green','pink']) 

### figure 3
df1 = df.loc[[4,6],:]
figure('C fate')
df1.plot(x='time to switch',kind='barh', stacked=True, mark_right=True,figsize=(10,8))
axvline(0, color='black')
xlabel('C (molC molC$^{-1}$) used for 12 hours')
ylabel('')
title('C fate')
savefig('C fate0614_H2.png',dpi=300)

df1 = df.loc[[0,2],:]
figure('C fate')
df1.plot(x='time to switch',kind='barh', stacked=True, mark_right=True, figsize=(10,8))
xlabel('C (molC molC$^{-1}$) used for 12 hours')
ylabel('')
title('C fate')
savefig('C fate0614_H1.png',dpi=300)

dfn1 = dfn.loc[[4,6],:]
dfn1.plot(x='time to switch', kind='barh', stacked=True, mark_right=True,figsize=(10,8))
xlabel('N (molN molC$^{-1}$) used for 12 hours')
ylabel('')
title('N fate')
legend()
savefig('N fate0614_H2.png',dpi=300)

dfn1 = dfn.loc[[0,2],:]
dfn2 = dfn.loc[[1,3],:]
dfn1.plot(x='time to switch', kind='barh', stacked=True, mark_right=True,figsize=(10,8))
xlabel('N (molN molC$^{-1}$) used for 12 hours')
ylabel('')
title('N fate')
legend()
savefig('N fate0614_H1.png',dpi=300)

show()
