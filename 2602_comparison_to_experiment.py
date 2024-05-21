'''
Created on Sep 6, 2023

@author: 31417
'''
from numpy import *
from time import time
from matplotlib.pyplot import *
from pandas import *


df = read_csv("TrichoMuN2fix.csv", header=1, skiprows=0)
print(df)
mu_df = df.loc[:, 'mu']
N2_df = df.loc[:, 'N2fix']


rcParams.update({'font.size': 15})
rcParams['figure.figsize']=8,6.5
figure()
scatter(0.19834784474587802,0.08139083567183798/2,color = 'orange',label = 'H1-Model')
scatter(0.20708969862128307,0.16204195383700198/2,color = 'red', label = 'H2-Model')
scatter(mu_df,N2_df,label = 'observation')
xlabel('growth rate (d$^{-1}$)')
ylabel('${N_2}$ fixation')
legend()
savefig('figure_comparison_with_experiment.png',dpi=300)

show()
