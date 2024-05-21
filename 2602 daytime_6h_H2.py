'''
Created on Mar 14, 2022

@author: 31417
'''
# coding:gbk
from pylab import *
from numpy import *
from time import time
from matplotlib.pyplot import *
from Solver_2D_2 import *

dt = 1/86400
tt = arange(0,0.5+dt,dt)

a=120#1 period
b=360#1h
#define parameters related to carbon flux
pmax = 0.0035*1.65#(mol C s-1 mol chl-1) carbon fixing rate (156-10) (156-15) for unit conversion) (from a704-16)
QC = 18333 #(molC / m3) carbon density in the cell
Cmax = 36666/QC#molC/molC max carbon concentrat
Chlmaxg=0.048 #(gchl gC-1) Chlorophyll max of the cell(186-40)
Chlmaxmol=Chlmaxg*12*55/868#molCchl mol c-1
Chl=Chlmaxmol*0.8

#Chlmaxmol=Chlmaxg*12*55/868 #(molCchl molC-1) Chlorophyll max of the cell in mol (186-40)
#Chl=Chlmaxmol*0.8*QC #(molCchl m-3) Chlorophyll concentration
Ai = 0.01 #umol-1 m2 s Light saturation coefficient9o
Csto = 0.1*5 # molc molc-1 initial carbon concentration
#define parameters related to nitrogen flux
Ync = 0.159#(molN / molC) ratio of N:C in nixtrogen fixation
Nsto = 0.1# molN/molC initial nitrogen concentration
#parameters used in calculate Fbio
Kc = QC/5/QC #mol c molc-1 half saturation concentration for carbon
Kn = Kc*0.159#mol n molc-1 half saturation concentration for nitrogen
Fbio_max = 0.053*4 #mol c m-3 s-1 max biomass synthesis rate
E = 0.4 # Y res:bio
# calculate O2 concentration
O2E = 0.213 # environment O2 mol m-3
O2 = 0.213 # O2 mol m-3 initial O2 concentration
A = 3.6 # s-1 diffusion coefficient
# calculate O2 concentration
###Nitrogen fixation state
QC =18333 #(molC m-3)

KO2 = 0.00002 
Ynfix_n_o2 = 1.04 # respiration to N2 fixation ratio
Ync = 1/6.3
Mumax=0.25*4#d-1
lmax=Mumax*QC
ResMax=lmax*10 # unit (molC m-3 d-1) Maximum respiration rate 
Fres_max = ResMax*1 # (molO2 m-3 d-1)
Fnfix_max = lmax*Ync*2.5/2 # maximum nitrogenfixation rate mol N m-3 d-1
O2cr = 0.1
# set up new arrays for Carbon and Nitrogen
C = zeros(size(tt))*nan# C storage
O = zeros(size(tt))*nan# O2 concentraion mol m-3
N = zeros(size(tt))*nan# N storage
mu = zeros(size(tt))*nan# growth rate d-1
mu_m = zeros(size(tt))*nan
cfix = zeros(size(tt))*nan# c fixation rate d-1
res = zeros(size(tt))*nan#respiration rate d-1
bio = zeros(size(tt))*nan#growth d-1
cNfix = zeros(size(tt))*nan# C usage in N fixation d-1

dc_dt = zeros(size(tt))*nan# C change d-1
dn_dt = zeros(size(tt))*nan# N change mol N mol C-1 d-1
Nfix = zeros(size(tt))*nan# N fixation rate mol N mol C-1 d-1
bion = zeros(size(tt))*nan# N to growth d-1

#set the initial values
C[0] = array([Csto])
O[0] = array([0.213])
N[0] = array([Nsto])
mu[0] = array([Fbio_max/QC*86400*min(Csto/(Csto+Kc),Nsto/(Nsto+Kn))/(1+Csto)])
mu[0]
#cfix[0] = array([((Cmax - Csto)/Cmax)*pmax*Chl*86400*(1-e**(-Ai*700))])
#res[0] = array([Fbio_max*86400*min(Csto/(Csto+Kc),Nsto/(Nsto+Kn))/QC*E])
#bio[0] = array([Fbio_max*86400*min(Csto/(Csto+Kc),Nsto/(Nsto+Kn))/QC])

dc_dt[0] = array([0])
cNfix[0] = array([0])

Nfix[0] = array([0])
#bion[0] = array([(Fbio_max*86400*min(Csto/(Csto+Kc),Nsto/(Nsto+Kn))/QC)*((Ync+Nsto)/(1+Csto))])

dn_dt[0] = array([0])

cfix[0] = array([0])
res[0] = array(0)
bio[0] = array([0])
bion[0] = array([0])

for i in range(size(tt)-1):
    I = 700
    if i < 21600:#photosynthetic state
    # Calculate C and N
        
        Fcfix_max = ((Cmax - Csto)/Cmax)*pmax*Chl*86400 #max carbon fixation rate d-1
        Fcfix = Fcfix_max*(1-e**(-Ai*I)) # carbon fixation rate d-1
        Fbio = Fbio_max*86400*min(Csto/(Csto+Kc),Nsto/(Nsto+Kn))/QC # biomass synthesis rate d-1
        Fnfix_c = Fnfix_max*(Csto/(Csto+Kc))*1# C usage in N fixation mol C m-3 d-1
        #Fo2 = A*86400*(O2E-O2)# O2 diffusivity molO2 m-3 d-1
        Fnfix = Fnfix_c*1/QC# N2 fixation rate mol N mol C-1 d-1
        dcdt = Fcfix - Fbio*(1+E) - Fnfix*1 - Fnfix*Ynfix_n_o2*1# carbon changing rate d-1
        dndt = Fnfix -Fbio*((Ync+Nsto)/(1+Csto))# nitrogen changing rate mol N mol C-1 d-1
        #dodt = A*86400*(O2E-O2)+Fcfix*QC*1-Fbio*E*QC*1-Fnfix_c*Ynfix_n_o2*1*1*1
        Csto = Csto+dcdt*dt # carbon concentration in cell molC molC-1
        Nsto = Nsto+dndt*dt # nitrogen concentration in cell molN molC-1
        O2 = (A*86400*O2E+Fcfix*QC*1-Fbio*E*QC*1-Fnfix_c*Ynfix_n_o2*1)/(A*86400)# O2 calculation mol m-3
        #O2 = O2+dodt*dt
        C[i+1] = Csto
        N[i+1] = Nsto
        O[i+1] = O2
        mu[i+1] = Fbio/(1+Csto)
        cfix[i+1] = Fcfix
        res[i+1] = Fbio*E + Fnfix*Ynfix_n_o2*1
        bio[i+1] = Fbio
        cNfix[i+1] = Fnfix*1
        dc_dt[i+1] = dcdt
        dn_dt[i+1] = dndt
        Nfix[i+1] = Fnfix
        bion[i+1] = Fbio*((Ync+Nsto)/(1+Csto))
        
    elif i < 43200:#non photosynthetic state
        Csto = C[i]
        Nsto = N[i]
        X1=solver_2D(A*86400,-A*86400*O2E-A*86400*KO2+Fres_max,-A*86400*O2E*KO2)
        X = X1.X
        O2 = X# O2 calculation mol m-3
        Fnfix_c = Fnfix_max*(Csto/(Csto+Kc))*1# C usage in N fixation mol C m-3 d-1
        Fres = Fres_max*(O2/(O2+KO2))# respiration rate (molO2 m-3 d-1)
        Fo2 = A*86400*(O2E-O2)# O2 diffusivity molO2 m-3 d-1
        Fnfix = Fnfix_c*1/QC# N2 fixation rate mol N mol C-1 d-1
        Fbio = 0# biomass synthesis rate molC d-1
        dcdt = -Fnfix*1 - Fres/(1*QC)-Fbio# carbon concentration in cell molC molC-1 d-1
        dndt = Fnfix-Fbio*((Ync+Nsto)/(1+Csto))# nitrogen concentration in cell molN molC-1 d-1
        #dodt = A*86400*(O2E-O2)+0*QC*1-Fbio*E*QC*1-Fnfix_c*Ynfix_n_o2
        #dodt = Fo2-Fres
        Csto=Csto+dcdt*dt
        Nsto=Nsto+dndt*dt
        #O2 = O2 + dodt*dt
        C[i+1] = Csto# carbon concentration in cell molC molC-1
        N[i+1] = Nsto# nitrogen concentration in cell molN molC-1
        O[i+1] = O2
        mu[i+1] = Fbio/(1+Csto)
        cfix[i+1] = 0
        res[i+1] = Fres/(1*QC)
        bio[i+1] = Fbio
        cNfix[i+1] = Fnfix*1
        dc_dt[i+1] = dcdt
        dn_dt[i+1] = dndt
        Nfix[i+1] = Fnfix
        bion[i+1] = Fbio*((Ync+Nsto)/(1+Csto))
 

#Calculate total C and N usage
CSTO = copy(dc_dt)
TC = bio+cNfix+res+CSTO

############
NSTO = copy(dn_dt)
TN = bion+NSTO

############
############ print each pathway usage precentage
print("res",sum(res)/sum(TC))
print("bio",sum(bio)/sum(TC))
print("cNifix" ,sum(cNfix)/sum(TC))
print("Csto" ,sum(CSTO)/sum(TC))

print("")

print("Nsto",sum(NSTO)/sum(TN))
print("Nbio",sum(bion)/sum(TN))


# make arrays to summarize C and N allocation
SC=array([["","C fate","C fate fraction"
],
              ["Csto1",sum(CSTO)*dt,sum(CSTO)/sum(TC)],
              ["bio",sum(bio)*dt,sum(bio)/sum(TC)],
              ["cNfix",sum(cNfix)*dt,sum(cNfix)/sum(TC)],
              ["res",sum(res)*dt,sum(res)/sum(TC)],
         ])
SN=array([["","N fate","N fate fraction"
],
              ["Nsto1",sum(NSTO)*dt,sum(NSTO)/sum(TN)],
              ["bio",sum(bion)*dt,sum(bion)/sum(TN)],
         ])


print(mean(mu)) ## calculate average growth rate

figure('mu')           
plot(tt,mu)
xlabel('t (day)')
ylabel('growth rate')

rcParams.update({'font.size': 15})
figure('O2')           
plot(tt*24,O)
xlabel('hour')
ylabel('${O_2}$(mol m$^{-3}$)')


show()

savetxt('mu_6h_H2_0207.csv', mu, delimiter=",", fmt='%.8e')
savetxt('SC6_6h_H2_0522.csv', SC, delimiter="," ,fmt ='%s')
savetxt('SN6_6h_H2_0522.csv', SN, delimiter=",", fmt="%s")
savetxt('Nfix_6h_H2.csv', Nfix, delimiter=",", fmt="%s")
