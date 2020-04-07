#!/usr/bin/python3

import math
line=200;
column=560;
shock=28;
delta_x=0.89/line;
R=0.25/delta_x;
center_x=0.5/delta_x;
center_y=0.89/2/delta_x;
def idx(y,x):
    return math.sqrt((x-center_x)*(x-center_x)+(y-center_y)*(y-center_y)) > R

eps = 1e-4;

phi_1=1.0;
rho_2=1.0;
u_2=0.0;
p_2=1.0;
phi_2=1.0;
rho_3=rho_2*0.4*0.72/0.249/0.365;
#rho_3=rho_2*0.287/0.091;
u_3=u_2;
p_3=p_2;
phi_3=0.0;

gamma=1.4;
M=1.22;
f=1.0/(2.0/(gamma+1)/M/M+(gamma-1.0)/(gamma+1.0));
g=2*gamma/(gamma+1)*M*M-(gamma-1.0)/(gamma+1.0);
rho_1=rho_2*f;
u_1=(1.0-1.0/f)*math.sqrt(u_2+(gamma*p_2/rho_2)*M)+u_2/f;
p_1=p_2*g;

import numpy as np
from scipy.integrate import dblquad

Z_a = np.zeros((line,column));
PHI = np.zeros((line,column));
RHO = np.zeros((line,column));
U   = np.zeros((line,column));
V   = np.zeros((line,column));
P   = np.zeros((line,column));
for j in range(0,line):
    for i in range(0,shock):
       Z_a[j,i] = 1.0-eps;
       RHO[j,i] = rho_1*(1.0-eps)+rho_3*eps;
       PHI[j,i] = rho_1*(1.0-eps)/RHO[j,i];
       U[j,i]   = u_1;
       V[j,i]   = 0.0;
       P[j,i]   = p_1;      

far_p =np.array([0,0]);
near_p=np.array([0,0]);
for j in range(0,line):
    for i in range(shock,column):
        if i<center_x+0.5:
           far_p[0]  = center_x-i+1;
           near_p[0] = center_x-i; 
        elif i>center_x+0.5:
           far_p[0]  = i-center_x;
           near_p[0] = i-center_x-1;
        if j<center_y+0.5:
           far_p[1]  = center_y-j+1;
           near_p[1] = center_y-j;
        elif j>center_y+0.5:
           far_p[1]  = j-center_y;
           near_p[1] = j-center_y-1;
        if (np.linalg.norm(near_p)<R) & (np.linalg.norm(far_p)>R):
           Z_a[j,i],err = dblquad(idx,i-1.0,i+0.0,lambda y:j-1.0,lambda y:j+0.0);
           Z_a[j,i] = max(Z_a[j,i],eps);
           Z_a[j,i] = min(Z_a[j,i],1.0-eps);
        elif np.linalg.norm(near_p)>=R:
           Z_a[j,i] = 1.0-eps;
        else:
           Z_a[j,i] = eps;

for j in range(0,line):
    for i in range(shock,column):
       RHO[j,i] = Z_a[j,i]*rho_2+(1.0-Z_a[j,i])*rho_3;
       PHI[j,i] = Z_a[j,i]*rho_2/RHO[j,i];
       U[j,i]   = u_2;
       V[j,i]   = 0.0;
       P[j,i]   = p_2;      

RHO=np.fliplr(RHO);
PHI=np.fliplr(PHI);
U=np.fliplr(-1.0*U);
V=np.fliplr(V);
P=np.fliplr(P);
Z_a=np.fliplr(Z_a);

np.savetxt('PHI.dat',PHI,fmt="%g",delimiter="\t");
np.savetxt('Z_a.dat',Z_a,fmt="%g",delimiter="\t");
np.savetxt('RHO.dat',RHO,fmt="%g",delimiter="\t");
np.savetxt('U.dat',  U,  fmt="%g",delimiter="\t");
np.savetxt('V.dat',  V,  fmt="%g",delimiter="\t");
np.savetxt('P.dat',  P,  fmt="%g",delimiter="\t");
