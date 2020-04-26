#!/usr/bin/python3

line=400;
column=1120;
shock=56;
delta_x=0.89/line;
R=0.25/delta_x;
center_x=0.5/delta_x-1;
center_y=0.89/2/delta_x-1;

Z_a_in  = 0.95;
Z_a_out = 0.05;
gamma_1 = 1.4;
gamma_2 = 1.67;

import math
def idx(y,x):
    return math.sqrt((x-center_x)**2+(y-center_y)**2) > R;

rho_a_2=1.0;
u_a_2=0.0;
p_a_2=1.0;
rho_b=0.1821;
u_b=u_a_2;
p_b=p_a_2;

M=1.22;
f=1.0/(2.0/(gamma_1+1.0)/M/M+(gamma_1-1.0)/(gamma_1+1.0));
g=2*gamma_1/(gamma_1+1.0)*M*M-(gamma_1-1.0)/(gamma_1+1.0);
rho_a_1=rho_a_2*f;
u_a_1=(1.0-1.0/f)*(u_a_2+math.sqrt(gamma_1*p_a_2/rho_a_2)*M)+u_a_2/f;
p_a_1=p_a_2*g;

import numpy as np
from scipy.integrate import dblquad

Z_a   = np.zeros((line,column));
RHO_a = np.zeros((line,column));
U_a   = np.zeros((line,column));
V_a   = np.zeros((line,column));
P_a   = np.zeros((line,column));
RHO_b = np.zeros((line,column));
U_b   = np.zeros((line,column));
V_b   = np.zeros((line,column));
P_b   = np.zeros((line,column));
err   = np.zeros((line,column));

for j in range(0,line):
    for i in range(0,shock):
       Z_a[j,i]   = Z_a_out;
       RHO_a[j,i] = rho_a_1;
       U_a[j,i]   = u_a_1;
       V_a[j,i]   = 0.0;
       V_b[j,i]   = 0.0;
       P_a[j,i]   = p_a_1;

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
           Z_a[j,i],err[j,i] = dblquad(idx,i-1.0,i-0.0,lambda y:j-1.0,lambda y:j-0.0);
           Z_a[j,i] = Z_a[j,i]*Z_a_out+(1-Z_a[j,i])*Z_a_in;
        elif np.linalg.norm(near_p)>=R:
           Z_a[j,i] = Z_a_out;
        else:
           Z_a[j,i] = Z_a_in;

for j in range(0,line):
    for i in range(shock,column):
       RHO_a[j,i] = rho_a_2;
       U_a[j,i]   = u_a_2;
       V_a[j,i]   = 0.0;
       P_a[j,i]   = p_a_2;     

for j in range(0,line):
    for i in range(0,column):
       RHO_b[j,i] = rho_b;
       U_b[j,i]   = u_b;
       V_b[j,i]   = 0.0;
       P_b[j,i]   = p_b;
       
Z_a   = np.fliplr(Z_a);
RHO_a = np.fliplr(RHO_a);
U_a   = np.fliplr(-1.0*U_a);
V_a   = np.fliplr(V_a);
P_a   = np.fliplr(P_a);
RHO_b = np.fliplr(RHO_b);
U_b   = np.fliplr(-1.0*U_b);
V_b   = np.fliplr(V_b);
P_b   = np.fliplr(P_b);

np.savetxt('Z_a.dat',  1.0-Z_a,  fmt="%.10g",delimiter="\t");
np.savetxt('RHO.dat',  RHO_a,fmt="%.10g",delimiter="\t");
np.savetxt('U.dat',    U_a,  fmt="%.10g",delimiter="\t");
np.savetxt('V.dat',    V_a,  fmt="%.10g",delimiter="\t");
np.savetxt('P.dat',    P_a,  fmt="%.10g",delimiter="\t");
np.savetxt('RHO_b.dat',RHO_b,fmt="%.10g",delimiter="\t");
np.savetxt('U_b.dat',  U_b,  fmt="%.10g",delimiter="\t");
np.savetxt('V_b.dat',  V_b,  fmt="%.10g",delimiter="\t");
np.savetxt('P_b.dat',  P_b,  fmt="%.10g",delimiter="\t");
