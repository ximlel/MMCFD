#!/usr/bin/python3

line=160;
column=600;
shock=590;
X_C=column/2;
Y_C=line/2;
d_y=0.08/line;
d_x=0.3/column;
R=0.02;
cen_x=0.1;
cen_y=0.0;

gamma_1 = 1.4;

gamma_2 = 5/3;
rho_b=0.167;
Z_a_out = 1e-3;

Z_a_in  = 1-Z_a_out;

import math
def idx(y,x):
    return math.sqrt((x-cen_x)**2+(y-cen_y)**2) <= R;

rho_a_2=1.29;
u_a_2=0.0;
p_a_2=0.101325;

u_b=u_a_2;
p_b=p_a_2;

M=1.5;
f=1.0/(2.0/(gamma_1+1.0)/M/M+(gamma_1-1.0)/(gamma_1+1.0));
g=2*gamma_1/(gamma_1+1.0)*M*M-(gamma_1-1.0)/(gamma_1+1.0);
rho_a_1=rho_a_2*f;
u_a_1=(1.0-1.0/f)*(u_a_2+math.sqrt(gamma_1*p_a_2/rho_a_2)*M)+u_a_2/f;
u_a_1=-u_a_1;
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
    for i in range(shock,column):
       Z_a[j,i]   = Z_a_out;
       RHO_a[j,i] = rho_a_1;
       U_a[j,i]   = u_a_1;
       V_a[j,i]   = 0.0;
       V_b[j,i]   = 0.0;
       P_a[j,i]   = p_a_1;
       RHO_b[j,i] = rho_b;
       U_b[j,i]   = u_a_1;
       V_b[j,i]   = 0.0;
       P_b[j,i]   = p_b;

for j in range(0,line):
    for i in range(0,shock):
       Z_a[j,i]   = Z_a_out;
       RHO_a[j,i] = rho_a_2;
       U_a[j,i]   = u_a_2;
       V_a[j,i]   = 0.0;
       P_a[j,i]   = p_a_2;
for j in range(0,line):
    for i in range(0,shock):
       RHO_b[j,i] = rho_b;
       U_b[j,i]   = u_b;
       V_b[j,i]   = 0.0;
       P_b[j,i]   = p_b;

for j in range(0,line):
    for i in range(0,shock):
        if idx((j-Y_C+0.5)*d_y, (i-X_C+0.5)*d_x) and idx((j-Y_C+1.5)*d_y, (i-X_C+0.5)*d_x) and idx((j-Y_C+0.5)*d_y, (i-X_C+1.5)*d_x) and idx((j-Y_C+1.5)*d_y, (i-X_C+1.5)*d_x):
           Z_a[j,i] = Z_a_in;
        elif idx((j-Y_C+0.5)*d_y, (i-X_C+0.5)*d_x) or idx((j-Y_C+1.5)*d_y, (i-X_C+0.5)*d_x) or idx((j-Y_C+0.5)*d_y, (i-X_C+1.5)*d_x) or idx((j-Y_C+1.5)*d_y, (i-X_C+1.5)*d_x):
           cell_int,err[j,i] = dblquad(idx,(i-X_C+0.5)*d_x,(i-X_C+1.5)*d_x,lambda y:(j-Y_C+0.5)*d_y,lambda y:(j-Y_C+1.5)*d_y);
           Z_a[j,i] = cell_int/d_x/d_y*Z_a_in+(1.0-cell_int/d_x/d_y)*Z_a_out;


np.savetxt('Z_a.dat',  Z_a,  fmt="%.10g",delimiter="\t");
np.savetxt('RHO.dat',  RHO_a,fmt="%.10g",delimiter="\t");
np.savetxt('U.dat',    U_a,  fmt="%.10g",delimiter="\t");
np.savetxt('V.dat',    V_a,  fmt="%.10g",delimiter="\t");
np.savetxt('P.dat',    P_a,  fmt="%.10g",delimiter="\t");
np.savetxt('RHO_b.dat',RHO_b,fmt="%.10g",delimiter="\t");
np.savetxt('U_b.dat',  U_b,  fmt="%.10g",delimiter="\t");
np.savetxt('V_b.dat',  V_b,  fmt="%.10g",delimiter="\t");
np.savetxt('P_b.dat',  P_b,  fmt="%.10g",delimiter="\t");
