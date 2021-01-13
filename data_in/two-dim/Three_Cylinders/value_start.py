#!/usr/bin/python3

line=100;
column=100;
X_C=line/2;
Y_C=column/2;
R0=1.0;
RC=0.2;
d_x=4.0/line;
d_y=4.0/column;
cen_x1=0.0;
cen_y1=1.0;
import math
cen_x2= math.sqrt(3)/2.0;
cen_y2=-0.5;
cen_x3=-math.sqrt(3)/2.0;
cen_y3=-0.5;

eps = 1e-4;
Z_a_in  = 1.0-eps;
Z_a_out = eps;
gamma_1 = 1.4;
gamma_2 = 1.67;

def idx(y,x):
    return math.sqrt((x-cen_x1)**2+(y-cen_y1)**2) <= RC or math.sqrt((x-cen_x2)**2+(y-cen_y2)**2) <= RC or math.sqrt((x-cen_x3)**2+(y-cen_y3)**2) <= RC
    # if math.sqrt((x-cen_x1)**2+(y-cen_y1)**2) > RC and math.sqrt((x-cen_x2)**2+(y-cen_y2)**2) > RC and math.sqrt((x-cen_x3)**2+(y-cen_y3)**2) > RC:
    #     return 1;
    # else:
    #     return 0;

rho_a = 1.0;
p_a = 1.0;
rho_b = 1.4;
p_b = 1.0;

import numpy as np
from scipy.integrate import dblquad

Z_a  = np.zeros((line,column));
RHO_a = np.zeros((line,column));
U_a  = np.zeros((line,column));
V_a  = np.zeros((line,column));
P_a  = np.zeros((line,column));
RHO_b = np.zeros((line,column));
U_b  = np.zeros((line,column));
V_b  = np.zeros((line,column));
P_b  = np.zeros((line,column));
err  = np.zeros((line,column));

for j in range(0,line):
    for i in range(0,column):
       Z_a[j,i]   = Z_a_out;
       RHO_a[j,i] = rho_a;
       P_a[j,i]   = p_a;
       RHO_b[j,i] = rho_b;
       U_b[j,i]   = 0.0;
       V_b[j,i]   = 0.0;
       P_b[j,i]   = p_b;

for j in range(0,line):
    for i in range(0,column):
        U_a[j,i]   =  3*(j-Y_C+0.5)*d_y;
        V_a[j,i]   = -3*(i-X_C+0.5)*d_x;
        if idx((j-Y_C+0.5)*d_y, (i-X_C+0.5)*d_x) or idx((j-Y_C+1.5)*d_y, (i-X_C+0.5)*d_x) or idx((j-Y_C+0.5)*d_y, (i-X_C+1.5)*d_x) or idx((j-Y_C+1.5)*d_y, (i-X_C+1.5)*d_x):
           cell_int,err[j,i] = dblquad(idx,(i-X_C+0.5)*d_x,(i-X_C+1.5)*d_x,lambda y:(j-Y_C+0.5)*d_y,lambda y:(j-Y_C+1.5)*d_y);
           Z_a[j,i] = cell_int/d_x/d_y*Z_a_in+(1.0-cell_int/d_x/d_y)*Z_a_out;
           

np.savetxt('Z_a.dat',  Z_a,  fmt="%.10g",delimiter="\t");
np.savetxt('RHO.dat',  RHO_b,fmt="%.10g",delimiter="\t");
np.savetxt('U.dat',    U_b,  fmt="%.10g",delimiter="\t");
np.savetxt('V.dat',    V_b,  fmt="%.10g",delimiter="\t");
np.savetxt('P.dat',    P_b,  fmt="%.10g",delimiter="\t");
np.savetxt('RHO_b.dat',RHO_a,fmt="%.10g",delimiter="\t");
np.savetxt('U_b.dat',  U_a,  fmt="%.10g",delimiter="\t");
np.savetxt('V_b.dat',  V_a,  fmt="%.10g",delimiter="\t");
np.savetxt('P_b.dat',  P_a,  fmt="%.10g",delimiter="\t");
