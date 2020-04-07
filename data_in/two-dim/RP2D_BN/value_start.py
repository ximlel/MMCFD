#!/usr/bin/python3

line=200;
column=200;
line_discon=100;
column_discon=100;

#gamma_1 = 1.4
#gamma_2 = 4.4

phi_1=1.0;
rho_1=0.5313;
u_1=0.0;
v_1=0.0;
p_1=0.4;

phi_2=1.0;
rho_2=1.0;
u_2=0.7276;
v_2=0.0;
p_2=1.0;

phi_3=1.0;
rho_3=0.8;
u_3=0.0;
v_3=0.0;
p_3=1.0;

phi_4=1.0;
rho_4=1.0;
u_4=0.0;
v_4=0.7276;
p_4=1.0;

import numpy as np

phi=np.zeros((line,column));
rho=np.zeros((line,column));
u  =np.zeros((line,column));
v  =np.zeros((line,column));
p  =np.zeros((line,column));

for j in range(0,line_discon):
    for i in range(0,column_discon):
        phi[j,i]=phi_3;
        rho[j,i]=rho_3;
        u[j,i]  =u_3;
        v[j,i]  =v_3;
        p[j,i]  =p_3;
    for i in range(column_discon,column):
        phi[j,i]=phi_4;
        rho[j,i]=rho_4;
        u[j,i]  =u_4;
        v[j,i]  =v_4;
        p[j,i]  =p_4;
for j in range(line_discon,line):
    for i in range(0,column_discon):
        phi[j,i]=phi_2;
        rho[j,i]=rho_2;
        u[j,i]  =u_2;
        v[j,i]  =v_2;
        p[j,i]  =p_2;
    for i in range(column_discon,column):
        phi[j,i]=phi_1;
        rho[j,i]=rho_1;
        u[j,i]  =u_1;
        v[j,i]  =v_1;
        p[j,i]  =p_1;
        
np.savetxt('PHI.dat',phi,fmt="%g",delimiter="\t");
np.savetxt('RHO.dat',rho,fmt="%g",delimiter="\t");
np.savetxt('U.dat',  u,  fmt="%g",delimiter="\t");
np.savetxt('V.dat',  v,  fmt="%g",delimiter="\t");
np.savetxt('P.dat',  p,  fmt="%g",delimiter="\t");

import shutil
shutil.copyfile('P.dat', 'Z_a.dat')
