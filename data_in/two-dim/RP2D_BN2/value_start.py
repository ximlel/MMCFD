#!/usr/bin/python3

line=200;
column=200;
line_discon=100;
column_discon=100;

gamma_1 = 1.6;
gamma_2 = 1.4;

z_a_1  =0.6;
rho_a_1=0.5;
u_a_1  =0.0;
v_a_1  =0.0;
p_a_1  =0.6;
rho_b_1=3.0;
u_b_1  =0.0;
v_b_1  =0.0;
p_b_1  =0.3;

z_a_2  =0.4;
rho_a_2=1.2;
u_a_2  =0.0;
v_a_2  =0.0;
p_a_2  =0.12;
rho_b_2=1.2;
u_b_2  =0.0;
v_b_2  =0.0;
p_b_2  =0.12;

import numpy as np

z_a  =np.zeros((line,column));
rho_a=np.zeros((line,column));
u_a  =np.zeros((line,column));
v_a  =np.zeros((line,column));
p_a  =np.zeros((line,column));
rho_b=np.zeros((line,column));
u_b  =np.zeros((line,column));
v_b  =np.zeros((line,column));
p_b  =np.zeros((line,column));

for j in range(0,line_discon):
    for i in range(0,column_discon):
        z_a[j,i]  =z_a_1;
        rho_a[j,i]=rho_a_1;
        u_a[j,i]  =u_a_1;
        v_a[j,i]  =v_a_1;
        p_a[j,i]  =p_a_1;
        rho_b[j,i]=rho_b_1;
        u_b[j,i]  =u_b_1;
        v_b[j,i]  =v_b_1;
        p_b[j,i]  =p_b_1;
    for i in range(column_discon,column):
        z_a[j,i]  =z_a_2;
        rho_a[j,i]=rho_a_2;
        u_a[j,i]  =u_a_2;
        v_a[j,i]  =v_a_2;
        p_a[j,i]  =p_a_2;
        rho_b[j,i]=rho_b_2;
        u_b[j,i]  =u_b_2;
        v_b[j,i]  =v_b_2;
        p_b[j,i]  =p_b_2;
for j in range(line_discon,line):
    for i in range(0,column_discon):
        z_a[j,i]  =z_a_2;
        rho_a[j,i]=rho_a_2;
        u_a[j,i]  =u_a_2;
        v_a[j,i]  =v_a_2;
        p_a[j,i]  =p_a_2;
        rho_b[j,i]=rho_b_2;
        u_b[j,i]  =u_b_2;
        v_b[j,i]  =v_b_2;
        p_b[j,i]  =p_b_2;
    for i in range(column_discon,column):
        z_a[j,i]  =z_a_1;
        rho_a[j,i]=rho_a_1;
        u_a[j,i]  =u_a_1;
        v_a[j,i]  =v_a_1;
        p_a[j,i]  =p_a_1;
        rho_b[j,i]=rho_b_1;
        u_b[j,i]  =u_b_1;
        v_b[j,i]  =v_b_1;
        p_b[j,i]  =p_b_1;
        
np.savetxt('Z_a.dat',  z_a,  fmt="%g",delimiter="\t");
np.savetxt('RHO.dat',  rho_a,fmt="%g",delimiter="\t");
np.savetxt('U.dat',    u_a,  fmt="%g",delimiter="\t");
np.savetxt('V.dat',    v_a,  fmt="%g",delimiter="\t");
np.savetxt('P.dat',    p_a,  fmt="%g",delimiter="\t");
np.savetxt('RHO_b.dat',rho_b,fmt="%g",delimiter="\t");
np.savetxt('U_b.dat',  u_b,  fmt="%g",delimiter="\t");
np.savetxt('V_b.dat',  v_b,  fmt="%g",delimiter="\t");
np.savetxt('P_b.dat',  p_b,  fmt="%g",delimiter="\t");
