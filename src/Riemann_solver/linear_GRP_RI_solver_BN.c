#include <math.h>
#include <stdio.h>

#include "../include/var_struc.h"
#include "../include/Riemann_solver.h"

void linear_GRP_RI_solver_BN
(double *D, const double D_z_s, const double z_s,
 const double *mid_g, const double *mid_s, 
 const struct GRP_RI_LR_var RI_L, const struct GRP_RI_LR_var RI_R,
 const struct GRP_LR_var GL, const struct GRP_LR_var GR,
 const double gamma_s, const double gamma_g,
 const double eps, int x_or_y)
{
	double rho_g = mid_g[0], u_g, p_g = mid_g[3];
	double rho_s = mid_s[0], u_s, p_s = mid_s[3];
	switch(x_or_y)
		{
		case 0:
			u_g = mid_g[1];
			u_s = mid_s[1];
		case 1:
			u_g = mid_g[2];
			u_s = mid_s[2];			
		}
	double c_s, c_g;
	c_s = sqrt(gamma_s * p_s / rho_s);
	c_g = sqrt(gamma_g * p_g / rho_g);
	double Lambda_v_p[7][7]={0.0}, Lambda_v_m[7][7]={0.0};
	Lambda_v_p[0][0]=fmax(u_s,0.0);
	Lambda_v_m[0][0]=fmin(u_s,0.0);
	Lambda_v_p[1][1]=fmax(u_s-c_s,0.0);
	Lambda_v_m[1][1]=fmin(u_s-c_s,0.0);
	Lambda_v_p[2][2]=fmax(u_s,0.0);
	Lambda_v_m[2][2]=fmin(u_s,0.0);
	Lambda_v_p[3][3]=fmax(u_s+c_s,0.0);
	Lambda_v_m[3][3]=fmin(u_s+c_s,0.0);
	Lambda_v_p[4][4]=fmax(u_g-c_g,0.0);
	Lambda_v_m[4][4]=fmin(u_g-c_g,0.0);
	Lambda_v_p[5][5]=fmax(u_g,0.0);
	Lambda_v_m[5][5]=fmin(u_g,0.0);
	Lambda_v_p[6][6]=fmax(u_g+c_g,0.0);
	Lambda_v_m[6][6]=fmin(u_g+c_g,0.0);
	double D_L[7], D_R[7], R[7][7], R_inv[7][7];
	D_L[0] = D_z_s;
	D_R[0] = D_z_s;
	switch(x_or_y)
		{
		case 0:
			D_L[1] = GL.rho_sx;
			D_L[2] = GL.u_sx;
			D_L[3] = RI_L.Px;
			D_L[4] = RI_L.Qx;
			D_L[5] = RI_L.Hx;
			D_L[6] = RI_L.eta_gx;
			D_R[1] = GR.rho_sx;
			D_R[2] = GR.u_sx;
			D_R[3] = RI_R.Px;
			D_R[4] = RI_R.Qx;
			D_R[5] = RI_R.Hx;
			D_R[6] = RI_R.eta_gx;			
		case 1:
			D_L[1] = GL.rho_sy;
			D_L[2] = GL.u_sy;
			D_L[3] = RI_L.Py;
			D_L[4] = RI_L.Qy;
			D_L[5] = RI_L.Hy;
			D_L[6] = RI_L.eta_gy;
			D_R[1] = GR.rho_sy;
			D_R[2] = GR.u_sy;
			D_R[3] = RI_R.Py;
			D_R[4] = RI_R.Qy;
			D_R[5] = RI_R.Hy;
			D_R[6] = RI_R.eta_gy;					
		}
}
