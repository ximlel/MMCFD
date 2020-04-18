#include <stdio.h>
#include <math.h>

#include "../include/var_struc.h"



void cons_qty_copy_cv2ifv(struct i_f_var * ifv, const struct cell_var * cv, const int c)
{
	const int dim = (int)config[0];
	
	ifv->U_rho = cv->U_rho[c];
	ifv->U_e = cv->U_e[c];
	ifv->U_u = cv->U_u[c];
	if (dim > 1)
		ifv->U_v = cv->U_v[c];
	if (dim > 2)
		ifv->U_w = cv->U_w[c];
	if ((int)config[2] == 2)
		{					
			ifv->U_phi = cv->U_phi[c];
			ifv->U_e_a = cv->U_e_a[c];
		}
	ifv->U_gamma = cv->U_gamma[c];
}

void cons_qty_copy_ifv2cv(const struct i_f_var * ifv, struct cell_var * cv, const int c)
{
	const int dim = (int)config[0];
	
	cv->U_rho[c] = ifv->U_rho;
	cv->U_e[c] = ifv->U_e;
	cv->U_u[c] = ifv->U_u;
	if (dim > 1)
		cv->U_v[c] = ifv->U_v;
	if (dim > 2)
		cv->U_w[c] = ifv->U_w;
	if ((int)config[2] == 2)
		{					
			cv->U_phi[c] = ifv->U_phi;
			cv->U_e_a[c] = ifv->U_e_a;
		}
	cv->U_gamma[c] = ifv->U_gamma;
}

void prim_var_copy_ifv2FV(const struct i_f_var * ifv, const struct flu_var * FV, const int c)
{
	const int dim = (int)config[0];

	FV->RHO[c] = ifv->RHO;
	FV->P[c]   = ifv->P;
	FV->U[c]   = ifv->U;
	if (dim > 1)
		FV->V[c] = ifv->V;
	if (dim > 2)
		FV->W[c] = ifv->W;
	if ((int)config[2] == 2)
		{					
			FV->PHI[c] = ifv->PHI;
			FV->Z_a[c] = ifv->Z_a;
		}
	FV->gamma[c] = ifv->gamma;
}

void flux_copy_ifv2cv(const struct i_f_var * ifv, const struct cell_var * cv, const int k, const int j)
{
	const int dim = (int)config[0];
	const int order = (int)config[9];
	
	cv->F_rho[k][j] = ifv->F_rho;
	cv->F_e[k][j]   = ifv->F_e;
	cv->F_u[k][j]   = ifv->F_u;
	if (dim > 1)
		cv->F_v[k][j] = ifv->F_v;
	if (dim > 2)
		cv->F_w[k][j] = ifv->F_w;
	if ((int)config[2] == 2)
		{					
			cv->F_phi[k][j] = ifv->F_phi;
			cv->F_e_a[k][j] = ifv->F_e_a;
		}
	if (!isinf(config[60]))
		cv->F_gamma[k][j] = ifv->F_gamma;

//	cv->RHO_p[k][j] = ifv->RHO;
//	cv->U_p[k][j]   = ifv->U;
//	if (dim > 1)
//		cv->V_p[k][j] = ifv->V;
//	cv->P_p[k][j] = ifv->P;	
	cv->RHO_p[k][j] = ifv->RHO_int;
	cv->U_p[k][j]   = ifv->U_int;
	cv->V_p[k][j]   = ifv->V_int;
	cv->P_p[k][j]   = ifv->P_int;	
	if ((int)config[2] == 2)
		{
			cv->PHI_p[k][j] = ifv->PHI;
			cv->Z_a_p[k][j] = ifv->Z_a;
		}
	cv->gamma_p[k][j] = ifv->gamma;


	cv->RHO_star[k][j]     = ifv->RHO_star;
	cv->P_star[k][j]       = ifv->P_star;
	cv->U_qt_star[k][j]    = ifv->U_qt_star;
	cv->V_qt_star[k][j]    = ifv->V_qt_star;
	cv->gamma_star[k][j]   = ifv->gamma_star;

	cv->RHO_add_c[k][j]    = ifv->RHO_add_c;
	cv->P_add_c[k][j]      = ifv->P_add_c;
	cv->U_qt_add_c[k][j]   = ifv->U_qt_add_c;
	cv->V_qt_add_c[k][j]   = ifv->V_qt_add_c;
	cv->gamma_add_c[k][j]  = ifv->gamma_add_c;

	cv->RHO_minus_c[k][j]  = ifv->RHO_minus_c;
	cv->P_minus_c[k][j]    = ifv->P_minus_c;
	cv->U_qt_minus_c[k][j] = ifv->U_qt_minus_c;
	cv->V_qt_minus_c[k][j] = ifv->V_qt_minus_c;
	cv->gamma_minus_c[k][j] = ifv->gamma_minus_c;

	cv->u_star[k][j]  = ifv->u_star;
	cv->u_minus_c[k][j] = ifv->u_minus_c;
	cv->u_add_c[k][j] = ifv->u_add_c;
}

void flux_add_ifv2cv(const struct i_f_var * ifv, const struct cell_var * cv, const int k, const int j)
{
	const int dim = (int)config[0];
	const int order = (int)config[9];

	cv->F_rho[k][j] = 0.5*(cv->F_rho[k][j] + ifv->F_rho);
	cv->F_e[k][j]   = 0.5*(cv->F_e[k][j]   + ifv->F_e);
	cv->F_u[k][j]   = 0.5*(cv->F_u[k][j]   + ifv->F_u);
	if (dim > 1)
		cv->F_v[k][j] = 0.5*(cv->F_v[k][j] + ifv->F_v);
	if (dim > 2)
		cv->F_w[k][j] = 0.5*(cv->F_w[k][j] + ifv->F_w);
	if ((int)config[2] == 2)
		{					
			cv->F_phi[k][j] = 0.5*(cv->F_phi[k][j] + ifv->F_phi);
			cv->F_e_a[k][j] = 0.5*(cv->F_e_a[k][j] + ifv->F_e_a);
		}
	if (!isinf(config[60]))
		cv->F_gamma[k][j] = 0.5*(cv->F_gamma[k][j] + ifv->F_gamma);

	cv->RHO_p[k][j] = 0.5*(cv->RHO_p[k][j] + ifv->RHO);
	cv->U_p[k][j]   = 0.5*(cv->U_p[k][j] + ifv->U);
	if (dim > 1)
		cv->V_p[k][j] = 0.5*(cv->V_p[k][j] + ifv->V);
	cv->P_p[k][j] = 0.5*(cv->P_p[k][j] + ifv->P);	
	if ((int)config[2] == 2)
		{
			cv->PHI_p[k][j] = 0.5*(cv->PHI_p[k][j] + ifv->PHI);
			cv->Z_a_p[k][j] = 0.5*(cv->Z_a_p[k][j] + ifv->Z_a);
		}
	cv->gamma_p[k][j] = 0.5*(cv->gamma_p[k][j] + ifv->gamma);
}
