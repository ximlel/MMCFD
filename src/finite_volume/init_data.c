#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/var_struc.h"

void init_data_1(double * Z_l, double * U_RHO_g, double * U_U_g, double * U_V_g, double * U_E_g,
                 double * U_RHO_l, double * U_U_l, double * U_V_l, double * U_E_l)
{
	double RHO_g[4], U_g[4], V_g[4], P_g[4], RHO_l[4], U_l[4], V_l[4], P_l[4];
	config[6] = 1.4;
	config[106] = 1.4;
	const double gama_g = config[6], gama_l = config[106];
	Z_l[0] = 0.8;
	Z_l[1] = 0.4;
	Z_l[2] = 0.8;
	Z_l[3] = 0.4;
	RHO_g[0] = 1.5;
	RHO_g[1] = 0.5;
	RHO_g[2] = 1.5;
	RHO_g[3] = 0.5;
	U_g[0] = 0.0;
	U_g[1] = 0.0;
	U_g[2] = 0.0;
	U_g[3] = 0.0;
	V_g[0] = 0.0;
	V_g[1] = 0.0;
	V_g[2] = 0.0;
	V_g[3] = 0.0;
	P_g[0] = 2.0;
	P_g[1] = 1.0;
	P_g[2] = 2.0;
	P_g[3] = 1.0;
	RHO_l[0] = 2.0;
	RHO_l[1] = 1.0;
	RHO_l[2] = 2.0;
	RHO_l[3] = 1.0;
	U_l[0] = 0.0;
	U_l[1] = 0.0;
	U_l[2] = 0.0;
	U_l[3] = 0.0;
	V_l[0] = 0.0;
	V_l[1] = 0.0;
	V_l[2] = 0.0;
	V_l[3] = 0.0;
	P_l[0] = 2.0;
	P_l[1] = 1.0;
	P_l[2] = 2.0;
	P_l[3] = 1.0;
	double Z_g;
	for (int i=0;i<4;i++)
		{
			Z_g = 1.0-Z_l[i];
			U_RHO_g[i] = RHO_g[i]*Z_g;
			U_U_g[i] = U_RHO_g[i]*U_g[i];
			U_V_g[i] = U_RHO_g[i]*V_g[i];
			U_E_g[i] = U_RHO_g[i]*(P_g[i]/RHO_g[i]/(gama_g-1.0)+0.5*U_g[i]*U_g[i]+0.5*V_g[i]*V_g[i]);
			U_RHO_l[i] = RHO_l[i]*Z_l[i];
			U_U_l[i] = U_RHO_l[i]*U_l[i];
			U_V_l[i] = U_RHO_l[i]*V_l[i];
			U_E_l[i] = U_RHO_l[i]*(P_l[i]/RHO_l[i]/(gama_l-1.0)+0.5*U_l[i]*U_l[i]+0.5*V_l[i]*V_l[i]);
		}
}

void init_data_2(double * Z_l, double * U_RHO_g, double * U_U_g, double * U_V_g, double * U_E_g,
                 double * U_RHO_l, double * U_U_l, double * U_V_l, double * U_E_l)
{
	double RHO_g[4], U_g[4], V_g[4], P_g[4], RHO_l[4], U_l[4], V_l[4], P_l[4];
	config[6] = 1.67;
	config[106] = 1.4;
	const double gama_g = config[6], gama_l = config[106];
	Z_l[0] = 0.8;
	Z_l[1] = 0.4;
	Z_l[2] = Z_l[0];
	Z_l[3] = Z_l[1];
	RHO_g[0] = 1.5;
	RHO_g[1] = 0.5;
	U_g[0] = 0.0;
	U_g[1] = 0.0;
	V_g[0] = 0.0;
	V_g[1] = 0.0;
	P_g[0] = 2.0;
	P_g[1] = 1.0;
	RHO_l[0] = 2.0;
	RHO_l[1] = 1.0;
	U_l[0] = 0.0;
	U_l[1] = 0.0;
	V_l[0] = 0.0;
	V_l[1] = 0.0;
	P_l[0] = 2.0;
	P_l[1] = 1.0;
	double Z_g;
	int i;
	for (int j=0;j<4;j++)
		{
			if (j<2)
				i=j;
			else
				i=j-2;	
			Z_g = 1.0-Z_l[i];
			U_RHO_g[j] = RHO_g[i]*Z_g;
			U_U_g[j] = U_RHO_g[j]*U_g[i];
			U_V_g[j] = U_RHO_g[j]*V_g[i];
			U_E_g[j] = U_RHO_g[j]*(P_g[i]/RHO_g[i]/(gama_g-1.0)+0.5*U_g[i]*U_g[i]+0.5*V_g[i]*V_g[i]);
			U_RHO_l[j] = RHO_l[i]*Z_l[i];
			U_U_l[j] = U_RHO_l[j]*U_l[i];
			U_V_l[j] = U_RHO_l[j]*V_l[i];
			U_E_l[j] = U_RHO_l[j]*(P_l[i]/RHO_l[i]/(gama_l-1.0)+0.5*U_l[i]*U_l[i]+0.5*V_l[i]*V_l[i]);
		}
}

void init_data_3(double * Z_l, double * U_RHO_g, double * U_U_g, double * U_V_g, double * U_E_g,
                 double * U_RHO_l, double * U_U_l, double * U_V_l, double * U_E_l)
{
	double RHO_g[4], U_g[4], V_g[4], P_g[4], RHO_l[4], U_l[4], V_l[4], P_l[4];
	config[6] = 1.4;
	config[106] = 1.6;
	const double gama_g = config[6], gama_l = config[106];
	Z_l[0] = 0.6;
	Z_l[1] = 0.4;
	Z_l[2] = Z_l[0];
	Z_l[3] = Z_l[1];
	RHO_g[0] = 3.0;
	RHO_g[1] = 1.2;
	U_g[0] = 0.0;
	U_g[1] = 0.0;
	V_g[0] = 0.0;
	V_g[1] = 0.0;
	P_g[0] = 0.3;
	P_g[1] = 0.12;
	RHO_l[0] = 0.5;
	RHO_l[1] = 1.2;
	U_l[0] = 0.0;
	U_l[1] = 0.0;
	V_l[0] = 0.0;
	V_l[1] = 0.0;
	P_l[0] = 0.6;
	P_l[1] = 0.12;
	double Z_g;
	int i;
	for (int j=0;j<4;j++)
		{
			if (j<2)
				i=j;
			else
				i=j-2;	
			Z_g = 1.0-Z_l[i];
			U_RHO_g[j] = RHO_g[i]*Z_g;
			U_U_g[j] = U_RHO_g[j]*U_g[i];
			U_V_g[j] = U_RHO_g[j]*V_g[i];
			U_E_g[j] = U_RHO_g[j]*(P_g[i]/RHO_g[i]/(gama_g-1.0)+0.5*U_g[i]*U_g[i]+0.5*V_g[i]*V_g[i]);
			U_RHO_l[j] = RHO_l[i]*Z_l[i];
			U_U_l[j] = U_RHO_l[j]*U_l[i];
			U_V_l[j] = U_RHO_l[j]*V_l[i];
			U_E_l[j] = U_RHO_l[j]*(P_l[i]/RHO_l[i]/(gama_l-1.0)+0.5*U_l[i]*U_l[i]+0.5*V_l[i]*V_l[i]);
		}
}


void U_init(int n_x, double ZRHO_gC[][n_x], double RHO_U_gC[][n_x], double RHO_V_gC[][n_x], double E_gC[][n_x],
            double ZRHO_lC[][n_x], double RHO_U_lC[][n_x], double RHO_V_lC[][n_x], double E_lC[][n_x], int i, int j, 
            double * U_RHO_g, double * U_U_g, double * U_V_g, double * U_E_g, 
            double * U_RHO_l, double * U_U_l, double * U_V_l, double * U_E_l, int a, int b, int c, int d)
{
	ZRHO_gC[i][j]  = 0.25*(U_RHO_g[a]+U_RHO_g[b]+U_RHO_g[c]+U_RHO_g[d]);
	RHO_U_gC[i][j] = 0.25*(U_U_g[a]+U_U_g[b]+U_U_g[c]+U_U_g[d]);
	RHO_V_gC[i][j] = 0.25*(U_V_g[a]+U_V_g[b]+U_V_g[c]+U_V_g[d]);
	E_gC[i][j]     = 0.25*(U_E_g[a]+U_E_g[b]+U_E_g[c]+U_E_g[d]);
	ZRHO_lC[i][j]  = 0.25*(U_RHO_l[a]+U_RHO_l[b]+U_RHO_l[c]+U_RHO_l[d]);
	RHO_U_lC[i][j] = 0.25*(U_U_l[a]+U_U_l[b]+U_U_l[c]+U_U_l[d]);
	RHO_V_lC[i][j] = 0.25*(U_V_l[a]+U_V_l[b]+U_V_l[c]+U_V_l[d]);
	E_lC[i][j]     = 0.25*(U_E_l[a]+U_E_l[b]+U_E_l[c]+U_E_l[d]);
}

