#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/var_struc.h"

void FV_2_C_init(struct center_var C, struct flu_var FV)
{
	const int n_y = (int)config[14], n_x = (int)config[13];
	const double gamma_g=config[6], gamma_s=config[106];
	int i,j,i_1,j_1;
	double Z_g;
    double U_RHO_g[n_y][n_x], U_U_g[n_y][n_x], U_V_g[n_y][n_x], U_E_g[n_y][n_x];
    double U_RHO_s[n_y][n_x], U_U_s[n_y][n_x], U_V_s[n_y][n_x], U_E_s[n_y][n_x];
	for(i = 0; i < n_y; ++i)
		for(j = 0; j < n_x; ++j)
			{
				C.Z_sC[i][j]  =FV.Z_a[i*n_x+j];
				C.RHO_gC[i][j]=FV.RHO[i*n_x+j];
				C.U_gC[i][j]  =FV.U[i*n_x+j];
				C.V_gC[i][j]  =FV.V[i*n_x+j];
				C.P_gC[i][j]  =FV.P[i*n_x+j];
				C.RHO_sC[i][j]=FV.RHO_b[i*n_x+j];
				C.U_sC[i][j]  =FV.U_b[i*n_x+j];
				C.V_sC[i][j]  =FV.V_b[i*n_x+j];
				C.P_sC[i][j]  =FV.P_b[i*n_x+j];
                Z_g = 1.0-C.Z_sC[i][j];
                U_RHO_g[i][j] = C.RHO_gC[i][j]*Z_g;
                U_U_g[i][j] = U_RHO_g[i][j]*C.U_gC[i][j];
                U_V_g[i][j] = U_RHO_g[i][j]*C.V_gC[i][j];
                U_E_g[i][j] = U_RHO_g[i][j]*(C.P_gC[i][j]/C.RHO_gC[i][j]/(gamma_g-1.0)+0.5*pow(C.U_gC[i][j],2)+0.5*pow(C.V_gC[i][j],2));
                U_RHO_s[i][j] = C.RHO_sC[i][j]*C.Z_sC[i][j];
                U_U_s[i][j] = U_RHO_s[i][j]*C.U_sC[i][j];
                U_V_s[i][j] = U_RHO_s[i][j]*C.V_sC[i][j];
                U_E_s[i][j] = U_RHO_s[i][j]*(C.P_sC[i][j]/C.RHO_sC[i][j]/(gamma_s-1.0)+0.5*pow(C.U_sC[i][j],2)+0.5*pow(C.V_sC[i][j],2));
			}
	for(i = 0; i < n_y; ++i)
		for(j = 0; j < n_x; ++j) //ignore "n_x-1, n_y-1" also is OK.
			{
			    i_1=fmax(i-1,0);
			    j_1=fmax(j-1,0);
                C.ZRHO_gC[i][j]  = 0.25*(U_RHO_g[i_1][j_1]+U_RHO_g[i_1][j]+U_RHO_g[i][j_1]+U_RHO_g[i][j]);
                C.RHO_U_gC[i][j] = 0.25*(U_U_g[i_1][j_1]  +U_U_g[i_1][j]  +U_U_g[i][j_1]  +U_U_g[i][j]);
                C.RHO_V_gC[i][j] = 0.25*(U_V_g[i_1][j_1]  +U_V_g[i_1][j]  +U_V_g[i][j_1]  +U_V_g[i][j]);
                C.E_gC[i][j]     = 0.25*(U_E_g[i_1][j_1]  +U_E_g[i_1][j]  +U_E_g[i][j_1]  +U_E_g[i][j]);
                C.ZRHO_sC[i][j]  = 0.25*(U_RHO_s[i_1][j_1]+U_RHO_s[i_1][j]+U_RHO_s[i][j_1]+U_RHO_s[i][j]);
                C.RHO_U_sC[i][j] = 0.25*(U_U_s[i_1][j_1]  +U_U_s[i_1][j]  +U_U_s[i][j_1]  +U_U_s[i][j]);
                C.RHO_V_sC[i][j] = 0.25*(U_V_s[i_1][j_1]  +U_V_s[i_1][j]  +U_V_s[i][j_1]  +U_V_s[i][j]);
                C.E_sC[i][j]     = 0.25*(U_E_s[i_1][j_1]  +U_E_s[i_1][j]  +U_E_s[i][j_1]  +U_E_s[i][j]);
			}
}
