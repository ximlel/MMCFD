#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "../include/var_struc.h"
#include "../include/tools.h"


//center_var 2 U
void BN_C2U(struct center_var C, double *U, int i, int j, int U_or_V)
{	
	U[0] = C.ZRHO_gC[i][j];
	U[1] = C.RHO_U_gC[i][j];
	U[2] = C.RHO_V_gC[i][j];
	U[4] = C.ZRHO_sC[i][j];
	U[5] = C.RHO_U_sC[i][j];
	U[6] = C.RHO_V_sC[i][j];
	switch(U_or_V)
		{		   		
		case 0: //x direction :-V		
			U[3] = C.E_gC[i][j]-0.5*pow(C.RHO_V_gC[i][j],2)/C.ZRHO_gC[i][j];
			U[7] = C.E_sC[i][j]-0.5*pow(C.RHO_V_sC[i][j],2)/C.ZRHO_sC[i][j];
			break;
		case 1: //y direction :-U		
			U[3] = C.E_gC[i][j]-0.5*pow(C.RHO_U_gC[i][j],2)/C.ZRHO_gC[i][j];
			U[7] = C.E_sC[i][j]-0.5*pow(C.RHO_U_sC[i][j],2)/C.ZRHO_sC[i][j];
			break;
		}
}

//U_L + U_R 2 primitive_var in center_var
void BN_ULR2prim(struct U_var U_L[][(int)config[13]], struct U_var U_R[][(int)config[13]], struct center_var C, int i, int j)
{	
	C.RHO_gC[i][j] = 0.5*(U_L[i][j].rho_g+U_R[i][j].rho_g);
	C.U_gC[i][j]   = 0.5*(U_L[i][j].u_g +U_R[i][j].u_g);
	C.V_gC[i][j]   = 0.5*(U_L[i][j].v_g +U_R[i][j].v_g);
	C.P_gC[i][j]   = 0.5*(U_L[i][j].p_g +U_R[i][j].p_g);
	C.RHO_sC[i][j] = 0.5*(U_L[i][j].rho_s+U_R[i][j].rho_s);
	C.U_sC[i][j]   = 0.5*(U_L[i][j].u_s +U_R[i][j].u_s);
	C.V_sC[i][j]   = 0.5*(U_L[i][j].v_s +U_R[i][j].v_s);
	C.P_sC[i][j]   = 0.5*(U_L[i][j].p_s +U_R[i][j].p_s);
}

//U_L + U_R 2 conservative_var in center_var
void BN_ULR2cons(struct U_var U_L[][(int)config[13]], struct U_var U_R[][(int)config[13]], struct center_var C, int i, int j)
{	
	C.ZRHO_gC[i][j]  = 0.5*(U_L[i][j].U_rho_g+U_R[i][j].U_rho_g);
	C.RHO_U_gC[i][j] = 0.5*(U_L[i][j].U_v_g +U_R[i][j].U_v_g);
	C.RHO_V_gC[i][j] = 0.5*(U_L[i][j].U_u_g +U_R[i][j].U_u_g);
	C.E_gC[i][j]     = 0.5*(U_L[i][j].U_e_g +U_R[i][j].U_e_g);
	C.ZRHO_sC[i][j]  = 0.5*(U_L[i][j].U_rho_s+U_R[i][j].U_rho_s);
	C.RHO_U_sC[i][j] = 0.5*(U_L[i][j].U_v_s +U_R[i][j].U_v_s);
	C.RHO_V_sC[i][j] = 0.5*(U_L[i][j].U_u_s +U_R[i][j].U_u_s);
	C.E_sC[i][j]     = 0.5*(U_L[i][j].U_e_s +U_R[i][j].U_e_s);
}

void BN_RI2Cx(struct RI_var RI, struct center_var C, int i, int j)
{
	C.Q_xd[i][j]=RI.Q;
 	C.P_xd[i][j]=RI.P;
	C.H_xd[i][j]=RI.H;
	C.eta_g_xd[i][j]=RI.eta_g;
}

void BN_RI2Cy(struct RI_var RI, struct center_var C, int i, int j)
{
	C.Q_yd[i][j]=RI.Q;
 	C.P_yd[i][j]=RI.P;
	C.H_yd[i][j]=RI.H;
	C.eta_g_yd[i][j]=RI.eta_g;
}

void GRP_var_init(struct GRP_LR_var *G, struct slope_var SV, struct U_var U[][(int)config[13]], double d, int i, int j, int p_or_m)
{
	G->rho_gx=SV.RHO_gx[i][j];
	G->p_gx=SV.P_gx[i][j];
	G->u_gx=SV.U_gx[i][j];
	G->v_gx=SV.V_gx[i][j];
	G->rho_gy=SV.RHO_gy[i][j];
	G->p_gy=SV.P_gy[i][j];
	G->u_gy=SV.U_gy[i][j];
	G->v_gy=SV.V_gy[i][j];
	G->rho_sx=SV.RHO_sx[i][j];
	G->p_sx=SV.P_sx[i][j];
	G->u_sx=SV.U_sx[i][j];
	G->v_sx=SV.V_sx[i][j];
	G->rho_sy=SV.RHO_sy[i][j];
	G->p_sy=SV.P_sy[i][j];
	G->u_sy=SV.U_sy[i][j];
	G->v_sy=SV.V_sy[i][j];

	switch(p_or_m)
		{
		case 0:
			G->rho_g =U[i][j].rho_g+d/2*G->rho_gx;
			G->p_g =U[i][j].p_g+d/2*G->p_gx;
			G->u_g =U[i][j].u_g+d/2*G->u_gx;
			G->v_g =U[i][j].v_g+d/2*G->v_gx;
			G->rho_s =U[i][j].rho_s+d/2*G->rho_sx;
			G->p_s =U[i][j].p_s+d/2*G->p_sx;
			G->u_s =U[i][j].u_s+d/2*G->u_sx;
			G->v_s =U[i][j].v_s+d/2*G->v_sx;
		case 1:
			G->rho_g =U[i][j].rho_g-d/2*G->rho_gx;
			G->p_g =U[i][j].p_g-d/2*G->p_gx;
			G->u_g =U[i][j].u_g-d/2*G->u_gx;
			G->v_g =U[i][j].v_g-d/2*G->v_gx;
			G->rho_s =U[i][j].rho_s-d/2*G->rho_sx;
			G->p_s =U[i][j].p_s-d/2*G->p_sx;
			G->u_s =U[i][j].u_s-d/2*G->u_sx;
			G->v_s =U[i][j].v_s-d/2*G->v_sx;
			break;			
		}
}
