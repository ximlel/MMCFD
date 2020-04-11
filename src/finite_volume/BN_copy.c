#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/var_struc.h"
#include "../include/tools.h"
#include "../include/finite_volume.h"

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

void GRP_var_init(struct GRP_LR_var *G, struct slope_var SV, struct U_var U[][(int)config[13]], double d, int i, int j, int pm_xy)
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

	switch(pm_xy)
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
		case 2:
			G->rho_g =U[i][j].rho_g+d/2*G->rho_gy;
			G->p_g =U[i][j].p_g+d/2*G->p_gy;
			G->u_g =U[i][j].u_g+d/2*G->u_gy;
			G->v_g =U[i][j].v_g+d/2*G->v_gy;
			G->rho_s =U[i][j].rho_s+d/2*G->rho_sy;
			G->p_s =U[i][j].p_s+d/2*G->p_sy;
			G->u_s =U[i][j].u_s+d/2*G->u_sy;
			G->v_s =U[i][j].v_s+d/2*G->v_sy;
			break;
		case 3:
			G->rho_g =U[i][j].rho_g-d/2*G->rho_gy;
			G->p_g =U[i][j].p_g-d/2*G->p_gy;
			G->u_g =U[i][j].u_g-d/2*G->u_gy;
			G->v_g =U[i][j].v_g-d/2*G->v_gy;
			G->rho_s =U[i][j].rho_s-d/2*G->rho_sy;
			G->p_s =U[i][j].p_s-d/2*G->p_sy;
			G->u_s =U[i][j].u_s-d/2*G->u_sy;
			G->v_s =U[i][j].v_s-d/2*G->v_sy;
			break;
		}
}
	
void GRP_RI_var_init(struct GRP_RI_LR_var *GRI, struct slope_var SV, struct center_var C, double d, int i, int j, int pm_xy)
{
	GRI->Qx=SV.Q_x[i][j];
	GRI->Px=SV.P_x[i][j];
	GRI->Hx=SV.H_x[i][j];
	GRI->eta_gx=SV.eta_g_x[i][j];
	GRI->Qy=SV.Q_y[i][j];
	GRI->Py=SV.P_y[i][j];
	GRI->Hy=SV.H_y[i][j];
	GRI->eta_gy=SV.eta_g_y[i][j];
	switch(pm_xy)
		{
		case 0:
			GRI->Q     =C.Q_xd[i][j]+d/2*GRI->Qx;
			GRI->P     =C.P_xd[i][j]+d/2*GRI->Px;
			GRI->H     =C.H_xd[i][j]+d/2*GRI->Hx;
			GRI->eta_g =C.eta_g_xd[i][j]+d/2*GRI->eta_gx;
			break;
		case 1:
			GRI->Q     =C.Q_xd[i][j]-d/2*GRI->Qx;
			GRI->P     =C.P_xd[i][j]-d/2*GRI->Px;
			GRI->H     =C.H_xd[i][j]-d/2*GRI->Hx;
			GRI->eta_g =C.eta_g_xd[i][j]-d/2*GRI->eta_gx;			
			break;			
		case 2:
			GRI->Q     =C.Q_yd[i][j]+d/2*GRI->Qy;
			GRI->P     =C.P_yd[i][j]+d/2*GRI->Py;
			GRI->H     =C.H_yd[i][j]+d/2*GRI->Hy;
			GRI->eta_g =C.eta_g_yd[i][j]+d/2*GRI->eta_gy;
			break;
		case 3:
			GRI->Q     =C.Q_yd[i][j]-d/2*GRI->Qy;
			GRI->P     =C.P_yd[i][j]-d/2*GRI->Py;
			GRI->H     =C.H_yd[i][j]-d/2*GRI->Hy;
			GRI->eta_g =C.eta_g_yd[i][j]-d/2*GRI->eta_gy;			
			break;
		}	
}

void RI_LR2G_LR(const struct GRP_RI_LR_var *GRI, struct GRP_LR_var *G, double z_s, int x_or_y)
{	
	struct U_var U;
	struct RI_var RI;
	RI.Q=GRI->Q;
	RI.P=GRI->P;
	RI.H=GRI->H;
	RI.eta_g=GRI->eta_g;
	RI.z_s=z_s;
	RI.rho_s=G->rho_s;
	switch(x_or_y)
		{
		case 0:		
			RI.u_s=G->u_s;
			break;
		case 1:
			RI.u_s=G->v_s;
			break;
		}
	RI2U_cal(&U, &RI, z_s, G->rho_g);
	G->rho_g =U.rho_g;
	G->p_g =U.p_g;
	G->rho_s =U.rho_s;
	G->p_s =U.p_s;
	switch(x_or_y)
		{
		case 0:		
			G->u_g =U.u_g;
			G->v_g =G->v_g;
			G->u_s =U.u_g;
			G->v_s =G->v_g;
			break;
		case 1:
			G->v_g =U.u_g;
			G->u_g =G->u_g;
			G->v_s =U.u_g;
			G->u_s =G->u_g;
			break;			
		}
	/*	
	G->rho_gx =0.0;
	G->p_gx =0.0;
	G->u_gx =0.0;
	G->v_gx =0.0;
	G->rho_sx =0.0;
	G->p_sx =0.0;
	G->u_sx =0.0;
	G->v_sx =0.0;
	G->rho_gy =0.0;
	G->p_gy =0.0;
	G->u_gy =0.0;
	G->v_gy =0.0;
	G->rho_sy =0.0;
	G->p_sy =0.0;
	G->u_sy =0.0;
	G->v_sy =0.0;
	*/
}
