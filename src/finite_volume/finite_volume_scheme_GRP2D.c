#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


#include "../include/var_struc.h"
#include "../include/tools.h"
#include "../include/finite_volume.h"
#include "../include/Riemann_solver.h"
#include "../include/file_io.h"

#define pi (4.*atan(1.0))

static double minmod3(const double a, const double b, const double c)
{
	const double Alpha=1.5;
	if (a>0.0&&b>0.0&&c>0.0)
		return fmin(fmin(Alpha*a,Alpha*c),b);
	if (a<0.0&&b<0.0&&c<0.0)
		return fmax(fmax(Alpha*a,Alpha*c),b);
	else
		return 0.0;
}

static double minmod2(const double a, const double b)
{
	if (a>0.0&&b>0.0)
		return fmin(a,b);
	if (a<0.0&&b<0.0)
		return fmax(a,b);
	else
		return 0.0;
}

#define cent_lim(fvar)													\
	do {																\
		SV.fvar##x[i][j] = (CV.fvar##C[i][jR]-CV.fvar##C[i][jL])/config[10]/2; \
		SV.fvar##y[i][j] = (CV.fvar##C[iR][j]-CV.fvar##C[iL][j])/config[11]/2; \
		SV.fvar##x[i][j] = minmod3((CV.fvar##C[i][j]-CV.fvar##C[i][jL])/config[10],SV.fvar##x[i][j],(CV.fvar##C[i][jR]-CV.fvar##C[i][j])/config[10]); \
		SV.fvar##y[i][j] = minmod3((CV.fvar##C[i][j]-CV.fvar##C[iL][j])/config[11],SV.fvar##y[i][j],(CV.fvar##C[iR][j]-CV.fvar##C[i][j])/config[11]); \
	} while(0)

#define side_lim(fvar)													\
	do {																\
		SV.fvar##x[i][j] = minmod2((CV.fvar##C[i][j]-CV.fvar##C[i][jL])/config[10],(CV.fvar##C[i][jR]-CV.fvar##C[i][j])/config[10]); \
		SV.fvar##y[i][j] = minmod2((CV.fvar##C[i][j]-CV.fvar##C[iL][j])/config[11],(CV.fvar##C[iR][j]-CV.fvar##C[i][j])/config[11]); \
	} while(0)

#define side_lim_RI(fvar)													\
	do {																\
		SV.fvar##x[i][j] = minmod2((CV.fvar##xd[i][j]-CV.fvar##xd[i][jL])/config[10],(CV.fvar##xd[i][jR]-CV.fvar##xd[i][j])/config[10]); \
		SV.fvar##y[i][j] = minmod2((CV.fvar##yd[i][j]-CV.fvar##yd[iL][j])/config[11],(CV.fvar##yd[iR][j]-CV.fvar##yd[i][j])/config[11]); \
	} while(0)

static void slope_simiter3_GRP(int n_x, struct face_var FV, struct center_var CV, struct slope_var SV)
{
	int HN = 2;
	const int order = (int)config[9];
	const double eps = config[4];
	const int n_y = (int)config[14];
	int i,iL,iR, j,jL,jR, k;
	double **p;
	for(i = 0; i < n_y; ++i)
		for(j = 0; j < n_x; ++j)
			{
				if (i==0 || i==1)
					{
						iL=1; iR=1;
					}
				else if (i==n_y-1 || i==n_y-2)
					{
						iL=n_y-2; iR=n_y-2;
					}
				else
					{
						iL=i-1; iR=i+1;
					}
				if (j==0 || j==1)
					{
						jL=1; jR=1;
					}
				else if (j==n_x-1 || j==n_x-2)
					{
						jL=n_x-2; jR=n_x-2;
					}
				else
					{
						jL=j-1; jR=j+1;
					}
				/*
				  cent_lim(Z_s);
				  cent_lim(RHO_g);cent_lim(U_g);cent_lim(V_g);cent_lim(P_g);
				  cent_lim(RHO_s);cent_lim(U_s);cent_lim(V_s);cent_lim(P_s);
				*/
				side_lim(Z_s);
				side_lim(RHO_g);side_lim(U_g);side_lim(V_g);side_lim(P_g);
				side_lim(RHO_s);side_lim(U_s);side_lim(V_s);side_lim(P_s);
				side_lim_RI(Q_);side_lim_RI(P_);side_lim_RI(H_);side_lim_RI(eta_g_);
				side_lim_RI(Z_sS_);
				
				if (order == 1)
					for(k=0, p=SV.Z_sx; k<sizeof(struct slope_var)/sizeof(double **); k++)
						(p++)[i][j] = 0.0;
			}
	for(i = 1; i < n_y-1; ++i)
		for(j = 1; j < n_x-1; ++j)
			{
				iL=i-1; iR=i+1;
				jL=j-1;	jR=j+1;
				if (fabs(CV.Z_sC[iR][j]-CV.Z_sC[i][j])>eps || fabs(CV.Z_sC[i][j]-CV.Z_sC[iL][j])>eps ||
					fabs(CV.Z_sC[i][jR]-CV.Z_sC[i][j])>eps || fabs(CV.Z_sC[i][j]-CV.Z_sC[i][jL])>eps)
					for (k = -HN; k < HN; k++)
						{
							if (i+k >= 0 && i+k < n_y)
								for(k=0, p=SV.Z_sx; k<sizeof(struct slope_var)/sizeof(double **); k++)
									(p++)[i+k][j] = 0.0;
							if (j+k >= 0 && j+k < n_x)
								for(k=0, p=SV.Z_sx; k<sizeof(struct slope_var)/sizeof(double **); k++)
									(p++)[i][j+k] = 0.0;
						}
			}
}

void finite_volume_scheme_GRP2D(struct flu_var * FV, const struct mesh_var * mv, const char * scheme, const char * problem) //R-K
{
	clock_t start_clock;
	double cpu_time = 0.0;

	const int dim = (int)config[0];
	const int order = (int)config[9];
	const double eps = config[4], eps_big = 1e-10;
	const double gamma_g=config[6], gamma_s=config[106];
	const int n_y = (int)config[14], n_x = (int)config[13];
	const double dx= config[10], dy= config[11];
	struct cell_var cv = cell_mem_init(mv, FV);

	vol_comp(&cv, mv);
	cell_rel(&cv, mv);
	printf("Grid has been constructed.\n");

	double tau; // the length of the time step
	double t_all = 0.0;
	const double delta_plot_t = 0.01;
	double plot_t = 0.0;
	int i, j, k, l, stop_step = 0, stop_t = 0;

	double **p;
	struct center_var C;
	struct slope_var SV;
	struct face_var F;
	struct flux_var FX;
	for(k=0, p=SV.Z_sx; k<sizeof(struct slope_var)/sizeof(double **); k++,p++)
		{
			p = (double**)calloc(n_y,sizeof(double*));
			for(i=0; i<n_y; i++)
				p[i]=(double*)calloc(n_x,sizeof(double));
		}
	for(k=0, p=C.Z_sC; k<sizeof(struct center_var)/sizeof(double **); k++,p++)
		{
			p = (double**)calloc(n_y,sizeof(double*));
			for(i=0; i<n_y; i++)
				p[i]=(double*)calloc(n_x,sizeof(double));
		}
	for(k=0, p=F.Z_I_sxL; k<sizeof(struct face_var)/sizeof(double **); k++,p++)
		{
			p = (double**)calloc(n_y+1,sizeof(double*));
			for(i=0; i<n_y+1; i++)
				p[i]=(double*)calloc(n_x+1,sizeof(double));
		}
	for(k=0, p=FX.RHO_F_gx; k<sizeof(struct flux_var)/sizeof(double **); k++,p++)
		{
			p = (double**)calloc(n_y+1,sizeof(double*));
			for(i=0; i<n_y+1; i++)
				p[i]=(double*)calloc(n_x+1,sizeof(double));
		}
	FV_2_C_init(C, *FV);
	
	int iL,iR,jL,jR, i_1, j_1, ip1, jp1;
	double z_gmid, rho_gmid, p_gmid, u_gmid, v_gmid;
	double z_smid, rho_smid, p_smid, u_smid, v_smid;
	double z_sx_mid, z_sy_mid;
	double a_g, a_s, S_max, S_max_g, S_max_s;
	double wave_speed[2], dire[6], mid_s[6], mid_g[6], star[6];
	struct GRP_LR_var GL, GR;
	struct U_var U_L[n_y][n_x], U_R[n_y][n_x], U_tmp;
	struct RI_var RI;
	struct GRP_RI_LR_var RI_L, RI_R;
	double U[8];
	double S, S_tmp, area_L, area_R, RHO_s_cell, ZRHO_s_cell;
	double z_sL, z_sR, z_sxL, z_sxR, z_syL, z_syR;
	double stag_RHO_F_sx[n_y+1][n_x+1], stag_ZRHO_F_sx[n_y+1][n_x+1];
	double stag_RHO_F_sy[n_y+1][n_x+1], stag_ZRHO_F_sy[n_y+1][n_x+1];

	for(l = 0; l < (int)config[5] && stop_step != 1; ++l)
		{
			start_clock = clock();

			if (t_all >= plot_t)
				{
					for(i = 0; i < n_y; ++i)
						for(j = 0; j < n_x; ++j)
							{
								FV->RHO[i*n_x+j]=C.RHO_gC[i][j];
								FV->U[i*n_x+j]  =C.U_gC[i][j];
								FV->V[i*n_x+j]  =C.V_gC[i][j];
								FV->P[i*n_x+j]  =C.P_gC[i][j];
								FV->Z_a[i*n_x+j]=C.Z_sC[i][j];
							}
					file_write_TEC(*FV, *mv, problem, plot_t, dim);
					plot_t += delta_plot_t;
				}

			if (stop_step == 0)
				{
					for(i = 0; i < n_y; ++i)
						for(j = 0; j < n_x; ++j)
							{
								i_1 = i-1>0?i-1:0;
								j_1 = j-1>0?j-1:0;
								C.Z_sS_xd[i][j] = 0.5*(C.Z_sC[i_1][j]  +C.Z_sC[i][j]);
								BN_C2U(C,U,i,j,0);
								primitive_comp(U, &U_L[i][j], &U_R[i][j], C.Z_sS_xd[i][j_1], C.Z_sS_xd[i][j], C.Z_sS_xd[i][j_1], C.Z_sS_xd[i][j], 0.5, 0.5);
								BN_ULR2prim(U_L,U_R,C,i,j);
								U2RI_cal(&U_L[i][j], &RI);
								BN_RI2Cx(RI,C,i,j);
							}
				}
			else if (stop_step == 2)
				{
					for(j = 0; j < n_x; ++j)
						for(i = 0; i < n_y; ++i)
							{
								i_1 = i-1>0?i-1:0;
								j_1 = j-1>0?j-1:0;
								C.Z_sS_yd[i][j] = 0.5*(C.Z_sC[i][j_1]+C.Z_sC[i][j]);
								BN_C2U(C,U,i,j,1);
								primitive_comp(U, &U_L[i][j], &U_R[i][j], C.Z_sS_yd[i_1][j], C.Z_sS_yd[i][j], C.Z_sS_yd[i_1][j], C.Z_sS_yd[i][j], 0.5, 0.5);
								BN_ULR2prim(U_L,U_R,C,i,j);
								U2RI_cal(&U_L[i][j], &RI);
								BN_RI2Cy(RI,C,i,j);
							}
				}
			slope_simiter3_GRP(n_x, F, C, SV);

			if (stop_step == 0)
				{
					tau = 1e15;
					for(i = 0; i < n_y; ++i)
						for(j = 0; j < n_x; ++j)
							{
								a_g     = sqrt(gamma_g*C.P_gC[i][j]/C.RHO_gC[i][j]);
								a_s     = sqrt(gamma_s*C.P_sC[i][j]/C.RHO_sC[i][j]);
								S_max_g = fmax(fmax(fabs(C.U_gC[i][j]-a_g),fabs(C.U_gC[i][j]+a_g)),fmax(fabs(C.V_gC[i][j]-a_g),fabs(C.V_gC[i][j]+a_g)));
								S_max_s = fmax(fmax(fabs(C.U_sC[i][j]-a_s),fabs(C.U_sC[i][j]+a_s)),fmax(fabs(C.V_sC[i][j]-a_s),fabs(C.V_sC[i][j]+a_s)));
								S_max = fmax(S_max_g,S_max_s);
								tau   = fmin(tau,dx*config[7]/S_max);
							}
				}

			if((t_all + tau + eps) > config[1])
				{
					tau = config[1] - t_all;
					if (stop_step == 2)
						{
							printf("\nThe time is enough at step %d.\n",l);
							stop_t = 1;
						}
				} // Time
			if (stop_step == 2)
				t_all += tau;

			if (stop_step == 0)
				{
					for(i = 0; i <  n_y; ++i)
						for(j = 0; j <= n_x; ++j)
							{
								jL = j-1>0?j-1:0;
								jR = j<n_x?j:n_x-1;
								z_smid = U_R[i][jL].z_s;
								z_gmid = 1.0-z_smid;
								z_sx_mid = SV.Z_sS_x[i][j];
								GRP_var_init(&GL, SV, U_R, dx, i, jL, 0);
								GRP_var_init(&GR, SV, U_L, dx, i, jR, 1);
								GRP_RI_var_init(&RI_L, SV, C, dx, i, jL, 0);
								GRP_RI_var_init(&RI_R, SV, C, dx, i, jR, 1);
								RI_LR2G_LR(&RI_L,&GL,z_smid,0);
								RI_LR2G_LR(&RI_R,&GR,z_smid,0);
								linear_GRP_solver_Edir_Q1D(wave_speed, dire, mid_g, star, 0.0, 0.0, GL.rho_g, GR.rho_g, GL.rho_gx, GR.rho_gx, GL.rho_gy, GR.rho_gy, GL.u_g, GR.u_g, GL.u_gx, GR.u_gx, GL.u_gy, GR.u_gy, GL.v_g, GR.v_g, GL.v_gx, GR.v_gx, GL.v_gy, GR.v_gy, GL.p_g, GR.p_g, GL.p_gx, GR.p_gx, GL.p_gy, GR.p_gy, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, gamma_g, gamma_g, eps, eps);
								rho_gmid = mid_g[0] + 0.5*tau*dire[0];
								u_gmid   = mid_g[1] + 0.5*tau*dire[1];
								v_gmid   = mid_g[2] + 0.5*tau*dire[2];
								p_gmid   = mid_g[3] + 0.5*tau*dire[3];
								linear_GRP_solver_Edir_Q1D(wave_speed, dire, mid_s, star, 0.0, 0.0, GL.rho_s, GR.rho_s, GL.rho_sx, GR.rho_sx, GL.rho_sy, GR.rho_sy, GL.u_s, GR.u_s, GL.u_sx, GR.u_sx, GL.u_sy, GR.u_sy, GL.v_s, GR.v_s, GL.v_sx, GR.v_sx, GL.v_sy, GR.v_sy, GL.p_s, GR.p_s, GL.p_sx, GR.p_sx, GL.p_sy, GR.p_sy, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, gamma_s, gamma_s, eps, eps);
								rho_smid = mid_s[0] + 0.5*tau*dire[0];
								u_smid   = mid_s[1] + 0.5*tau*dire[1];
								v_smid   = mid_s[2] + 0.5*tau*dire[2];
								p_smid   = mid_s[3] + 0.5*tau*dire[3];
								linear_GRP_RI_solver_BN(&RI, z_sx_mid, z_smid, mid_g, mid_s, RI_L, RI_R, GL, GR, gamma_s, gamma_g, eps, tau, 0);
								RI2U_cal(&U_tmp, &RI, RI.z_s, rho_gmid);
								rho_gmid = U_tmp.rho_g;
								u_gmid   = U_tmp.u_g;
								v_gmid   = mid_g[2];
								p_gmid   = U_tmp.p_g;
								rho_smid = U_tmp.rho_s;
								u_smid   = U_tmp.u_s;
								v_smid   = mid_s[2];
								p_smid   = U_tmp.p_s;

								FX.RHO_F_gx[i][j] = z_gmid*rho_gmid*u_gmid;
								FX.U_F_gx[i][j]   = FX.RHO_F_gx[i][j]*u_gmid + z_gmid*p_gmid;
								FX.V_F_gx[i][j]   = FX.RHO_F_gx[i][j]*v_gmid;
								FX.E_F_gx[i][j]   = gamma_g/(gamma_g-1.0)*p_gmid/rho_gmid + 0.5*(u_gmid*u_gmid + v_gmid*v_gmid);
								FX.E_F_gx[i][j]   = FX.RHO_F_gx[i][j]*FX.E_F_gx[i][j];
								FX.RHO_F_sx[i][j] = z_smid*rho_smid*u_smid;
								FX.U_F_sx[i][j]   = FX.RHO_F_sx[i][j]*u_smid + z_smid*p_smid;
								FX.V_F_sx[i][j]   = FX.RHO_F_sx[i][j]*v_smid;
								FX.E_F_sx[i][j]   = gamma_s/(gamma_s-1.0)*p_smid/rho_smid + 0.5*(u_smid*u_smid + v_smid*v_smid);
								FX.E_F_sx[i][j]   = FX.RHO_F_sx[i][j]*FX.E_F_sx[i][j];
								FX.P_s_MIDx[i][j] = p_smid;
								FX.Z_s_MIDx[i][j] = z_smid;
							}
					for(i = 0; i < n_y; ++i)
						for(j = 0; j <=n_x; ++j)
							{
								jL = j-1>0?j-1:0;
								jR = j<n_x?j:n_x-1;
								z_sxL = SV.Z_sx[i][jL];
								z_sxR = SV.Z_sx[i][jR];
								z_syL = SV.Z_sy[i][jL];
								z_syR = SV.Z_sy[i][jR];
								z_sL = C.Z_sC[i][jL]+dx/2*z_sxL;
								z_sR = C.Z_sC[i][jR]-dx/2*z_sxR;
								ip1 = i<n_y-1?i+1:n_y-1;
								stag_RHO_F_sx[i][j] = 0.5*(U_R[i][j].rho_s*U_R[i][j].u_s+U_R[ip1][j].rho_s*U_R[ip1][j].u_s);
								if ((U_R[i][j].u_s + U_R[ip1][j].u_s) > 0.0)
									stag_ZRHO_F_sx[i][j] = stag_RHO_F_sx[i][j]*(z_sL - 0.25*tau*((U_R[i][j].u_s + U_R[ip1][j].u_s)*z_sxL+(U_R[i][j].v_s + U_R[ip1][j].v_s)*z_syL));
								else
									stag_ZRHO_F_sx[i][j] = stag_RHO_F_sx[i][j]*(z_sR - 0.25*tau*((U_R[i][j].u_s + U_R[ip1][j].u_s)*z_sxR+(U_R[i][j].v_s + U_R[ip1][j].v_s)*z_syR));
							}
					for(i = 0; i < n_y; ++i)
						for(j = 0; j < n_x; ++j)
							{
								ip1 = i<n_y-1?i+1:n_y-1;
								jp1 = j<n_x-1?j+1:n_x-1;
								RHO_s_cell  = 0.25*(C.RHO_sC[i][j]+C.RHO_sC[ip1][j]+C.RHO_sC[i][jp1]+C.RHO_sC[ip1][jp1]);
								//ZRHO_s_cell = 0.25*(ZRHO_sC[i][j]+ZRHO_sC[ip1][j]+ZRHO_sC[i][jp1]+ZRHO_sC[ip1][jp1]);
								ZRHO_s_cell = RHO_s_cell*C.Z_sC[i][j];
								RHO_s_cell  -=tau*(stag_RHO_F_sx[i][j+1] -stag_RHO_F_sx[i][j])/dx;
								ZRHO_s_cell -=tau*(stag_ZRHO_F_sx[i][j+1]-stag_ZRHO_F_sx[i][j])/dx;
								C.Z_sC[i][j]   =ZRHO_s_cell/RHO_s_cell;
							}
					for(i = 0; i < n_y; ++i)
						for(j = 0; j < n_x; ++j)
							{
								i_1=i-1>0?i-1:0;
								C.Z_sS_xd[i][j] = 0.5*(C.Z_sC[i_1][j]+C.Z_sC[i][j]);
							}
					for(i = 0; i < n_y; ++i)
						for(j = 0; j < n_x; ++j)
							{
								j_1=j-1>0?j-1:0;
								if(fabs(U_R[i][j].z_s-U_L[i][j].z_s)<eps)
									S = C.P_gC[i][j]*(U_R[i][j].z_s-U_L[i][j].z_s);
								else
									S = U_R[i][j].z_s*FX.P_s_MIDx[i][j+1]-U_L[i][j].z_s*FX.P_s_MIDx[i][j];
								C.ZRHO_gC[i][j] -= tau*(FX.RHO_F_gx[i][j+1]-FX.RHO_F_gx[i][j])/dx;
								C.RHO_U_gC[i][j]-= tau*(FX.U_F_gx[i][j+1]  -FX.U_F_gx[i][j])  /dx-tau/dx*S;
								C.RHO_V_gC[i][j]-= tau*(FX.V_F_gx[i][j+1]  -FX.V_F_gx[i][j])  /dx;
								C.E_gC[i][j]    -= tau*(FX.E_F_gx[i][j+1]  -FX.E_F_gx[i][j])  /dx-tau/dx*S*C.U_sC[i][j];
								C.ZRHO_sC[i][j] -= tau*(FX.RHO_F_sx[i][j+1]-FX.RHO_F_sx[i][j])/dx;
								C.RHO_U_sC[i][j]-= tau*(FX.U_F_sx[i][j+1]  -FX.U_F_sx[i][j])  /dx+tau/dx*S;
								C.RHO_V_sC[i][j]-= tau*(FX.V_F_sx[i][j+1]  -FX.V_F_sx[i][j])  /dx;
								C.E_sC[i][j]    -= tau*(FX.E_F_sx[i][j+1]  -FX.E_F_sx[i][j])  /dx+tau/dx*S*C.U_sC[i][j];
								area_L=0.5+C.U_sC[i][j]*tau/dx;
								area_R=1.0-area_L;
								BN_C2U(C,U,i,j,0);
								primitive_comp(U, &U_L[i][j], &U_R[i][j], U_L[i][j].z_s, U_R[i][j].z_s, C.Z_sS_xd[i][j_1], C.Z_sS_xd[i][j], area_L, area_R);
								BN_ULR2cons(U_L,U_R,C,i,j);
							}
				}
			else if (stop_step == 2)
				{
					for(j = 0; j <  n_x; ++j)
						for(i = 0; i <= n_y; ++i)
							{
								iL = i-1>0?i-1:0;
								iR = i<n_y?i:n_y-1;
								z_smid = U_R[iL][j].z_s;
								z_gmid = 1.0-z_smid;
								z_sy_mid = SV.Z_sS_y[i][j];
								GRP_var_init(&GL, SV, U_R, dy, iL, j, 2);
								GRP_var_init(&GR, SV, U_L, dy, iR, j, 3);
								GRP_RI_var_init(&RI_L, SV, C, dx, iL, j, 2);
								GRP_RI_var_init(&RI_R, SV, C, dx, iR, j, 3);
								RI_LR2G_LR(&RI_L,&GL,z_smid,1);
								RI_LR2G_LR(&RI_R,&GR,z_smid,1);
								linear_GRP_solver_Edir_Q1D(wave_speed, dire, mid_g, star, 0.0, 0.0, GL.rho_g, GR.rho_g, GL.rho_gy, GR.rho_gy, -GL.rho_gx, -GR.rho_gx, GL.v_g, GR.v_g, GL.v_gy, GR.v_gy, -GL.v_gx, -GR.v_gx, -GL.u_g, -GR.u_g, -GL.u_gy, -GR.u_gy, GL.u_gx, GR.u_gx, GL.p_g, GR.p_g, GL.p_gy, GR.p_gy, -GL.p_gx, -GR.p_gx, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, gamma_g, gamma_g, eps, eps);
								rho_gmid = mid_g[0] + 0.5*tau*dire[0];
								u_gmid   =-mid_g[2] - 0.5*tau*dire[2];
								v_gmid   = mid_g[1] + 0.5*tau*dire[1];
								p_gmid   = mid_g[3] + 0.5*tau*dire[3];
								linear_GRP_solver_Edir_Q1D(wave_speed, dire, mid_s, star, 0.0, 0.0, GL.rho_s, GR.rho_s, GL.rho_sy, GR.rho_sy, -GL.rho_sx, -GR.rho_sx, GL.v_s, GR.v_s, GL.v_sy, GR.v_sy, -GL.v_sx, -GR.v_sx, -GL.u_s, -GR.u_s, -GL.u_sy, -GR.u_sy, GL.u_sx, GR.u_sx, GL.p_s, GR.p_s, GL.p_sy, GR.p_sy, -GL.p_sx, -GR.p_sx, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, gamma_s, gamma_s, eps, eps);
								rho_smid = mid_s[0] + 0.5*tau*dire[0];
								u_smid   =-mid_s[2] - 0.5*tau*dire[2];
								v_smid   = mid_s[1] + 0.5*tau*dire[1];
								p_smid   = mid_s[3] + 0.5*tau*dire[3];
								linear_GRP_RI_solver_BN(&RI, z_sy_mid, z_smid, mid_g, mid_s, RI_L, RI_R, GL, GR, gamma_s, gamma_g, eps, tau, 1);
								RI2U_cal(&U_tmp, &RI, RI.z_s, rho_gmid);
								rho_gmid = U_tmp.rho_g;
								u_gmid   =-mid_g[2];
								v_gmid   = U_tmp.u_g;
								p_gmid   = U_tmp.p_g;
								rho_smid = U_tmp.rho_s;
								u_smid   =-mid_s[2];
								v_smid   = U_tmp.u_s;
								p_smid   = U_tmp.p_s;

								FX.RHO_F_gy[i][j] = z_gmid*rho_gmid*v_gmid;
								FX.U_F_gy[i][j]   = FX.RHO_F_gy[i][j]*u_gmid;
								FX.V_F_gy[i][j]   = FX.RHO_F_gy[i][j]*v_gmid + z_gmid*p_gmid;
								FX.E_F_gy[i][j]   = gamma_g/(gamma_g-1.0)*p_gmid/rho_gmid + 0.5*(u_gmid*u_gmid + v_gmid*v_gmid);
								FX.E_F_gy[i][j]   = FX.RHO_F_gy[i][j]*FX.E_F_gy[i][j];
								FX.RHO_F_sy[i][j] = z_smid*rho_smid*v_smid;
								FX.U_F_sy[i][j]   = FX.RHO_F_sy[i][j]*u_smid;
								FX.V_F_sy[i][j]   = FX.RHO_F_sy[i][j]*v_smid + z_smid*p_smid;
								FX.E_F_sy[i][j]   = gamma_s/(gamma_s-1.0)*p_smid/rho_smid + 0.5*(u_smid*u_smid + v_smid*v_smid);
								FX.E_F_sy[i][j]   = FX.RHO_F_sy[i][j]*FX.E_F_sy[i][j];
								FX.P_s_MIDy[i][j] = p_smid;
								FX.Z_s_MIDy[i][j] = z_smid;
							}
					for(j = 0; j <  n_x; ++j)
						for(i = 0; i <= n_y; ++i)
							{
								iL = i-1>0?i-1:0;
								iR = i<n_y?i:n_y-1;
								z_sxL = SV.Z_sx[iL][j];
								z_sxR = SV.Z_sx[iR][j];
								z_syL = SV.Z_sy[iL][j];
								z_syR = SV.Z_sy[iR][j];
								z_sL = C.Z_sC[iL][j]+dy/2*z_syL;
								z_sR = C.Z_sC[iR][j]-dy/2*z_syR;
								jp1 = j<n_x-1?j+1:n_x-1;
								stag_RHO_F_sy[i][j] = 0.5*(U_R[i][j].rho_s*U_R[i][j].u_s+U_R[i][jp1].rho_s*U_R[i][jp1].u_s);
								if ((U_R[i][j].u_s + U_R[i][jp1].u_s) > 0.0)
									stag_ZRHO_F_sy[i][j] = stag_RHO_F_sy[i][j]*(z_sL - 0.25*tau*(-(U_R[i][j].v_s + U_R[i][jp1].v_s)*z_sxL+(U_R[i][j].u_s + U_R[i][jp1].u_s)*z_syL));
								else
									stag_ZRHO_F_sy[i][j] = stag_RHO_F_sy[i][j]*(z_sR - 0.25*tau*(-(U_R[i][j].v_s + U_R[i][jp1].v_s)*z_sxR+(U_R[i][j].u_s + U_R[i][jp1].u_s)*z_syR));
							}
					for(j = 0; j < n_x; ++j)
						for(i = 0; i < n_y; ++i)
							{
								ip1 = i<n_y-1?i+1:n_y-1;
								jp1 = j<n_x-1?j+1:n_x-1;
								RHO_s_cell  = 0.25*(C.RHO_sC[i][j]+C.RHO_sC[ip1][j]+C.RHO_sC[i][jp1]+C.RHO_sC[ip1][jp1]);
								//ZRHO_s_cell = 0.25*(ZRHO_sC[i][j]+ZRHO_sC[ip1][j]+ZRHO_sC[i][jp1]+ZRHO_sC[ip1][jp1]);
								ZRHO_s_cell = RHO_s_cell*C.Z_sC[i][j];
								RHO_s_cell  -=tau*(stag_RHO_F_sy[i+1][j] -stag_RHO_F_sy[i][j])/dy;
								ZRHO_s_cell -=tau*(stag_ZRHO_F_sy[i+1][j]-stag_ZRHO_F_sy[i][j])/dy;
								C.Z_sC[i][j]   =ZRHO_s_cell/RHO_s_cell;
							}
					for(j = 0; j < n_x; ++j)
						for(i = 0; i < n_y; ++i)
							{
								j_1=j-1>0?j-1:0;
								C.Z_sS_yd[i][j] = 0.5*(C.Z_sC[i][j_1]+C.Z_sC[i][j]);
							}
					for(j = 0; j < n_x; ++j)
						for(i = 0; i < n_y; ++i)
							{
								i_1=i-1>0?i-1:0;
								if(fabs(U_R[i][j].z_s-U_L[i][j].z_s)<eps)
									S = C.P_gC[i][j]*(U_R[i][j].z_s-U_L[i][j].z_s);
								else
									S = U_R[i][j].z_s*FX.P_s_MIDy[i+1][j]-U_L[i][j].z_s*FX.P_s_MIDy[i][j];
								C.ZRHO_gC[i][j] -= tau*(FX.RHO_F_gy[i+1][j]-FX.RHO_F_gy[i][j])/dy;
								C.RHO_U_gC[i][j]-= tau*(FX.U_F_gy[i+1][j]  -FX.U_F_gy[i][j])  /dy;
								C.RHO_V_gC[i][j]-= tau*(FX.V_F_gy[i+1][j]  -FX.V_F_gy[i][j])  /dy-tau/dx*S;
								C.E_gC[i][j]    -= tau*(FX.E_F_gy[i+1][j]  -FX.E_F_gy[i][j])  /dy-tau/dx*S*C.V_sC[i][j];
								C.ZRHO_sC[i][j] -= tau*(FX.RHO_F_sy[i+1][j]-FX.RHO_F_sy[i][j])/dy;
								C.RHO_U_sC[i][j]-= tau*(FX.U_F_sy[i+1][j]  -FX.U_F_sy[i][j])  /dy;
								C.RHO_V_sC[i][j]-= tau*(FX.V_F_sy[i+1][j]  -FX.V_F_sy[i][j])  /dy+tau/dx*S;
								C.E_sC[i][j]    -= tau*(FX.E_F_sy[i+1][j]  -FX.E_F_sy[i][j])  /dy+tau/dx*S*C.V_sC[i][j];
								area_L=0.5+C.V_sC[i][j]*tau/dx;
								area_R=1.0-area_L;
								BN_C2U(C,U,i,j,1);
								primitive_comp(U, &U_L[i][j], &U_R[i][j], U_L[i][j].z_s, U_R[i][j].z_s, C.Z_sS_yd[i_1][j], C.Z_sS_yd[i][j], area_L, area_R);
								BN_ULR2cons(U_L,U_R,C,i,j);
							}
				}
			for(j = 0; j < n_x; ++j)
				{
					C.Z_sC[n_y-1][j]     = C.Z_sC[n_y-2][j];
					C.ZRHO_gC[n_y-1][j]  = C.ZRHO_gC[n_y-2][j];
					C.RHO_U_gC[n_y-1][j] = C.RHO_U_gC[n_y-2][j];
					C.RHO_V_gC[n_y-1][j] = C.RHO_V_gC[n_y-2][j];
					C.E_gC[n_y-1][j]     = C.E_gC[n_y-2][j];
					C.ZRHO_sC[n_y-1][j]  = C.ZRHO_sC[n_y-2][j];
					C.RHO_U_sC[n_y-1][j] = C.RHO_U_sC[n_y-2][j];
					C.RHO_V_sC[n_y-1][j] = C.RHO_V_sC[n_y-2][j];
					C.E_sC[n_y-1][j]     = C.E_sC[n_y-2][j];
				}
			for(i = 0; i < n_y; ++i)
				{
					C.Z_sC[i][n_x-1]     = C.Z_sC[i][n_x-2];
					C.ZRHO_gC[i][n_x-1]  = C.ZRHO_gC[i][n_x-2];
					C.RHO_U_gC[i][n_x-1] = C.RHO_U_gC[i][n_x-2];
					C.RHO_V_gC[i][n_x-1] = C.RHO_V_gC[i][n_x-2];
					C.E_gC[i][n_x-1]     = C.E_gC[i][n_x-2];
					C.ZRHO_sC[i][n_x-1]  = C.ZRHO_sC[i][n_x-2];
					C.RHO_U_sC[i][n_x-1] = C.RHO_U_sC[i][n_x-2];
					C.RHO_V_sC[i][n_x-1] = C.RHO_V_sC[i][n_x-2];
					C.E_sC[i][n_x-1]     = C.E_sC[i][n_x-2];
				}
			if (stop_step == 0)
				stop_step = 2;
			else if (stop_step == 2)
				stop_step = 0;
			if (stop_t == 1)
				break;

			DispPro(t_all*100.0/config[1], l);
			cpu_time += (clock() - start_clock) / (double)CLOCKS_PER_SEC;
		}

	for(i = 0; i < n_y; ++i)
		for(j = 0; j < n_x; ++j)
			{
				FV->RHO[i*n_x+j]=C.RHO_gC[i][j];
				FV->U[i*n_x+j]  =C.U_gC[i][j];
				FV->V[i*n_x+j]  =C.V_gC[i][j];
				FV->P[i*n_x+j]  =C.P_gC[i][j];
				FV->Z_a[i*n_x+j]=C.Z_sC[i][j];
			}

	printf("\nThe cost of CPU time for the Eulerian method is %g seconds.\n", cpu_time);
}
