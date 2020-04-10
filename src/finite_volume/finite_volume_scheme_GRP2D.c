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

				if (order == 1)
					for(k=0, p=SV.Z_sx; k<sizeof(struct slope_var)/sizeof(double **); k++)
						(p++)[i][j] = 0.0;
			}
	/*
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
	*/
}

void finite_volume_scheme_GRP2D(struct flu_var * FV, const struct mesh_var * mv, const char * scheme, const char * problem) //R-K
{
	clock_t start_clock;
	double cpu_time = 0.0;

	const int dim = (int)config[0];
	const int order = (int)config[9];
	const double eps = config[4];
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
	int i, j, k, stop_step = 0, stop_t = 0;
	double eps_big = 1e-10, eps_big2 = 1e-3, eps3=1e-10;

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

	double z_sL, z_sR, z_sxL, z_sxR, z_syL, z_syR;
	double rho_gL, rho_gR, p_gL ,p_gR ,u_gL ,u_gR ,v_gL ,v_gR;
	double rho_gxL, rho_gxR, p_gxL,p_gxR,u_gxL,u_gxR,v_gxL,v_gxR;
	double rho_gyL, rho_gyR, p_gyL,p_gyR,u_gyL,u_gyR,v_gyL,v_gyR;
	double rho_sL, rho_sR, p_sL ,p_sR ,u_sL ,u_sR ,v_sL ,v_sR;
	double rho_sxL, rho_sxR, p_sxL,p_sxR,u_sxL,u_sxR,v_sxL,v_sxR;
	double rho_syL, rho_syR, p_syL,p_syR,u_syL,u_syR,v_syL,v_syR;

	double z_gmid, rho_gmid, p_gmid, u_gmid, v_gmid;
	double z_smid, rho_smid, p_smid, u_smid, v_smid;

	int iL,iR,jL,jR, i_1, j_1, ip1, jp1;
	double a_g, a_s, S_max, S_max_g, S_max_s;
	double wave_speed[2], dire[6], mid[6], star[6];

	double U[8];
	struct U_var U_L[n_y][n_x], U_R[n_y][n_x];
	double Z_sL[n_y][n_x], Z_sR[n_y][n_x];
	struct RI_var RI;
	double S, S_tmp, area_L, area_R, RHO_s_cell, ZRHO_s_cell;
	int l;
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
								if (i == 0)
									i_1 = 0;
								else
									i_1 = i-1;
								if (j == 0)
									j_1 = 0;
								else
									j_1 = j-1;
								U[0] = C.ZRHO_gC[i][j];
								U[1] = C.RHO_U_gC[i][j];
								U[2] = C.RHO_V_gC[i][j];
								U[3] = C.E_gC[i][j]-0.5*pow(C.RHO_V_gC[i][j],2)/C.ZRHO_gC[i][j];
								U[4] = C.ZRHO_sC[i][j];
								U[5] = C.RHO_U_sC[i][j];
								U[6] = C.RHO_V_sC[i][j];
								U[7] = C.E_sC[i][j]-0.5*pow(C.RHO_V_sC[i][j],2)/C.ZRHO_sC[i][j];
								Z_sL[i][j] = 0.5*(C.Z_sC[i_1][j_1]+C.Z_sC[i][j_1]);
								Z_sR[i][j] = 0.5*(C.Z_sC[i_1][j]  +C.Z_sC[i][j]);
								primitive_comp(U, &U_L[i][j], &U_R[i][j], Z_sL[i][j], Z_sR[i][j], Z_sL[i][j], Z_sR[i][j], 0.5, 0.5);
								C.RHO_gC[i][j] = 0.5*(U_L[i][j].rho_g+U_R[i][j].rho_g);
								C.U_gC[i][j]   = 0.5*(U_L[i][j].u_g +U_R[i][j].u_g);
								C.V_gC[i][j]   = 0.5*(U_L[i][j].v_g +U_R[i][j].v_g);
								C.P_gC[i][j]   = 0.5*(U_L[i][j].p_g +U_R[i][j].p_g);
								C.RHO_sC[i][j] = 0.5*(U_L[i][j].rho_s+U_R[i][j].rho_s);
								C.U_sC[i][j]   = 0.5*(U_L[i][j].u_s +U_R[i][j].u_s);
								C.V_sC[i][j]   = 0.5*(U_L[i][j].v_s +U_R[i][j].v_s);
								C.P_sC[i][j]   = 0.5*(U_L[i][j].p_s +U_R[i][j].p_s);
							}
				}
			else if (stop_step == 2)
				{
					for(j = 0; j < n_x; ++j)
						for(i = 0; i < n_y; ++i)
							{
								if (i == 0)
									i_1 = 0;
								else
									i_1 = i-1;
								if (j == 0)
									j_1 = 0;
								else
									j_1 = j-1;
								U[0] = C.ZRHO_gC[i][j];
								U[1] = C.RHO_V_gC[i][j];
								U[2] = C.RHO_U_gC[i][j];
								U[3] = C.E_gC[i][j]-0.5*pow(C.RHO_U_gC[i][j],2)/C.ZRHO_gC[i][j];
								U[4] = C.ZRHO_sC[i][j];
								U[5] = C.RHO_V_sC[i][j];
								U[6] = C.RHO_U_sC[i][j];
								U[7] = C.E_sC[i][j]-0.5*pow(C.RHO_U_sC[i][j],2)/C.ZRHO_sC[i][j];
								Z_sL[i][j] = 0.5*(C.Z_sC[i_1][j_1]+C.Z_sC[i_1][j]);
								Z_sR[i][j] = 0.5*(C.Z_sC[i][j_1]+C.Z_sC[i][j]);
								primitive_comp(U, &U_L[i][j], &U_R[i][j], Z_sL[i][j], Z_sR[i][j], Z_sL[i][j], Z_sR[i][j], 0.5, 0.5);
								C.RHO_gC[i][j] = 0.5*(U_L[i][j].rho_g+U_R[i][j].rho_g);
								C.U_gC[i][j]   = 0.5*(U_L[i][j].v_g +U_R[i][j].v_g);
								C.V_gC[i][j]   = 0.5*(U_L[i][j].u_g +U_R[i][j].u_g);
								C.P_gC[i][j]   = 0.5*(U_L[i][j].p_g +U_R[i][j].p_g);
								C.RHO_sC[i][j] = 0.5*(U_L[i][j].rho_s+U_R[i][j].rho_s);
								C.U_sC[i][j]   = 0.5*(U_L[i][j].v_s +U_R[i][j].v_s);
								C.V_sC[i][j]   = 0.5*(U_L[i][j].u_s +U_R[i][j].u_s);
								C.P_sC[i][j]   = 0.5*(U_L[i][j].p_s +U_R[i][j].p_s);
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
								if (j==0)
									{
										jL=0;
										jR=0;
									}
								else if (j==n_x)
									{
										jL=n_x-1;
										jR=n_x-1;
									}
								else
									{
										jL=j-1;
										jR=j;
									}
								z_smid = U_R[i][jL].z_s;
								z_gmid = 1.0-z_smid;
								Z_s_MIDx[i][j] = z_smid;

								rho_gxL=RHO_gx[i][jL];
								rho_gxR=RHO_gx[i][jR];
								p_gxL=P_gx[i][jL];
								p_gxR=P_gx[i][jR];
								u_gxL=U_gx[i][jL];
								u_gxR=U_gx[i][jR];
								v_gxL=V_gx[i][jL];
								v_gxR=V_gx[i][jR];
								rho_gyL=RHO_gy[i][jL];
								rho_gyR=RHO_gy[i][jR];
								p_gyL=P_gy[i][jL];
								p_gyR=P_gy[i][jR];
								u_gyL=U_gy[i][jL];
								u_gyR=U_gy[i][jR];
								v_gyL=V_gy[i][jL];
								v_gyR=V_gy[i][jR];
								rho_gL =U_R[i][jL].rho_g+dx/2*rho_gxL;
								rho_gR =U_L[i][jR].rho_g-dx/2*rho_gxR;
								p_gL =U_R[i][jL].p_g+dx/2*p_gxL;
								p_gR =U_L[i][jR].p_g-dx/2*p_gxR;
								u_gL =U_R[i][jL].u_g+dx/2*u_gxL;
								u_gR =U_L[i][jR].u_g-dx/2*u_gxR;
								v_gL =U_R[i][jL].v_g+dx/2*v_gxL;
								v_gR =U_L[i][jR].v_g-dx/2*v_gxR;

								rho_sxL=RHO_sx[i][jL];
								rho_sxR=RHO_sx[i][jR];
								p_sxL=P_sx[i][jL];
								p_sxR=P_sx[i][jR];
								u_sxL=U_sx[i][jL];
								u_sxR=U_sx[i][jR];
								v_sxL=V_sx[i][jL];
								v_sxR=V_sx[i][jR];
								rho_syL=RHO_sy[i][jL];
								rho_syR=RHO_sy[i][jR];
								p_syL=P_sy[i][jL];
								p_syR=P_sy[i][jR];
								u_syL=U_sy[i][jL];
								u_syR=U_sy[i][jR];
								v_syL=V_sy[i][jL];
								v_syR=V_sy[i][jR];
								rho_sL =U_R[i][jL].rho_s+dx/2*rho_sxL;
								rho_sR =U_L[i][jR].rho_s-dx/2*rho_sxR;
								p_sL =U_R[i][jL].p_s+dx/2*p_sxL;
								p_sR =U_L[i][jR].p_s-dx/2*p_sxR;
								u_sL =U_R[i][jL].u_s+dx/2*u_sxL;
								u_sR =U_L[i][jR].u_s-dx/2*u_sxR;
								v_sL =U_R[i][jL].v_s+dx/2*v_sxL;
								v_sR =U_L[i][jR].v_s-dx/2*v_sxR;

								linear_GRP_solver_Edir_Q1D(wave_speed, dire, mid, star, 0.0, 0.0, rho_gL, rho_gR, rho_gxL, rho_gxR, rho_gyL, rho_gyR, u_gL, u_gR, u_gxL, u_gxR, u_gyL, u_gyR, v_gL, v_gR, v_gxL, v_gxR, v_gyL, v_gyR, p_gL, p_gR, p_gxL, p_gxR, p_gyL, p_gyR, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, gamma_g, gamma_g, eps, eps);
								rho_gmid = mid[0] + 0.5*tau*dire[0];
								u_gmid   = mid[1] + 0.5*tau*dire[1];
								v_gmid   = mid[2] + 0.5*tau*dire[2];
								p_gmid   = mid[3] + 0.5*tau*dire[3];
								RHO_F_gx[i][j] = z_gmid*rho_gmid*u_gmid;
								U_F_gx[i][j]   = RHO_F_gx[i][j]*u_gmid + z_gmid*p_gmid;
								V_F_gx[i][j]   = RHO_F_gx[i][j]*v_gmid;
								E_F_gx[i][j]   = gamma_g/(gamma_g-1.0)*p_gmid/rho_gmid + 0.5*(u_gmid*u_gmid + v_gmid*v_gmid);
								E_F_gx[i][j]   = RHO_F_gx[i][j]*E_F_gx[i][j];

								linear_GRP_solver_Edir_Q1D(wave_speed, dire, mid, star, 0.0, 0.0, rho_sL, rho_sR, rho_sxL, rho_sxR, rho_syL, rho_syR, u_sL, u_sR, u_sxL, u_sxR, u_syL, u_syR, v_sL, v_sR, v_sxL, v_sxR, v_syL, v_syR, p_sL, p_sR, p_sxL, p_sxR, p_syL, p_syR, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, gamma_s, gamma_s, eps, eps);
								rho_smid = mid[0] + 0.5*tau*dire[0];
								u_smid   = mid[1] + 0.5*tau*dire[1];
								v_smid   = mid[2] + 0.5*tau*dire[2];
								p_smid   = mid[3] + 0.5*tau*dire[3];
								P_s_MIDx[i][j] = p_smid;
								RHO_F_sx[i][j] = z_smid*rho_smid*u_smid;
								U_F_sx[i][j]   = RHO_F_sx[i][j]*u_smid + z_smid*p_smid;
								V_F_sx[i][j]   = RHO_F_sx[i][j]*v_smid;
								E_F_sx[i][j]   = gamma_s/(gamma_s-1.0)*p_smid/rho_smid + 0.5*(u_smid*u_smid + v_smid*v_smid);
								E_F_sx[i][j]   = RHO_F_sx[i][j]*E_F_sx[i][j];
							}
					for(i = 0; i < n_y; ++i)
						for(j = 0; j <=n_x; ++j)
							{
								if (j==0)
									{
										jL=0;
										jR=0;
									}
								else if (j==n_x)
									{
										jL=n_x-1;
										jR=n_x-1;
									}
								else
									{
										jL=j-1;
										jR=j;
									}
								z_sxL = Z_sx[i][jL];
								z_sxR = Z_sx[i][jR];
								z_syL = Z_sy[i][jL];
								z_syR = Z_sy[i][jR];
								z_sL =Z_sC[i][jL]+dx/2*z_sxL;
								z_sR =Z_sC[i][jR]-dx/2*z_sxR;
								if (i == n_y-1)
									ip1 = n_y-1;
								else
									ip1 = i+1;
								stag_RHO_F_sx[i][j] = 0.5*(U_R[i][j].rho_s*U_R[i][j].u_s+U_R[ip1][j].rho_s*U_R[ip1][j].u_s);
								if ((U_R[i][j].u_s + U_R[ip1][j].u_s) > 0.0)
									stag_ZRHO_F_sx[i][j] = stag_RHO_F_sx[i][j]*(z_sL - 0.25*tau*((U_R[i][j].u_s + U_R[ip1][j].u_s)*z_sxL+(U_R[i][j].v_s + U_R[ip1][j].v_s)*z_syL));
								else
									stag_ZRHO_F_sx[i][j] = stag_RHO_F_sx[i][j]*(z_sR - 0.25*tau*((U_R[i][j].u_s + U_R[ip1][j].u_s)*z_sxR+(U_R[i][j].v_s + U_R[ip1][j].v_s)*z_syR));
							}
					for(i = 0; i < n_y; ++i)
						for(j = 0; j < n_x; ++j)
							{
								if (i == n_y-1)
									ip1 = n_y-1;
								else
									ip1 = i+1;
								if (j == n_x-1)
									jp1 = n_x-1;
								else
									jp1 = j+1;
								RHO_s_cell  = 0.25*(RHO_sC[i][j]+RHO_sC[ip1][j]+RHO_sC[i][jp1]+RHO_sC[ip1][jp1]);
								//ZRHO_s_cell = 0.25*(ZRHO_sC[i][j]+ZRHO_sC[ip1][j]+ZRHO_sC[i][jp1]+ZRHO_sC[ip1][jp1]);
								ZRHO_s_cell = RHO_s_cell*Z_sC[i][j];
								RHO_s_cell  -=tau*(stag_RHO_F_sx[i][j+1] -stag_RHO_F_sx[i][j])/dx;
								ZRHO_s_cell -=tau*(stag_ZRHO_F_sx[i][j+1]-stag_ZRHO_F_sx[i][j])/dx;
								Z_sC[i][j]   =ZRHO_s_cell/RHO_s_cell;
							}
					for(i = 0; i < n_y; ++i)
						for(j = 0; j < n_x; ++j)
							{
								if (i == 0)
									i_1 = 0;
								else
									i_1 = i-1;
								if (j == 0)
									j_1 = 0;
								else
									j_1 = j-1;
								Z_sL[i][j] = 0.5*(Z_sC[i_1][j_1]+Z_sC[i][j_1]);
								Z_sR[i][j] = 0.5*(Z_sC[i_1][j]+Z_sC[i][j]);
							}
					for(i = 0; i < n_y; ++i)
						for(j = 0; j < n_x; ++j)
							{
								if(fabs(U_R[i][j].z_s-U_L[i][j].z_s)<eps)
									S = P_gC[i][j]*(U_R[i][j].z_s-U_L[i][j].z_s);
								else
									S = U_R[i][j].z_s*P_s_MIDx[i][j+1]-U_L[i][j].z_s*P_s_MIDx[i][j];
								ZRHO_gC[i][j] -= tau*(RHO_F_gx[i][j+1]-RHO_F_gx[i][j])/dx;
								RHO_U_gC[i][j]-= tau*(U_F_gx[i][j+1]  -U_F_gx[i][j])  /dx-tau/dx*S;
								RHO_V_gC[i][j]-= tau*(V_F_gx[i][j+1]  -V_F_gx[i][j])  /dx;
								E_gC[i][j]    -= tau*(E_F_gx[i][j+1]  -E_F_gx[i][j])  /dx-tau/dx*S*U_sC[i][j];
								ZRHO_sC[i][j] -= tau*(RHO_F_sx[i][j+1]-RHO_F_sx[i][j])/dx;
								RHO_U_sC[i][j]-= tau*(U_F_sx[i][j+1]  -U_F_sx[i][j])  /dx+tau/dx*S;
								RHO_V_sC[i][j]-= tau*(V_F_sx[i][j+1]  -V_F_sx[i][j])  /dx;
								E_sC[i][j]    -= tau*(E_F_sx[i][j+1]  -E_F_sx[i][j])  /dx+tau/dx*S*U_sC[i][j];
								area_L=0.5+U_sC[i][j]*tau/dx;
								area_R=1.0-area_L;
								U[0] = ZRHO_gC[i][j];
								U[1] = RHO_U_gC[i][j];
								U[2] = RHO_V_gC[i][j];
								U[3] = E_gC[i][j]-0.5*RHO_V_gC[i][j]*RHO_V_gC[i][j]/ZRHO_gC[i][j];
								U[4] = ZRHO_sC[i][j];
								U[5] = RHO_U_sC[i][j];
								U[6] = RHO_V_sC[i][j];
								U[7] = E_sC[i][j]-0.5*RHO_V_sC[i][j]*RHO_V_sC[i][j]/ZRHO_sC[i][j];
								primitive_comp(U, &U_L[i][j], &U_R[i][j], U_L[i][j].z_s, U_R[i][j].z_s, Z_sL[i][j], Z_sR[i][j], area_L, area_R);
								ZRHO_gC[i][j]  = 0.5*(U_L[i][j].U_rho_g+U_R[i][j].U_rho_g);
								RHO_U_gC[i][j] = 0.5*(U_L[i][j].U_u_g +U_R[i][j].U_u_g);
								RHO_V_gC[i][j] = 0.5*(U_L[i][j].U_v_g +U_R[i][j].U_v_g);
								E_gC[i][j]     = 0.5*(U_L[i][j].U_e_g +U_R[i][j].U_e_g);
								ZRHO_sC[i][j]  = 0.5*(U_L[i][j].U_rho_s+U_R[i][j].U_rho_s);
								RHO_U_sC[i][j] = 0.5*(U_L[i][j].U_u_s +U_R[i][j].U_u_s);
								RHO_V_sC[i][j] = 0.5*(U_L[i][j].U_v_s +U_R[i][j].U_v_s);
								E_sC[i][j]     = 0.5*(U_L[i][j].U_e_s +U_R[i][j].U_e_s);
							}
				}
			else if (stop_step == 2)
				{
					for(j = 0; j <  n_x; ++j)
						for(i = 0; i <= n_y; ++i)
							{
								if (i==0)
									{
										iL=0;
										iR=0;
									}
								else if (i==n_y)
									{
										iL=n_y-1;
										iR=n_y-1;
									}
								else
									{
										iL=i-1;
										iR=i;
									}
								z_smid = U_R[iL][j].z_s;
								z_gmid = 1.0-z_smid;
								Z_s_MIDy[i][j] = z_smid;

								rho_gxL=RHO_gx[iL][j];
								rho_gxR=RHO_gx[iR][j];
								p_gxL=P_gx[iL][j];
								p_gxR=P_gx[iR][j];
								u_gxL=U_gx[iL][j];
								u_gxR=U_gx[iR][j];
								v_gxL=V_gx[iL][j];
								v_gxR=V_gx[iR][j];
								rho_gyL=RHO_gy[iL][j];
								rho_gyR=RHO_gy[iR][j];
								p_gyL=P_gy[iL][j];
								p_gyR=P_gy[iR][j];
								u_gyL=U_gy[iL][j];
								u_gyR=U_gy[iR][j];
								v_gyL=V_gy[iL][j];
								v_gyR=V_gy[iR][j];
								rho_gL =U_R[iL][j].rho_g+dy/2*rho_gyL;
								rho_gR =U_L[iR][j].rho_g-dy/2*rho_gyR;
								p_gL =U_R[iL][j].p_g+dy/2*p_gyL;
								p_gR =U_L[iR][j].p_g-dy/2*p_gyR;
								u_gL =U_R[iL][j].v_g+dy/2*u_gyL;
								u_gR =U_L[iR][j].v_g-dy/2*u_gyR;
								v_gL =U_R[iL][j].u_g+dy/2*v_gyL;
								v_gR =U_L[iR][j].u_g-dy/2*v_gyR;

								rho_sxL=RHO_sx[iL][j];
								rho_sxR=RHO_sx[iR][j];
								p_sxL=P_sx[iL][j];
								p_sxR=P_sx[iR][j];
								u_sxL=U_sx[iL][j];
								u_sxR=U_sx[iR][j];
								v_sxL=V_sx[iL][j];
								v_sxR=V_sx[iR][j];
								rho_syL=RHO_sy[iL][j];
								rho_syR=RHO_sy[iR][j];
								p_syL=P_sy[iL][j];
								p_syR=P_sy[iR][j];
								u_syL=U_sy[iL][j];
								u_syR=U_sy[iR][j];
								v_syL=V_sy[iL][j];
								v_syR=V_sy[iR][j];
								rho_sL =U_R[iL][j].rho_s+dy/2*rho_syL;
								rho_sR =U_L[iR][j].rho_s-dy/2*rho_syR;
								p_sL =U_R[iL][j].p_s+dy/2*p_syL;
								p_sR =U_L[iR][j].p_s-dy/2*p_syR;
								u_sL =U_R[iL][j].v_s+dy/2*u_syL;
								u_sR =U_L[iR][j].v_s-dy/2*u_syR;
								v_sL =U_R[iL][j].u_s+dy/2*v_syL;
								v_sR =U_L[iR][j].u_s-dy/2*v_syR;

								linear_GRP_solver_Edir_Q1D(wave_speed, dire, mid, star, 0.0, 0.0, rho_gL, rho_gR, rho_gyL, rho_gyR, -rho_gxL, -rho_gxR, v_gL, v_gR, v_gyL, v_gyR, -v_gxL, -v_gxR, -u_gL, -u_gR, -u_gyL, -u_gyR, u_gxL, u_gxR, p_gL, p_gR, p_gyL, p_gyR, -p_gxL, -p_gxR, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, gamma_g, gamma_g, eps, eps);
								rho_gmid = mid[0] + 0.5*tau*dire[0];
								u_gmid   =-mid[2] - 0.5*tau*dire[2];
								v_gmid   = mid[1] + 0.5*tau*dire[1];
								p_gmid   = mid[3] + 0.5*tau*dire[3];
								RHO_F_gy[i][j] = z_gmid*rho_gmid*v_gmid;
								U_F_gy[i][j]   = RHO_F_gy[i][j]*u_gmid;
								V_F_gy[i][j]   = RHO_F_gy[i][j]*v_gmid + z_gmid*p_gmid;
								E_F_gy[i][j]   = gamma_g/(gamma_g-1.0)*p_gmid/rho_gmid + 0.5*(u_gmid*u_gmid + v_gmid*v_gmid);
								E_F_gy[i][j]   = RHO_F_gy[i][j]*E_F_gy[i][j];

								linear_GRP_solver_Edir_Q1D(wave_speed, dire, mid, star, 0.0, 0.0, rho_sL, rho_sR, rho_syL, rho_syR, -rho_sxL, -rho_sxR, v_sL, v_sR, v_syL, v_syR, -v_sxL, -v_sxR, -u_sL, -u_sR, -u_syL, -u_syR, u_sxL, u_sxR, p_sL, p_sR, p_syL, p_syR, -p_sxL, -p_sxR, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, gamma_s, gamma_s, eps, eps);
								rho_smid = mid[0] + 0.5*tau*dire[0];
								u_smid   =-mid[2] - 0.5*tau*dire[2];
								v_smid   = mid[1] + 0.5*tau*dire[1];
								p_smid   = mid[3] + 0.5*tau*dire[3];
								P_s_MIDy[i][j] = p_smid;
								RHO_F_sy[i][j] = z_smid*rho_smid*v_smid;
								U_F_sy[i][j]   = RHO_F_sy[i][j]*u_smid;
								V_F_sy[i][j]   = RHO_F_sy[i][j]*v_smid + z_smid*p_smid;
								E_F_sy[i][j]   = gamma_s/(gamma_s-1.0)*p_smid/rho_smid + 0.5*(u_smid*u_smid + v_smid*v_smid);
								E_F_sy[i][j]   = RHO_F_sy[i][j]*E_F_sy[i][j];
							}
					for(j = 0; j <  n_x; ++j)
						for(i = 0; i <= n_y; ++i)
							{
								if (i==0)
									{
										iL=0;
										iR=0;
									}
								else if (i==n_y)
									{
										iL=n_y-1;
										iR=n_y-1;
									}
								else
									{
										iL=i-1;
										iR=i;
									}
								z_sxL = Z_sx[iL][j];
								z_sxR = Z_sx[iR][j];
								z_syL = Z_sy[iL][j];
								z_syR = Z_sy[iR][j];
								z_sL =Z_sC[iL][j]+dy/2*z_syL;
								z_sR =Z_sC[iR][j]-dy/2*z_syR;
								if (j == n_x-1)
									jp1 = n_x-1;
								else
									jp1 = j+1;
								stag_RHO_F_sy[i][j] = 0.5*(U_R[i][j].rho_s*U_R[i][j].u_s+U_R[i][jp1].rho_s*U_R[i][jp1].u_s);
								if ((U_R[i][j].u_s + U_R[i][jp1].u_s) > 0.0)
									stag_ZRHO_F_sy[i][j] = stag_RHO_F_sy[i][j]*(z_sL - 0.25*tau*(-(U_R[i][j].v_s + U_R[i][jp1].v_s)*z_sxL+(U_R[i][j].u_s + U_R[i][jp1].u_s)*z_syL));
								else
									stag_ZRHO_F_sy[i][j] = stag_RHO_F_sy[i][j]*(z_sR - 0.25*tau*(-(U_R[i][j].v_s + U_R[i][jp1].v_s)*z_sxR+(U_R[i][j].u_s + U_R[i][jp1].u_s)*z_syR));
							}
					for(j = 0; j < n_x; ++j)
						for(i = 0; i < n_y; ++i)
							{
								if (i == n_y-1)
									ip1 = n_y-1;
								else
									ip1 = i+1;
								if (j == n_x-1)
									jp1 = n_x-1;
								else
									jp1 = j+1;
								RHO_s_cell  = 0.25*(RHO_sC[i][j]+RHO_sC[ip1][j]+RHO_sC[i][jp1]+RHO_sC[ip1][jp1]);
								//ZRHO_s_cell = 0.25*(ZRHO_sC[i][j]+ZRHO_sC[ip1][j]+ZRHO_sC[i][jp1]+ZRHO_sC[ip1][jp1]);
								ZRHO_s_cell = RHO_s_cell*Z_sC[i][j];
								RHO_s_cell  -=tau*(stag_RHO_F_sy[i+1][j] -stag_RHO_F_sy[i][j])/dy;
								ZRHO_s_cell -=tau*(stag_ZRHO_F_sy[i+1][j]-stag_ZRHO_F_sy[i][j])/dy;
								Z_sC[i][j]   =ZRHO_s_cell/RHO_s_cell;
							}
					for(j = 0; j < n_x; ++j)
						for(i = 0; i < n_y; ++i)
							{
								if (i == 0)
									i_1 = 0;
								else
									i_1 = i-1;
								if (j == 0)
									j_1 = 0;
								else
									j_1 = j-1;
								Z_sL[i][j] = 0.5*(Z_sC[i_1][j_1]+Z_sC[i_1][j]);
								Z_sR[i][j] = 0.5*(Z_sC[i][j_1]+Z_sC[i][j]);
							}
					for(j = 0; j < n_x; ++j)
						for(i = 0; i < n_y; ++i)
							{
								if(fabs(U_R[i][j].z_s-U_L[i][j].z_s)<eps)
									S = P_gC[i][j]*(U_R[i][j].z_s-U_L[i][j].z_s);
								else
									S = U_R[i][j].z_s*P_s_MIDy[i+1][j]-U_L[i][j].z_s*P_s_MIDy[i][j];
								ZRHO_gC[i][j] -= tau*(RHO_F_gy[i+1][j]-RHO_F_gy[i][j])/dy;
								RHO_U_gC[i][j]-= tau*(U_F_gy[i+1][j]  -U_F_gy[i][j])  /dy;
								RHO_V_gC[i][j]-= tau*(V_F_gy[i+1][j]  -V_F_gy[i][j])  /dy-tau/dx*S;
								E_gC[i][j]    -= tau*(E_F_gy[i+1][j]  -E_F_gy[i][j])  /dy-tau/dx*S*V_sC[i][j];
								ZRHO_sC[i][j] -= tau*(RHO_F_sy[i+1][j]-RHO_F_sy[i][j])/dy;
								RHO_U_sC[i][j]-= tau*(U_F_sy[i+1][j]  -U_F_sy[i][j])  /dy;
								RHO_V_sC[i][j]-= tau*(V_F_sy[i+1][j]  -V_F_sy[i][j])  /dy+tau/dx*S;
								E_sC[i][j]    -= tau*(E_F_sy[i+1][j]  -E_F_sy[i][j])  /dy+tau/dx*S*V_sC[i][j];
								area_L=0.5+V_sC[i][j]*tau/dx;
								area_R=1.0-area_L;
								U[0] = ZRHO_gC[i][j];
								U[1] = RHO_V_gC[i][j];
								U[2] = RHO_U_gC[i][j];
								U[3] = E_gC[i][j]-0.5*RHO_U_gC[i][j]*RHO_U_gC[i][j]/ZRHO_gC[i][j];
								U[4] = ZRHO_sC[i][j];
								U[5] = RHO_V_sC[i][j];
								U[6] = RHO_U_sC[i][j];
								U[7] = E_sC[i][j]-0.5*RHO_U_sC[i][j]*RHO_U_sC[i][j]/ZRHO_sC[i][j];
								primitive_comp(U, &U_L[i][j], &U_R[i][j], U_L[i][j].z_s, U_R[i][j].z_s, Z_sL[i][j], Z_sR[i][j], area_L, area_R);
								ZRHO_gC[i][j]  = 0.5*(U_L[i][j].U_rho_g+U_R[i][j].U_rho_g);
								RHO_U_gC[i][j] = 0.5*(U_L[i][j].U_v_g +U_R[i][j].U_v_g);
								RHO_V_gC[i][j] = 0.5*(U_L[i][j].U_u_g +U_R[i][j].U_u_g);
								E_gC[i][j]     = 0.5*(U_L[i][j].U_e_g +U_R[i][j].U_e_g);
								ZRHO_sC[i][j]  = 0.5*(U_L[i][j].U_rho_s+U_R[i][j].U_rho_s);
								RHO_U_sC[i][j] = 0.5*(U_L[i][j].U_v_s +U_R[i][j].U_v_s);
								RHO_V_sC[i][j] = 0.5*(U_L[i][j].U_u_s +U_R[i][j].U_u_s);
								E_sC[i][j]     = 0.5*(U_L[i][j].U_e_s +U_R[i][j].U_e_s);
							}
				}
			for(j = 0; j < n_x; ++j)
				{
					Z_sC[n_y-1][j]     = Z_sC[n_y-2][j];
					ZRHO_gC[n_y-1][j]  = ZRHO_gC[n_y-2][j];
					RHO_U_gC[n_y-1][j] = RHO_U_gC[n_y-2][j];
					RHO_V_gC[n_y-1][j] = RHO_V_gC[n_y-2][j];
					E_gC[n_y-1][j]     = E_gC[n_y-2][j];
					ZRHO_sC[n_y-1][j]  = ZRHO_sC[n_y-2][j];
					RHO_U_sC[n_y-1][j] = RHO_U_sC[n_y-2][j];
					RHO_V_sC[n_y-1][j] = RHO_V_sC[n_y-2][j];
					E_sC[n_y-1][j]     = E_sC[n_y-2][j];
				}
			for(i = 0; i < n_y; ++i)
				{
					Z_sC[i][n_x-1]     = Z_sC[i][n_x-2];
					ZRHO_gC[i][n_x-1]  = ZRHO_gC[i][n_x-2];
					RHO_U_gC[i][n_x-1] = RHO_U_gC[i][n_x-2];
					RHO_V_gC[i][n_x-1] = RHO_V_gC[i][n_x-2];
					E_gC[i][n_x-1]     = E_gC[i][n_x-2];
					ZRHO_sC[i][n_x-1]  = ZRHO_sC[i][n_x-2];
					RHO_U_sC[i][n_x-1] = RHO_U_sC[i][n_x-2];
					RHO_V_sC[i][n_x-1] = RHO_V_sC[i][n_x-2];
					E_sC[i][n_x-1]     = E_sC[i][n_x-2];
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
				FV->RHO[i*n_x+j]=RHO_gC[i][j];
				FV->U[i*n_x+j]  =U_gC[i][j];
				FV->V[i*n_x+j]  =V_gC[i][j];
				FV->P[i*n_x+j]  =P_gC[i][j];
				FV->Z_a[i*n_x+j]=Z_sC[i][j];
			}

	printf("\nThe cost of CPU time for the Eulerian method is %g seconds.\n", cpu_time);
}
