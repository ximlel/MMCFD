#ifndef Riemann_solver_H
#define Riemann_solver_H


void Roe_solver(double *F, double gamma, double P_L, double RHO_L, double U_L, double P_R, double RHO_R, double U_R, double *lambda_max, double delta);
double Riemann_solver_exact(double * U_star, double * P_star, double gammaL, double gammaR, double u_L, double u_R, double p_L, double p_R, double c_L, double c_R, int * CRW, double eps, double tol, int N);

void Roe_2D_solver(double *F, double gamma, double P_L, double RHO_L, double U_L, double V_L,double n_x, double n_y, double P_R, double RHO_R, double U_R,double V_R,double *lambda_max, double delta);
void HLL_solver(double *F, double gamma, double P_L, double RHO_L, double U_L, double V_L,double n_x, double n_y, double P_R, double RHO_R, double U_R,double V_R,double *lambda_max);
void Roe_Goundov_solver(double *F, double gamma, double P_L, double RHO_L, double U_L, double P_R, double RHO_R, double U_R, double *lambda_max, double delta);
void Roe_HLL_solver(double *V_mk, double *F, double gamma, double P_L, double RHO_L, double U_L, double V_L, double P_R, double RHO_R, double U_R, double V_R, double *lambda_max, double delta);

void linear_GRP_solver_Edir_Q1D
(double *wave_speed, double *D, double *U, double *U_star, const double lambda_u, const double lambda_v,
 const double rho_L, const double rho_R, const double d_rho_L, const double d_rho_R, const double t_rho_L, const double t_rho_R,
 const double   u_L, const double   u_R, const double   d_u_L, const double   d_u_R, const double   t_u_L, const double   t_u_R,
 const double   v_L, const double   v_R, const double   d_v_L, const double   d_v_R, const double   t_v_L, const double   t_v_R,
 const double   p_L, const double   p_R, const double   d_p_L, const double   d_p_R, const double   t_p_L, const double   t_p_R,
 const double   z_L, const double   z_R, const double   d_z_L, const double   d_z_R, const double   t_z_L, const double   t_z_R,
 const double phi_L, const double phi_R, const double d_phi_L, const double d_phi_R, const double t_phi_L, const double t_phi_R,
 const double gammaL, const double gammaR, const double  eps, const double  atc);
void linear_GRP_solver_Edir_G2D
(double *wave_speed, double *D, double *U, double *U_star, const double lambda_u, const double lambda_v,
 const double rho_L, const double rho_R, const double d_rho_L, const double d_rho_R, const double t_rho_L, const double t_rho_R,
 const double   u_L, const double   u_R, const double   d_u_L, const double   d_u_R, const double   t_u_L, const double   t_u_R,
 const double   v_L, const double   v_R, const double   d_v_L, const double   d_v_R, const double   t_v_L, const double   t_v_R,
 const double   p_L, const double   p_R, const double   d_p_L, const double   d_p_R, const double   t_p_L, const double   t_p_R,
 const double   z_L, const double   z_R, const double   d_z_L, const double   d_z_R, const double   t_z_L, const double   t_z_R,
 const double phi_L, const double phi_R, const double d_phi_L, const double d_phi_R, const double t_phi_L, const double t_phi_R,
 const double gammaL, const double gammaR, const double  eps, const double  atc);

int Roe_GRP_solver_BN
(double *Dt_U_all, double *U_all,
 const double   z_g_L, const double   z_g_R, const double   d_z_g_L, const double   d_z_g_R, const double   t_z_g_L, const double   t_z_g_R,
 const double rho_g_L, const double rho_g_R, const double d_rho_g_L, const double d_rho_g_R, const double t_rho_g_L, const double t_rho_g_R,
 const double   u_g_L, const double   u_g_R, const double   d_u_g_L, const double   d_u_g_R, const double   t_u_g_L, const double   t_u_g_R,
 const double   v_g_L, const double   v_g_R, const double   d_v_g_L, const double   d_v_g_R, const double   t_v_g_L, const double   t_v_g_R,
 const double   p_g_L, const double   p_g_R, const double   d_p_g_L, const double   d_p_g_R, const double   t_p_g_L, const double   t_p_g_R,
 const double rho_l_L, const double rho_l_R, const double d_rho_l_L, const double d_rho_l_R, const double t_rho_l_L, const double t_rho_l_R,
 const double   u_l_L, const double   u_l_R, const double   d_u_l_L, const double   d_u_l_R, const double   t_u_l_L, const double   t_u_l_R,
 const double   v_l_L, const double   v_l_R, const double   d_v_l_L, const double   d_v_l_R, const double   t_v_l_L, const double   t_v_l_R,
 const double   p_l_L, const double   p_l_R, const double   d_p_l_L, const double   d_p_l_R, const double   t_p_l_L, const double   t_p_l_R,
 const double gamma_g, const double gamma_l, const double  eps);

int Roe_GRP_solver_BN_1D
(double *Dt_U_all, double *U_all,
 const double   z_g_L, const double   z_g_R, const double   d_z_g_L, const double   d_z_g_R, const double   t_z_g_L, const double   t_z_g_R,
 const double rho_g_L, const double rho_g_R, const double d_rho_g_L, const double d_rho_g_R, const double t_rho_g_L, const double t_rho_g_R,
 const double   u_g_L, const double   u_g_R, const double   d_u_g_L, const double   d_u_g_R, const double   t_u_g_L, const double   t_u_g_R,
 const double   v_g_L, const double   v_g_R, const double   d_v_g_L, const double   d_v_g_R, const double   t_v_g_L, const double   t_v_g_R,
 const double   p_g_L, const double   p_g_R, const double   d_p_g_L, const double   d_p_g_R, const double   t_p_g_L, const double   t_p_g_R,
 const double rho_l_L, const double rho_l_R, const double d_rho_l_L, const double d_rho_l_R, const double t_rho_l_L, const double t_rho_l_R,
 const double   u_l_L, const double   u_l_R, const double   d_u_l_L, const double   d_u_l_R, const double   t_u_l_L, const double   t_u_l_R,
 const double   v_l_L, const double   v_l_R, const double   d_v_l_L, const double   d_v_l_R, const double   t_v_l_L, const double   t_v_l_R,
 const double   p_l_L, const double   p_l_R, const double   d_p_l_L, const double   d_p_l_R, const double   t_p_l_L, const double   t_p_l_R,
 const double gamma_g, const double gamma_l, const double  eps);

int linear_GRP_RI_solver_BN
(struct RI_var *RI, const double D_z_s, const double z_s, const double *mid_g, const double *mid_s, 
 const struct GRP_LR_var GL, const struct GRP_LR_var GR,
 const double gamma_s, const double gamma_g, const double eps, const double tau, const int x_or_y);

#endif
