#ifndef VARSTRUC_H
#define VARSTRUC_H


#define EPS 0.000000000000001


#define N_CONF 400

extern double config[];

#define CONF_INI(i,j) printf("%3d-th configuration = %g .\n", i, j)

#define CONF_ERR(n)														\
	do {																\
		fprintf(stderr, "Error in the %d-th value of the configuration!\n", n); \
		exit(2);														\
	} while (0)


//fluid
struct flu_var {
	double *RHO, *U, *V, *W, *P, *PHI, *Z_a, *gamma, *RHO_b, *U_b, *V_b, *P_b;
};

//cell
struct cell_var {
	int **cell_cell;
	double **n_x, **n_y, **n_z;
	double **F_rho, **F_e, **F_gamma, **F_phi, **F_u, **F_v, **F_w;
	double  *U_rho,  *U_e,  *U_gamma,  *U_phi,  *U_u,  *U_v,  *U_w;
	double  **F_p_x,  **F_p_y, **RHO_p, **U_p, **V_p, **P_p, **PHI_p, **Z_a_p, **gamma_p;
	double **dt_U_p, **dt_V_p, **dt_F_p_x, **dt_F_p_y;
	double *X_c, *Y_c, *Z_c;
	double *vol, *c, *dist_p;
	double *gradx_rho,   *grady_rho,   *gradz_rho;
	double *gradx_phi,   *grady_phi,   *gradz_phi;
	double *gradx_gamma, *grady_gamma, *gradz_gamma;
	double *gradx_e,     *grady_e,     *gradz_e;
	double *gradx_u,     *grady_u,     *gradz_u;
	double *gradx_v,     *grady_v,     *gradz_v;
	double *gradx_w,     *grady_w,     *gradz_w;
	double **F_e_a, *U_e_a, *gradx_z_a, *grady_z_a, *gradz_z_a;
	double **RHO_star,    **P_star,    **U_qt_star,    **V_qt_star,    **gamma_star;
	double **RHO_minus_c, **P_minus_c, **U_qt_minus_c, **V_qt_minus_c, **gamma_minus_c;
	double **RHO_add_c,   **P_add_c,   **U_qt_add_c,   **V_qt_add_c,   **gamma_add_c;
	double **u_star, **u_minus_c, **u_add_c;
};

//interface
struct i_f_var {
	double n_x, n_y, n_z;
	double delta_x, delta_y, delta_z;
	double F_rho, F_e, F_gamma, F_phi, F_u, F_v, F_w;
	double U_rho, U_e, U_gamma, U_phi, U_u, U_v, U_w;
	double   RHO,   P,   gamma,   PHI,   U,   V,   W;
	double d_rho, d_phi, d_gamma, d_e, d_p, d_u, d_v, d_w;
	double t_rho, t_phi, t_gamma, t_e, t_p, t_u, t_v, t_w;
	double F_e_a, U_e_a, Z_a, d_z_a, t_z_a;
	double length;
	double RHO_star,    P_star,    U_qt_star,    V_qt_star,    gamma_star;
	double RHO_minus_c, P_minus_c, U_qt_minus_c, V_qt_minus_c, gamma_minus_c;
	double RHO_add_c,   P_add_c,   U_qt_add_c,   V_qt_add_c,   gamma_add_c;
	double u_star, u_minus_c, u_add_c;
	double RHO_int,  P_int,  U_int,  V_int;
};

//mesh
struct mesh_var {
	int num_pt, num_ghost, *cell_type, **cell_pt;
	int num_border[10], *border_pt, *border_cond, *peri_cell, *normal_v;
	double *X, *Y, *Z;
	void (*bc)(struct cell_var * cv, struct mesh_var mv, struct flu_var * FV, double t);
};

//===============================BN model===============================

//primitive and conservative variables
struct U_var {
	double     z_s;
	double   rho_s,   u_s,   v_s,   p_s,   rho_g,   u_g,   v_g,   p_g;
	double U_rho_s, U_u_s, U_v_s, U_e_s, U_rho_g, U_u_g, U_v_g, U_e_g;
};
//Riemann invariants
struct RI_var {
	double z_s, rho_s, u_s, Q, P, H, eta_g;
};
//variables at cell centers (including staggered solid cells)
struct center_var {
	//solid cell
	double **Z_sC;
	//gaseous cell
	double **RHO_gC, **P_gC, **U_gC, **V_gC;
	double **RHO_sC, **P_sC, **U_sC, **V_sC;
	double ** ZRHO_gC, ** ZRHO_sC;
	double **RHO_U_gC, **RHO_V_gC, **E_gC;
	double **RHO_U_sC, **RHO_V_sC, **E_sC;
	//gaseous cell in x/y direction
	double **Z_sS_xd, **Z_sS_yd; //S-staggered, d-direction
	double **Q_xd, **P_xd, **H_xd, **eta_g_xd;
	double **Q_yd, **P_yd, **H_yd, **eta_g_yd;
};
//slopes at centers of cells (including staggered solid cells)
struct slope_var {
	double **Z_sx, **Z_sy;
	double **RHO_gx, **P_gx, **U_gx, **V_gx;
	double **RHO_gy, **P_gy, **U_gy, **V_gy;
	double **RHO_sx, **P_sx, **U_sx, **V_sx;
	double **RHO_sy, **P_sy, **U_sy, **V_sy;
	double **Z_sS_x, **Z_sS_y;
	double **Q_x, **P_x, **H_x, **eta_g_x;
	double **Q_y, **P_y, **H_y, **eta_g_y;
};
//interfacial variables at t_{n+1} in GRP used to calculate slopes
struct face_var {
	double **Z_I_sxL, Z_I_sxR, Z_I_syL, Z_I_syR;
	double **RHO_I_gxL, **P_I_gxL, **U_I_gxL, **V_I_gxL;
	double **RHO_I_gxR, **P_I_gxR, **U_I_gxR, **V_I_gxR;
	double **RHO_I_gyL, **P_I_gyL, **U_I_gyL, **V_I_gyL;
	double **RHO_I_gyR, **P_I_gyR, **U_I_gyR, **V_I_gyR;
	double **RHO_I_sxL, **P_I_sxL, **U_I_sxL, **V_I_sxL;
	double **RHO_I_sxR, **P_I_sxR, **U_I_sxR, **V_I_sxR;
	double **RHO_I_syL, **P_I_syL, **U_I_syL, **V_I_syL;
	double **RHO_I_syR, **P_I_syR, **U_I_syR, **V_I_syR;
};
//left and right variables for GRO solver (both Edir and BN-Riemann-invariant solver)
struct GRP_LR_var {
	double rho_g,  p_g,  u_g,  v_g;
	double rho_gx, p_gx, u_gx, v_gx;
	double rho_gy, p_gy, u_gy, v_gy;
	double rho_s,  p_s,  u_s,  v_s;
	double rho_sx, p_sx, u_sx, v_sx;
	double rho_sy, p_sy, u_sy, v_sy;
	double Q,  P,  H,  eta_g;
	double Qx, Px, Hx, eta_gx;
	double Qy, Py, Hy, eta_gy;
};
//flux and midpoint (t_{n+1/2}) variables at cell interfaces (including solid cell interfaces)
struct flux_var {
	//flux
	double **ZRHO_F_gx, **E_F_gx, **U_F_gx, **V_F_gx;
	double **ZRHO_F_gy, **E_F_gy, **U_F_gy, **V_F_gy;
	double **ZRHO_F_sx, **E_F_sx, **U_F_sx, **V_F_sx;
	double **ZRHO_F_sy, **E_F_sy, **U_F_sy, **V_F_sy;
	//flux (solid-cell)
	double **stag_RHO_F_sx, **stag_ZRHO_F_sx;
	double **stag_RHO_F_sy, **stag_ZRHO_F_sy;
	//midpoint
  	double **Z_s_MIDx, **Z_s_MIDy;
  	double **P_g_MIDx, **U_g_MIDx, **V_g_MIDx;
  	double **P_s_MIDx, **U_s_MIDx, **V_s_MIDx;
	double **P_g_MIDy, **U_g_MIDy, **V_g_MIDy;
	double **P_s_MIDy, **U_s_MIDy, **V_s_MIDy;
};

#endif
