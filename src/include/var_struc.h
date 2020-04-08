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
	double   RHO_int,   P_int,    U_int,   V_int;
};

//Riemann invariants
struct RI_var {
	double rho_s, u_s, Q, P, H, eta_g;
};
//Riemann invariants
struct U_var {
	double z_s, rho_s, u_s, v_s, p_s, z_g, rho_g, u_g, v_g, p_g;
	double U_rho_s, U_u_s, U_v_s, U_e_s, U_rho_g, U_u_g, U_v_g, U_e_g;
};

//mesh
struct mesh_var {
	int num_pt, num_ghost, *cell_type, **cell_pt;
	int num_border[10], *border_pt, *border_cond, *peri_cell, *normal_v;
	double *X, *Y, *Z;
	void (*bc)(struct cell_var * cv, struct mesh_var mv, struct flu_var * FV, double t);
};
	
#endif
