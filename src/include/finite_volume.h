void cons_qty_init(const struct cell_var * cv, const struct flu_var * FV);
int cons2prim(struct i_f_var * ifv);
int cons_qty_update(const struct cell_var * cv, const struct mesh_var * mv,
					const struct flu_var *  FV, const double tau);
int cons_qty_update_corr_ave_P(struct cell_var * cv, const struct mesh_var * mv,
							   const struct flu_var * FV, double tau, const int stop_step);


struct cell_var cell_mem_init(const struct mesh_var * mv, struct flu_var * FV);
void vol_comp(const struct cell_var * cv, const struct mesh_var * mv);
void cell_pt_clockwise(const struct mesh_var * mv);
void cell_rel(const struct cell_var * cv, const struct mesh_var * mv);
void cell_centroid(const struct cell_var * cv, const struct mesh_var * mv);


void slope_limiter(const struct cell_var * cv,const struct mesh_var * mv, const struct flu_var * FV);
void slope_limiter2(const struct cell_var * cv,const struct mesh_var * mv, const struct flu_var * FV);


void cons_qty_copy_cv2ifv(struct i_f_var * ifv, const struct cell_var * cv, const int c);
void cons_qty_copy_ifv2cv(const struct i_f_var * ifv, struct cell_var * cv, const int c);
void prim_var_copy_ifv2FV(const struct i_f_var * ifv, const struct flu_var * FV,const int c);
void flux_copy_ifv2cv(const struct i_f_var * ifv, const struct cell_var *cv, const int k, const int j);
void flux_add_ifv2cv(const struct i_f_var * ifv, const struct cell_var * cv, const int k, const int j);


int fluid_var_update(struct flu_var *FV, struct cell_var *cv);
int interface_var_init(const struct cell_var * cv, const struct mesh_var * mv,
					   struct i_f_var * ifv, struct i_f_var * ifv_R,
					   const int k, const int j, const int i, const double gauss);
double tau_calc(const struct cell_var * cv, const struct mesh_var * mv);


void Roe_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R);
void HLL_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R);
void Riemann_exact_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R);
void GRP_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R, const double tau);
void GRP_2D_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R, const double tau);


void finite_volume_scheme(struct flu_var * FV, const struct mesh_var * mv, const char * scheme, const char * problem);
void finite_volume_scheme_GRP2D(struct flu_var * FV, const struct mesh_var * mv, const char * scheme, const char * problem);


void init_data_1(double * Z_l, double * U_RHO_g, double * U_U_g, double * U_V_g, double * U_E_g,
                 double * U_RHO_l, double * U_U_l, double * U_V_l, double * U_E_l);
void init_data_2(double * Z_l, double * U_RHO_g, double * U_U_g, double * U_V_g, double * U_E_g,
                 double * U_RHO_l, double * U_U_l, double * U_V_l, double * U_E_l);
void init_data_3(double * Z_l, double * U_RHO_g, double * U_U_g, double * U_V_g, double * U_E_g,
                 double * U_RHO_l, double * U_U_l, double * U_V_l, double * U_E_l);
void U_init(int n_x, double ZRHO_gC[][n_x], double RHO_U_gC[][n_x], double RHO_V_gC[][n_x], double E_gC[][n_x],
            double ZRHO_lC[][n_x], double RHO_U_lC[][n_x], double RHO_V_lC[][n_x], double E_lC[][n_x], int i, int j, 
            double * U_RHO_g, double * U_U_g, double * U_V_g, double * U_E_g, 
            double * U_RHO_l, double * U_U_l, double * U_V_l, double * U_E_l, int a, int b, int c, int d);
void NewtonRapshon(double * x_star, double * err, double fun, double dfun, double x0, double eps);
void NewtonRapshon_matrix(double * x_star, double * err, double * fun, double * dfun, double * x0, double eps);
void RI2U_cal(struct U_var * U, const struct RI_var * RI, double phi_s, const double lo_g_start);
void U2RI_cal(const struct U_var * U, struct RI_var * RI);
void primitive_comp(double * U, struct U_var * U_L, struct U_var * U_R, double phi_sL, double phi_sR, double phi_sL_out, double phi_sR_out, double area_L, double area_R);
