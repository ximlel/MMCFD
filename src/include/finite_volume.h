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


void FV_2_C_init(struct center_var C, struct flu_var FV);
void NewtonRapshon(double * x_star, double * err, double fun, double dfun, double eps);
void NewtonRapshon_matrix(double * x_star, double * err, double * fun, double * dfun, double eps);
void RI2U_cal(struct U_var * U, const struct RI_var * RI, double z_s, const double rho_g_start);
void U2RI_cal(const struct U_var * U, struct RI_var * RI);
void primitive_comp(double * U, struct U_var * U_L, struct U_var * U_R, double z_sL, double z_sR, double z_sL_out, double z_sR_out, double area_L, double area_R);

void BN_C2U(struct center_var C, double *U, int i, int j, int x_or_y);
void BN_ULR2prim(struct U_var U_L, struct U_var U_R, struct center_var C, int i, int j, int x_or_y);
void BN_ULR2cons(struct U_var U_L, struct U_var U_R, struct center_var C, int i, int j, int x_or_y);
void BN_RI2Cx(struct RI_var RI, struct center_var C, int i, int j);
void BN_RI2Cy(struct RI_var RI, struct center_var C, int i, int j);
void GRP_var_init(struct GRP_LR_var *G, struct slope_var SV, struct U_var U, double d, int i, int j, int pm_xy);
void GRP_RI_var_init(struct GRP_LR_var *G, struct slope_var SV, struct center_var C, double d, int i, int j, int pm_xy);
void G_LR_RI2U(struct GRP_LR_var *G, double z_s, int x_or_y);
void boundary_cond_x(struct center_var C, int cond);
void boundary_cond_y(struct center_var C, int cond);
