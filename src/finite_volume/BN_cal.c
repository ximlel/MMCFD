#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "../include/var_struc.h"
#include "../include/tools.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>


//Newton-Rapshon iteration
void NewtonRapshon(double * x_star, double * err, double fun, double dfun, double eps)
{
    double d;
    if (fabs(fun) <= eps)
	d = 0.0;
    else
	d = -fun/dfun;
    * x_star = * x_star + d;
    * err = fabs(d);
}

static inline double V_norm(double * x)
{
    return sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3]);
}

void NewtonRapshon_matrix(double * x_star, double * err, double * fun, double * dfun, double eps)
{
    double d[4]={0.0};
    rinv(dfun,4); //Matrix inv of dfun
    int i,j;
    if (V_norm(fun) > eps) {
	for(i=0; i<4; i++)
	    for(j=0; j<4; j++)
		d[i]-=fun[j]*(*(dfun+i+j*4));
    }
    //V_add(x_star, x_star, d);
    for(i=0; i<4; i++)
	x_star[i] = x_star[i]+d[i];
    * err = V_norm(d);
}

//From Riemann invariants to calculate var U.
void RI2U_cal(struct U_var * U, const struct RI_var * RI, double z_s, const double rho_g_start)
{
    const double eps = config[4];
    const double gama_g = config[106];
    double z_g = 1.0-z_s;
    double rho_s=RI->rho_s;
    double u_s=RI->u_s;
    double Q=RI->Q;
    double P=RI->P;
    double H=RI->H;
    double eta_g=RI->eta_g;
    int it_max = 500, k = 0;
    double err1 = 1e50;
    double rho_g=rho_g_start;
    double fun, dfun;
    while (k<it_max && err1>eps) {					
	fun  = H-0.5*pow(Q/z_g,2)/pow(rho_g,2)-gama_g/(gama_g-1.0)*eta_g*pow(rho_g,gama_g-1.0);
	dfun = pow(Q/z_g,2)/pow(rho_g,3)-gama_g*eta_g*pow(rho_g,gama_g-2.0);
	NewtonRapshon(&rho_g, &err1, fun,dfun,eps);
	k=k+1;
    }
    if(k>=it_max)
        printf("\nRI2P_err:%lf! %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n", err1, z_s,rho_s,u_s,Q,P,H,eta_g,rho_g);
    U->p_g = pow(rho_g,gama_g)*eta_g;
    U->u_g = Q/z_g/rho_g+u_s;
    U->p_s = (P-Q*(U->u_g-u_s)-z_g*U->p_g)/z_s;
    U->rho_g = rho_g;
    U->rho_s = rho_s;
    U->u_s = u_s;
}

//From var U to calculate Riemann invariants.
void U2RI_cal(const struct U_var * U, struct RI_var * RI)
{
    const double eps = config[4];
    const double gama_g = config[106];
    double z_s = U->z_s, rho_s = U->rho_s, u_s = U->u_s, p_s = U->p_s, rho_g = U->rho_g, u_g = U->u_g, p_g = U->p_g;
    double z_g = 1.0-z_s;
	
    RI->eta_g=p_g/pow(rho_g,gama_g);
    RI->Q=z_g*rho_g*(u_g-u_s);
    RI->P=z_g*rho_g*pow(u_g-u_s,2)+z_g*p_g+z_s*p_s;
    RI->H=0.5*pow(u_g-u_s,2)+gama_g/(gama_g-1.0)*p_g/rho_g;
    RI->rho_s=rho_s;
    RI->u_s=u_s;
}

//compute primitive var
void primitive_comp(double * U, struct U_var * U_L, struct U_var * U_R, double z_sL, double z_sR, double z_sL_out, double z_sR_out, double area_L, double area_R)
{   
    double z_gL=1-z_sL;
    double z_gR=1-z_sR;
    const double gama_g = config[106], gama_s = config[6];
    double eps = config[4];
    double z_s = area_L*z_sL+area_R*z_sR;
    double z_g = 1.0-z_s;
    double U1=U[0], U2=U[1], U3=U[3], U4=U[4], U5=U[5], U6=U[7];
    double rho_gR = U1/z_g;
    double u_gR  = U2/U1;
    double p_gR  = (U3/z_g - 0.5*rho_gR*pow(u_gR,2))*(gama_g-1.0);
    double rho_s  = U4/z_s;
    double u_s   = U5/U4;
    double p_sR  = (U6/z_s - 0.5*rho_s*pow(u_s,2))*(gama_s-1.0);
    U_L->v_g = U[2]/U1;
    U_R->v_g = U_L->v_g;	
    U_L->v_s = U[6]/U4;
    U_R->v_s = U_L->v_s;
    U_L->z_s = z_sL;
    U_R->z_s = z_sR;
    double fun[4], dfun[4][4], x_star[4];
    int it_max = 5000, k = 0;
    double err2 = 1e50;
    while (k<it_max && err2>eps && fabs(z_sL-z_sR)>eps) {			
	fun[0] = U3-area_L*z_gL*(0.5*((U1-area_R*z_gR*rho_gR)/area_L/z_gL)*pow((U2-area_R*z_gR*rho_gR*u_gR)/(U1-area_R*z_gR*rho_gR),2.0)+(p_gR/pow(rho_gR,gama_g)*pow((U1-area_R*z_gR*rho_gR)/area_L/z_gL,gama_g))/(gama_g-1.0))-area_R*z_gR*(0.5*rho_gR*pow(u_gR,2.0)+p_gR/(gama_g-1.0));
	fun[1] = z_gL*((U1-area_R*z_gR*rho_gR)/area_L/z_gL)*(((U2-area_R*z_gR*rho_gR*u_gR)/(U1-area_R*z_gR*rho_gR))-u_s)-z_gR*rho_gR*(u_gR-u_s);
	fun[2] = z_gL*((U1-area_R*z_gR*rho_gR)/area_L/z_gL)*pow(((U2-area_R*z_gR*rho_gR*u_gR)/(U1-area_R*z_gR*rho_gR))-u_s,2.0)+z_gL*(p_gR/pow(rho_gR,gama_g)*pow((U1-area_R*z_gR*rho_gR)/area_L/z_gL,gama_g))+z_sL*(((U6-0.5*z_s*rho_s*pow(u_s,2.0))*(gama_s-1.0)-area_R*z_sR*p_sR)/area_L/z_sL)-z_gR*rho_gR*pow(u_gR-u_s,2.0)-z_gR*p_gR-z_sR*p_sR;
	fun[3] = 0.5*pow(((U2-area_R*z_gR*rho_gR*u_gR)/(U1-area_R*z_gR*rho_gR))-u_s,2.0)+gama_g/(gama_g-1.0)*(p_gR/pow(rho_gR,gama_g)*pow((U1-area_R*z_gR*rho_gR)/area_L/z_gL,gama_g))/((U1-area_R*z_gR*rho_gR)/area_L/z_gL)-0.5*pow(u_gR-u_s,2.0)-gama_g/(gama_g-1.0)*p_gR/rho_gR;
	dfun[0][0] = area_L*z_gL*((gama_g*p_gR*pow((U1 - area_R*rho_gR*z_gR)/(area_L*z_gL),gama_g))/(pow(rho_gR,gama_g + 1.0)*(gama_g - 1.0)) - (area_R*z_gR*pow(U2 - area_R*rho_gR*z_gR*u_gR,2.0))/(2.0*area_L*z_gL*pow(U1 - area_R*rho_gR*z_gR,2.0)) + (area_R*z_gR*u_gR*(U2 - area_R*rho_gR*z_gR*u_gR))/(area_L*z_gL*(U1 - area_R*rho_gR*z_gR)) + (area_R*gama_g*p_gR*z_gR*pow((U1 - area_R*rho_gR*z_gR)/(area_L*z_gL),gama_g - 1.0))/(area_L*pow(rho_gR,gama_g)*z_gL*(gama_g - 1.0))) - (area_R*z_gR*pow(u_gR,2.0))/2.0;
	dfun[1][0] = - (area_R*z_gR)/(gama_g - 1.0) - (area_L*z_gL*pow((U1 - area_R*rho_gR*z_gR)/(area_L*z_gL),gama_g))/(pow(rho_gR,gama_g)*(gama_g - 1));
	dfun[2][0] = (area_R*rho_gR*z_gR*(U2 - area_R*rho_gR*z_gR*u_gR))/(U1 - area_R*rho_gR*z_gR) - area_R*rho_gR*z_gR*u_gR;
	dfun[3][0] = 0.0;
	dfun[0][1] = (area_R*z_gR*(u_s - (U2 - area_R*rho_gR*z_gR*u_gR)/(U1 - area_R*rho_gR*z_gR)))/area_L - ((U1 - area_R*rho_gR*z_gR)*((area_R*z_gR*u_gR)/(U1 - area_R*rho_gR*z_gR) - (area_R*z_gR*(U2 - area_R*rho_gR*z_gR*u_gR))/pow(U1 - area_R*rho_gR*z_gR,2.0)))/area_L - z_gR*(u_gR - u_s);
	dfun[1][1] = 0.0;
	dfun[2][1] = - rho_gR*z_gR - (area_R*rho_gR*z_gR)/area_L;
	dfun[3][1] = 0.0;
	dfun[0][2] = (2*(U1 - area_R*rho_gR*z_gR)*((area_R*z_gR*u_gR)/(U1 - area_R*rho_gR*z_gR) - (area_R*z_gR*(U2 - area_R*rho_gR*z_gR*u_gR))/pow(U1 - area_R*rho_gR*z_gR,2.0))*(u_s - (U2 - area_R*rho_gR*z_gR*u_gR)/(U1 - area_R*rho_gR*z_gR)))/area_L - z_gR*pow(u_gR - u_s,2.0) - (area_R*z_gR*pow(u_s - (U2 - area_R*rho_gR*z_gR*u_gR)/(U1 - area_R*rho_gR*z_gR),2))/area_L - (gama_g*p_gR*z_gL*pow((U1 - area_R*rho_gR*z_gR)/(area_L*z_gL),gama_g))/pow(rho_gR,gama_g + 1.0) - (area_R*gama_g*p_gR*z_gR*pow((U1 - area_R*rho_gR*z_gR)/(area_L*z_gL),gama_g - 1.0))/(area_L*pow(rho_gR,gama_g));
	dfun[1][2] = (z_gL*pow((U1 - area_R*rho_gR*z_gR)/(area_L*z_gL),gama_g))/pow(rho_gR,gama_g) - z_gR;
	dfun[2][2] = (2.0*area_R*rho_gR*z_gR*(u_s - (U2 - area_R*rho_gR*z_gR*u_gR)/(U1 - area_R*rho_gR*z_gR)))/area_L - rho_gR*z_gR*(2.0*u_gR - 2.0*u_s);
	dfun[3][2] = - z_sR - (area_R*z_sR)/area_L;
	dfun[0][3] = ((area_R*z_gR*u_gR)/(U1 - area_R*rho_gR*z_gR) - (area_R*z_gR*(U2 - area_R*rho_gR*z_gR*u_gR))/pow(U1 - area_R*rho_gR*z_gR,2.0))*(u_s - (U2 - area_R*rho_gR*z_gR*u_gR)/(U1 - area_R*rho_gR*z_gR)) + (gama_g*p_gR)/(pow(rho_gR,2.0)*(gama_g - 1.0)) - (area_L*pow(gama_g,2.0)*p_gR*z_gL*pow((U1 - area_R*rho_gR*z_gR)/(area_L*z_gL),gama_g))/(pow(rho_gR,gama_g + 1.0)*(U1 - area_R*rho_gR*z_gR)*(gama_g - 1.0)) - (area_R*pow(gama_g,2.0)*p_gR*z_gR*pow((U1 - area_R*rho_gR*z_gR)/(area_L*z_gL),gama_g - 1.0))/(pow(rho_gR,gama_g)*(U1 - area_R*rho_gR*z_gR)*(gama_g - 1.0)) + (area_L*area_R*gama_g*p_gR*z_gL*z_gR*pow((U1 - area_R*rho_gR*z_gR)/(area_L*z_gL),gama_g))/(pow(rho_gR,gama_g)*pow(U1 - area_R*rho_gR*z_gR,2.0)*(gama_g - 1.0));
	dfun[1][3] = (area_L*gama_g*z_gL*pow((U1 - area_R*rho_gR*z_gR)/(area_L*z_gL),gama_g))/(pow(rho_gR,gama_g)*(U1 - area_R*rho_gR*z_gR)*(gama_g - 1)) - gama_g/(rho_gR*(gama_g - 1.0));
	dfun[2][3] = u_s - u_gR + (area_R*rho_gR*z_gR*(u_s - (U2 - area_R*rho_gR*z_gR*u_gR)/(U1 - area_R*rho_gR*z_gR)))/(U1 - area_R*rho_gR*z_gR);
	dfun[3][3] = 0.0;
	x_star[0] = rho_gR;
	x_star[1] = p_gR;
	x_star[2] = u_gR;
	x_star[3] = p_sR;
	NewtonRapshon_matrix(x_star, &err2, fun, dfun[0], eps);
	rho_gR=fmax(x_star[0],eps);
	rho_gR=fmin(rho_gR,U1/area_R/z_gR-eps);
	p_gR =fmax(x_star[1],eps);
	u_gR =x_star[2];
	p_sR =fmax(x_star[3],eps);
	p_sR =fmin(p_sR,(U6-0.5*z_s*rho_s*pow(u_s,2))*(gama_s-1.0)/area_R/z_sR-eps);
	k=k+1;
    }
    if (k>=it_max)
        printf("\nRIeq_err:%lf! %lf, %lf, %lf, %lf, %lf, %lf\n",err2,z_sL,z_sR,rho_gR,p_gR,u_gR,p_sR);

    U_L->rho_g = (U1-area_R*z_gR*rho_gR)/area_L/z_gL;
    U_R->rho_g = rho_gR;
    U_L->u_g = (U2-area_R*z_gR*rho_gR*u_gR)/(U1-area_R*z_gR*rho_gR);
    U_R->u_g = u_gR;
    U_L->p_s = ((U6-0.5*z_s*rho_s*pow(u_s,2.0))*(gama_s-1)-area_R*z_sR*p_sR)/area_L/z_sL;
    U_L->p_g = p_gR/pow(rho_gR,gama_g)*pow((U1-area_R*z_gR*rho_gR)/area_L/z_gL,gama_g);
    U_R->p_s = p_sR;
    U_R->p_g = p_gR;	
    U_L->rho_s= rho_s;
    U_R->rho_s= rho_s;
    U_L->u_s = u_s;
    U_R->u_s = u_s;
    
    struct RI_var RI_L, RI_R;
    U2RI_cal(U_L, &RI_L);
    RI2U_cal(U_L, &RI_L, z_sL_out, U_L->rho_g);
    U2RI_cal(U_R, &RI_R);
    RI2U_cal(U_R, &RI_R, z_sR_out, U_R->rho_g);
    
    U_L->U_rho_g = (1.0-z_sL_out)*U_L->rho_g;
    U_L->U_u_g  = U_L->U_rho_g*U_L->u_g;
    U_L->U_v_g  = U_L->U_rho_g*U_L->v_g;
    U_L->U_e_g  = U_L->U_rho_g*(U_L->p_g/U_L->rho_g/(gama_g-1.0)+0.5*U_L->u_g*U_L->u_g+0.5*U_L->v_g*U_L->v_g);
    U_L->U_rho_s = z_sL_out*U_L->rho_s;
    U_L->U_u_s  = U_L->U_rho_s*U_L->u_s;
    U_L->U_v_s  = U_L->U_rho_s*U_L->v_s;
    U_L->U_e_s  = U_L->U_rho_s*(U_L->p_s/U_L->rho_s/(gama_s-1.0)+0.5*U_L->u_s*U_L->u_s+0.5*U_L->v_s*U_L->v_s);
    U_R->U_rho_g = (1.0-z_sR_out)*U_R->rho_g;
    U_R->U_u_g  = U_R->U_rho_g*U_R->u_g;
    U_R->U_v_g  = U_R->U_rho_g*U_R->v_g;
    U_R->U_e_g  = U_R->U_rho_g*(U_R->p_g/U_R->rho_g/(gama_g-1.0)+0.5*U_R->u_g*U_R->u_g+0.5*U_R->v_g*U_R->v_g);
    U_R->U_rho_s = z_sR_out*U_R->rho_s;
    U_R->U_u_s  = U_R->U_rho_s*U_R->u_s;
    U_R->U_v_s  = U_R->U_rho_s*U_R->v_s;
    U_R->U_e_s  = U_R->U_rho_s*(U_R->p_s/U_R->rho_s/(gama_s-1.0)+0.5*U_R->u_s*U_R->u_s+0.5*U_R->v_s*U_R->v_s);	
}

void Lagrangian_with_Multiplier()
{
    const int N_x = 4, N_lmd = 3;
    int i,j,l;
    gsl_matrix *D2_xx_L_c_k = gsl_matrix_alloc(N_x, N_x), *H_k = gsl_matrix_alloc(N_x, N_x);
    gsl_matrix *N_k = gsl_matrix_alloc(N_x, N_lmd);
    // H_k = D2_xx_L, N_k = D_h
    gsl_vector *D_x_L_c_k = gsl_vector_alloc(N_x);
    gsl_vector *D_L = gsl_vector_alloc(N_x+N_lmd);
    gsl_vector *lambda_k = gsl_vector_alloc(N_lmd), *lambda_k_b = gsl_vector_alloc(N_lmd);
    gsl_vector *x_k = gsl_vector_alloc(N_x), *x_k_b = gsl_vector_alloc(N_x);
    gsl_vector *d_k = gsl_vector_alloc(N_x), *h_k = gsl_vector_alloc(N_x), *c_k_h = gsl_vector_alloc(N_x);
    gsl_permutation *per = gsl_permutation_alloc(N_x);
    double c_k, vareps_k, omega_k, m_k;
    int m_idx;
    double L_c_k, L_c_k_beta, ddot;
    const double gamma = 0.5, r = 2, beta = 0.5, sigma = 0.25;
    
    for (i = 0; i < N_x; i++)
	for (j = 0; j < N_x; j++)
	    gsl_matrix_set(D2_xx_L_c_k, i, j, 10086);
    for (i = 0; i < N_lmd; i++)
	gsl_vector_set(D_x_L_c_k, i, 10086);    
    gsl_vector_scale(D_x_L_c_k,-1.0);
    /* Modified Cholesky */
    gsl_linalg_mcholesky_decomp(D2_xx_L_c_k, per, NULL);
    gsl_linalg_mcholesky_solve(D2_xx_L_c_k, per, D_x_L_c_k, d_k);

    gsl_vector_memcpy(x_k_b, x_k);
    gsl_vector_add(x_k_b, d_k);
    gsl_vector_memcpy(lambda_k_b, lambda_k);
    gsl_vector_memcpy(c_k_h, h_k);
    gsl_vector_scale(c_k_h, c_k);
    // gsl_vector_add(lambda_k_b, c_h_k);


    
    if (pow(gsl_blas_dnrm2(D_L),2) < omega_k) {
	gsl_vector_memcpy(x_k, x_k_b);
	gsl_vector_memcpy(lambda_k,lambda_k_b);
	omega_k = gamma*pow(gsl_blas_dnrm2(D_L),2);
    }
    else {
	m_k = 0;
	m_idx = 1;
	while (m_idx) {
	    gsl_blas_ddot(d_k, D_x_L_c_k, &ddot);
	    if ((L_c_k - L_c_k_beta) >= -(sigma*pow(beta,m_k)*ddot))
		m_idx = 0;
	    else
		m_k++;
	    if (m_k>50)
		printf("m_k is bigger than 50!\n");
	}
	if (gsl_blas_dnrm2(D_x_L_c_k) <= vareps_k)
	    {
		gsl_vector_add(lambda_k,c_k_h);
		vareps_k *= gamma;
		c_k *= r;
		omega_k = gamma*pow(gsl_blas_dnrm2(D_L),2);
	    }
    }
    
    gsl_matrix_free(D2_xx_L_c_k);
    gsl_matrix_free(H_k);
    gsl_matrix_free(N_k);
    gsl_vector_free(D_x_L_c_k);
    gsl_vector_free(D_L);
    gsl_vector_free(lambda_k); gsl_vector_free(lambda_k_b);
    gsl_vector_free(x_k); gsl_vector_free(x_k_b);
    gsl_vector_free(d_k); gsl_vector_free(h_k);
    gsl_permutation_free(per);
}
