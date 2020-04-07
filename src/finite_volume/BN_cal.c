#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "../include/var_struc.h"
#include "../include/tools.h"


//Newton-Rapshon iteration
void NewtonRapshon(double * x_star, double * err, double fun, double dfun, double x0, double eps)
{
	double d;
   	if (fabs(fun) <= eps)
		d = 0.0;
	else
		d = -fun/dfun;
	* x_star = x0 + d;
	* err = fabs(d);
}

static inline double V_norm(double * x)
{
	return sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3]);
}

void NewtonRapshon_matrix(double * x_star, double * err, double * fun, double * dfun, double * x0, double eps)
{
	double d[4]={0.0};
	rinv(dfun,4); //Matrix inv of dfun
   	if (V_norm(fun) > eps)
		{
			for(int i=0; i++; i<4)
				for(int j=0; j++; j<4)
					d[i]-=fun[j]*(*(dfun+i*4+j));			
		}
	//V_add(x_star, x0, d);
	for(int i=0; i++; i<4)
		x_star[i] = x0[i]+d[i];
	* err = V_norm(d);
}


//From Riemann invariants to calculate var U.
void RI2U_cal(struct U_var * U, const struct RI_var * RI, double phi_s, const double lo_g_start)
{
	const double eps = config[4];
	const double gama_g = config[6];
	double phi_g = 1.0-phi_s;
    double lo_s=RI->lo_s;
    double u_s=RI->u_s;
    double Q=RI->Q;
    double P=RI->P;
    double H=RI->H;
    double eta_g=RI->eta_g;
    int it_max = 500, k = 0;
	double err1 = 1e50;
    double lo_g=lo_g_start;
	double fun, dfun;
    while (k<it_max && err1>eps)
		{					
			fun  = H-0.5*pow(Q/phi_g,2)/pow(lo_g,2)-gama_g/(gama_g-1.0)*eta_g*pow(lo_g,gama_g-1.0);
			dfun = pow(Q/phi_g,2)/pow(lo_g,3)-gama_g*eta_g*pow(lo_g,gama_g-2.0);
			NewtonRapshon(&lo_g, &err1, fun,dfun,lo_g,eps);
			k=k+1;
		}
    if(k>=it_max)
        printf("RI2U,err:%lf!",err1);
	U->p_g = pow(lo_g,gama_g)*eta_g;
	U->u_g = Q/phi_g/lo_g+u_s;
	U->p_s = (P-Q*(U->u_g-u_s)-phi_g*U->p_g)/phi_s;
	U->lo_g = lo_g;
	U->lo_s = lo_s;
	U->u_s = u_s;
}

//From var U to calculate Riemann invariants.
void U2RI_cal(const struct U_var * U, struct RI_var * RI)
{
	const double eps = config[4];
	const double gama_g = config[6];
	double phi_s = U->phi_s, lo_s = U->lo_s, u_s = U->u_s, p_s = U->p_s, lo_g = U->lo_g, u_g = U->u_g, p_g = U->p_g;
	double phi_g = 1.0-phi_s;
	
    RI->eta_g=p_g/pow(lo_g,gama_g);
    RI->Q=phi_g*lo_g*(u_g-u_s);
    RI->P=phi_g*lo_g*pow(u_g-u_s,2)+phi_g*p_g+phi_s*p_s;
    RI->H=0.5*pow(u_g-u_s,2)+gama_g/(gama_g-1.0)*p_g/lo_g;
    RI->lo_s=lo_s;
    RI->u_s=u_s;
}


//compute primitive var
void primitive_comp(double * U, struct U_var * U_L, struct U_var * U_R, double phi_sL, double phi_sR, double phi_sL_out, double phi_sR_out, double area_L, double area_R)
{   
	double phi_gL=1-phi_sL;
	double phi_gR=1-phi_sR;
	const double gama_g = config[6], gama_s = config[106];
	double eps = config[4];
	double phi_s = area_L*phi_sL+area_R*phi_sR;
	double phi_g = 1.0-phi_s;
	double U1=U[0], U2=U[1], U3=U[3], U4=U[4], U5=U[5], U6=U[7];
	double lo_gR = U1/phi_g;
	double u_gR  = U2/U1;
	double p_gR  = (U3/phi_g - 0.5*lo_gR*pow(u_gR,2))*(gama_g-1.0);
	double lo_s  = U4/phi_s;
	double u_s   = U5/U4;
	double p_sR  = (U6/phi_s - 0.5*lo_s*pow(u_s,2))*(gama_s-1.0);
	U_L->v_g = U[2]/U1;
	U_R->v_g = U_L->v_g;	
	U_L->v_s = U[6]/U4;
	U_R->v_s = U_L->v_s;
	U_L->phi_s = phi_sL;
	U_R->phi_s = phi_sR;
	U_L->phi_g = phi_gL;
	U_R->phi_g = phi_gR;
	double fun[4], dfun[4][4], x_star[4], x_0[4];
	int it_max = 500, k = 0;
	double err2 = 1e50;
	while (k<it_max && err2>eps && abs(phi_sL-phi_sR)>eps)
		{			
			fun[0] = U3-area_L*phi_gL*(0.5*((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL)*pow((U2-area_R*phi_gR*lo_gR*u_gR)/(U1-area_R*phi_gR*lo_gR),2.0)+(p_gR/pow(lo_gR,gama_g)*pow((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL,gama_g))/(gama_g-1.0))-area_R*phi_gR*(0.5*lo_gR*pow(u_gR,2.0)+p_gR/(gama_g-1.0));
			fun[1] = phi_gL*((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL)*(((U2-area_R*phi_gR*lo_gR*u_gR)/(U1-area_R*phi_gR*lo_gR))-u_s)-phi_gR*lo_gR*(u_gR-u_s);
			fun[2] = phi_gL*((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL)*pow(((U2-area_R*phi_gR*lo_gR*u_gR)/(U1-area_R*phi_gR*lo_gR))-u_s,2.0)+phi_gL*(p_gR/pow(lo_gR,gama_g)*pow((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL,gama_g))+phi_sL*(((U6-0.5*phi_s*lo_s*pow(u_s,2.0))*(gama_s-1.0)-area_R*phi_sR*p_sR)/area_L/phi_sL)-phi_gR*lo_gR*pow(u_gR-u_s,2.0)-phi_gR*p_gR-phi_sR*p_sR;
			fun[3] = 0.5*pow(((U2-area_R*phi_gR*lo_gR*u_gR)/(U1-area_R*phi_gR*lo_gR))-u_s,2.0)+gama_g/(gama_g-1.0)*(p_gR/pow(lo_gR,gama_g)*pow((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL,gama_g))/((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL)-0.5*pow(u_gR-u_s,2.0)-gama_g/(gama_g-1.0)*p_gR/lo_gR;
			dfun[0][0] = area_L*phi_gL*((gama_g*p_gR*pow((U1 - area_R*lo_gR*phi_gR)/(area_L*phi_gL),gama_g))/(pow(lo_gR,gama_g + 1.0)*(gama_g - 1.0)) - (area_R*phi_gR*pow(U2 - area_R*lo_gR*phi_gR*u_gR,2.0))/(2.0*area_L*phi_gL*pow(U1 - area_R*lo_gR*phi_gR,2.0)) + (area_R*phi_gR*u_gR*(U2 - area_R*lo_gR*phi_gR*u_gR))/(area_L*phi_gL*(U1 - area_R*lo_gR*phi_gR)) + (area_R*gama_g*p_gR*phi_gR*pow((U1 - area_R*lo_gR*phi_gR)/(area_L*phi_gL),gama_g - 1.0))/(area_L*pow(lo_gR,gama_g)*phi_gL*(gama_g - 1.0))) - (area_R*phi_gR*pow(u_gR,2.0))/2.0;
			dfun[1][0] = - (area_R*phi_gR)/(gama_g - 1.0) - (area_L*phi_gL*pow((U1 - area_R*lo_gR*phi_gR)/(area_L*phi_gL),gama_g))/(pow(lo_gR,gama_g)*(gama_g - 1));
			dfun[2][0] = (area_R*lo_gR*phi_gR*(U2 - area_R*lo_gR*phi_gR*u_gR))/(U1 - area_R*lo_gR*phi_gR) - area_R*lo_gR*phi_gR*u_gR;
			dfun[3][0] = 0.0;
			dfun[0][1] = (area_R*phi_gR*(u_s - (U2 - area_R*lo_gR*phi_gR*u_gR)/(U1 - area_R*lo_gR*phi_gR)))/area_L - ((U1 - area_R*lo_gR*phi_gR)*((area_R*phi_gR*u_gR)/(U1 - area_R*lo_gR*phi_gR) - (area_R*phi_gR*(U2 - area_R*lo_gR*phi_gR*u_gR))/pow(U1 - area_R*lo_gR*phi_gR,2.0)))/area_L - phi_gR*(u_gR - u_s);
			dfun[1][1] = 0.0;
			dfun[2][1] = - lo_gR*phi_gR - (area_R*lo_gR*phi_gR)/area_L;
			dfun[3][1] = 0.0;
			dfun[0][2] = (2*(U1 - area_R*lo_gR*phi_gR)*((area_R*phi_gR*u_gR)/(U1 - area_R*lo_gR*phi_gR) - (area_R*phi_gR*(U2 - area_R*lo_gR*phi_gR*u_gR))/pow(U1 - area_R*lo_gR*phi_gR,2.0))*(u_s - (U2 - area_R*lo_gR*phi_gR*u_gR)/(U1 - area_R*lo_gR*phi_gR)))/area_L - phi_gR*pow(u_gR - u_s,2.0) - (area_R*phi_gR*pow(u_s - (U2 - area_R*lo_gR*phi_gR*u_gR)/(U1 - area_R*lo_gR*phi_gR),2))/area_L - (gama_g*p_gR*phi_gL*pow((U1 - area_R*lo_gR*phi_gR)/(area_L*phi_gL),gama_g))/pow(lo_gR,gama_g + 1.0) - (area_R*gama_g*p_gR*phi_gR*pow((U1 - area_R*lo_gR*phi_gR)/(area_L*phi_gL),gama_g - 1.0))/(area_L*pow(lo_gR,gama_g));
			dfun[1][2] = (phi_gL*pow((U1 - area_R*lo_gR*phi_gR)/(area_L*phi_gL),gama_g))/pow(lo_gR,gama_g) - phi_gR;
			dfun[2][2] = (2.0*area_R*lo_gR*phi_gR*(u_s - (U2 - area_R*lo_gR*phi_gR*u_gR)/(U1 - area_R*lo_gR*phi_gR)))/area_L - lo_gR*phi_gR*(2.0*u_gR - 2.0*u_s);
			dfun[3][2] = - phi_sR - (area_R*phi_sR)/area_L;
			dfun[0][3] = ((area_R*phi_gR*u_gR)/(U1 - area_R*lo_gR*phi_gR) - (area_R*phi_gR*(U2 - area_R*lo_gR*phi_gR*u_gR))/pow(U1 - area_R*lo_gR*phi_gR,2.0))*(u_s - (U2 - area_R*lo_gR*phi_gR*u_gR)/(U1 - area_R*lo_gR*phi_gR)) + (gama_g*p_gR)/(pow(lo_gR,2.0)*(gama_g - 1.0)) - (area_L*pow(gama_g,2.0)*p_gR*phi_gL*pow((U1 - area_R*lo_gR*phi_gR)/(area_L*phi_gL),gama_g))/(pow(lo_gR,gama_g + 1.0)*(U1 - area_R*lo_gR*phi_gR)*(gama_g - 1.0)) - (area_R*pow(gama_g,2.0)*p_gR*phi_gR*pow((U1 - area_R*lo_gR*phi_gR)/(area_L*phi_gL),gama_g - 1.0))/(pow(lo_gR,gama_g)*(U1 - area_R*lo_gR*phi_gR)*(gama_g - 1.0)) + (area_L*area_R*gama_g*p_gR*phi_gL*phi_gR*pow((U1 - area_R*lo_gR*phi_gR)/(area_L*phi_gL),gama_g))/(pow(lo_gR,gama_g)*pow(U1 - area_R*lo_gR*phi_gR,2.0)*(gama_g - 1.0));
			dfun[1][3] = (area_L*gama_g*phi_gL*pow((U1 - area_R*lo_gR*phi_gR)/(area_L*phi_gL),gama_g))/(pow(lo_gR,gama_g)*(U1 - area_R*lo_gR*phi_gR)*(gama_g - 1)) - gama_g/(lo_gR*(gama_g - 1.0));
			dfun[2][3] = u_s - u_gR + (area_R*lo_gR*phi_gR*(u_s - (U2 - area_R*lo_gR*phi_gR*u_gR)/(U1 - area_R*lo_gR*phi_gR)))/(U1 - area_R*lo_gR*phi_gR);
			dfun[3][3] = 0.0;
			x_0[0] = lo_gR;
			x_0[1] = p_gR;
			x_0[2] = u_gR;
			x_0[3] = p_sR;
			NewtonRapshon_matrix(x_star, &err2, fun, dfun[0], x_0, eps);
			lo_gR=fmax(x_star[0],eps);
			lo_gR=fmin(lo_gR,U1/area_R/phi_gR-eps);
			p_gR =fmax(x_star[1],eps);
			u_gR =x_star[2];
			p_sR =fmax(x_star[3],eps);
			p_sR =fmin(p_sR,(U6-0.5*phi_s*lo_s*pow(u_s,2))*(gama_s-1.0)/area_R/phi_sR-eps);
			k=k+1;
		}
	if (k>=it_max)
        printf("RI2U,err:%lf! %lf,%lf,%lf,%lf,%lf,%lf",err2,phi_sL,phi_sR,lo_gR,p_gR,u_gR,p_sR);

	U_L->lo_g = (U1-area_R*phi_gR*lo_gR)/area_L/phi_gL;
	U_R->lo_g = lo_gR;
	U_L->u_g = (U2-area_R*phi_gR*lo_gR*u_gR)/(U1-area_R*phi_gR*lo_gR);
	U_R->u_g = u_gR;
	U_L->p_s = ((U6-0.5*phi_s*lo_s*pow(u_s,2.0))*(gama_s-1)-area_R*phi_sR*p_sR)/area_L/phi_sL;
	U_L->p_g = p_gR/pow(lo_gR,gama_g)*pow((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL,gama_g);
	U_R->p_s = p_sR;
	U_R->p_g = p_gR;	
	U_L->lo_s= lo_s;
	U_R->lo_s= lo_s;
	U_L->u_s = u_s;
	U_R->u_s = u_s;

	U_L->U_lo_g = (1.0-phi_sL_out)*U_L->lo_g;
	U_L->U_u_g  = U_L->U_lo_g*U_L->u_g;
	U_L->U_v_g  = U_L->U_lo_g*U_L->v_g;
	U_L->U_e_g  = U_L->U_lo_g*(U_L->p_g/U_L->lo_g/(gama_g-1.0)+0.5*U_L->u_g*U_L->u_g+0.5*U_L->v_g*U_L->v_g);
	U_L->U_lo_s = phi_sL_out*U_L->lo_s;
	U_L->U_u_s  = U_L->U_lo_s*U_L->u_s;
	U_L->U_v_s  = U_L->U_lo_s*U_L->v_s;
	U_L->U_e_s  = U_L->U_lo_s*(U_L->p_s/U_L->lo_s/(gama_s-1.0)+0.5*U_L->u_s*U_L->u_s+0.5*U_L->v_s*U_L->v_s);
	U_R->U_lo_g = (1.0-phi_sR_out)*U_R->lo_g;
	U_R->U_u_g  = U_R->U_lo_g*U_R->u_g;
	U_R->U_v_g  = U_R->U_lo_g*U_R->v_g;
	U_R->U_e_g  = U_R->U_lo_g*(U_R->p_g/U_R->lo_g/(gama_g-1.0)+0.5*U_R->u_g*U_R->u_g+0.5*U_R->v_g*U_R->v_g);
	U_R->U_lo_s = phi_sR_out*U_R->lo_s;
	U_R->U_u_s  = U_R->U_lo_s*U_R->u_s;
	U_R->U_v_s  = U_R->U_lo_s*U_R->v_s;
	U_R->U_e_s  = U_R->U_lo_s*(U_R->p_s/U_R->lo_s/(gama_s-1.0)+0.5*U_R->u_s*U_R->u_s+0.5*U_R->v_s*U_R->v_s);	
}
