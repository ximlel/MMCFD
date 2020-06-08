#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>


void mat_add(double A[], double B[], double C[], int m, int n) 
{
	int i,j;
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			C[i*n+j] = A[i*n+j] + B[i*n+j];
}

void mat_sub(double A[], double B[], double C[], int m, int n)
{
	int i, j;
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			C[i*n+j] = A[i*n+j] - B[i*n+j];
}

void mat_mul(double A[], double B[], double C[], int m, int p, int n)
{
	int i,j,k;
	double A0[m][n], B0[n][p];
	for (i = 0; i < m; i++)
		for (k = 0; k < n; k++)
			A0[i][k] = A[i*n+k];
	for (k = 0; k < n; k++)
		for (j = 0; j < p; j++)
			B0[k][j] = B[k*p+j];
	for (i = 0; i < m; i++)
		for (j = 0; j < p; j++)
			{
				C[i*p+j] = 0.0;
				for (k = 0; k < n; k++)
					C[i*p+j] += A0[i][k] * B0[k][j];
			}
}

int rinv(double a[], int n)
{
    int *is,*js,i,j,k,l,u,v;
    double d,p;
    is=malloc(n*sizeof(int));
    js=malloc(n*sizeof(int));
    for (k=0; k<=n-1; k++)
		{
			d=0.0;
			for (i=k; i<=n-1; i++)
				for (j=k; j<=n-1; j++)
					{
						l=i*n+j;
						p=fabs(a[l]);
						if (p>d)
							{
								d=p;
								is[k]=i;
								js[k]=j;
							}
					}
			if (d+1.0==1.0)
				{
					free(is);
					free(js);
					fprintf(stderr, "Error: no inverse matrix!\n");
					return(0);
				}
			if (is[k]!=k)
				for (j=0; j<=n-1; j++)
					{
						u=k*n+j;
						v=is[k]*n+j;
						p=a[u];
						a[u]=a[v];
						a[v]=p;
					}
			if (js[k]!=k)
				for (i=0; i<=n-1; i++)
					{
						u=i*n+k;
						v=i*n+js[k];
						p=a[u];
						a[u]=a[v];
						a[v]=p;
					}
			l=k*n+k;
			a[l]=1.0/a[l];
			for (j=0; j<=n-1; j++)
				if (j!=k)
					{
						u=k*n+j;
						a[u]=a[u]*a[l];
					}
			for (i=0; i<=n-1; i++)
				if (i!=k)
					for (j=0; j<=n-1; j++)
						if (j!=k)
							{
								u=i*n+j;
								a[u]=a[u]-a[i*n+k]*a[k*n+j];
							}
			for (i=0; i<=n-1; i++)
				if (i!=k)
					{
						u=i*n+k;
						a[u]=-a[u]*a[l];
					}
		}
    for (k=n-1; k>=0; k--)
		{
			if (js[k]!=k)
				for (j=0; j<=n-1; j++)
					{
						u=k*n+j;
						v=js[k]*n+j;
						p=a[u];
						a[u]=a[v];
						a[v]=p;
					}
			if (is[k]!=k)
				for (i=0; i<=n-1; i++)
					{
						u=i*n+k;
						v=i*n+is[k];
						p=a[u];
						a[u]=a[v];
						a[v]=p;
					}
		}
    free(is); free(js);
    return(1);
}

void Gauss_elimination(int n, double (*a)[n+1], double *x)
{ 
	int i,j,k;
	double temp,s,l;

	for(i=0;i<n-1;i++)
		{
			//selected principal element of columns
			k=i;     
			for(j=i+1;j<n;j++)
				{
					if(fabs(a[j][i])>fabs(a[k][i]))
						k=j;
				}			
			//line feed   
			if(k!=i)
				for(j=i;j<=n;j++)
					{
						temp=a[i][j];
						a[i][j]=a[k][j];
						a[k][j]=temp;
					}
			//elimination
			for(j=i+1;j<n;j++)
				{ 
					l=a[j][i]/a[i][i];
					for(k=0;k<n+1;k++)
						a[j][k]=a[j][k]-a[i][k]*l;      
				}			
		}
	//back substitution
	x[n-1]=a[n-1][n]/a[n-1][n-1];
	for(i=n-2;i>=0;i--)
		{
			s=0.0;
			for(j=i;j<n;j++)
				{
					if(j==i)
						continue;
					s+=a[i][j]*x[j];					
				}
			x[i]=(a[i][n]-s)/a[i][i];
		}
}

void Lagrangian_with_Multiplier()
{
    const int N = 4, M = 3;
    int i,j,l;
    /* H_k = D2_xx_L, N_k = D_h */
    gsl_matrix *D2_xx_L_c_k = gsl_matrix_alloc(N, N), *H_k = gsl_matrix_alloc(N, N);
    gsl_matrix *N_k = gsl_matrix_alloc(N, M);
    int sign = 0;
    gsl_matrix *H_c_NN_k_inv = gsl_matrix_alloc(N, N);
    /* Temp Matrix */
    gsl_matrix *M_tmp_M_M = gsl_matrix_calloc(M, M), *M_tmp_M_N = gsl_matrix_calloc(M, N);
    /* Temp Vector*/
    gsl_vector *V_tmp_M = gsl_vector_calloc(M);
    
    gsl_vector *D_x_L_c_k = gsl_vector_alloc(N);
    gsl_vector *D_L = gsl_vector_alloc(N+M);
    gsl_vector *D_f = gsl_vector_alloc(N);
    gsl_vector *lambda_k = gsl_vector_alloc(M), *lambda_k_b = gsl_vector_alloc(M);
    gsl_vector *x_k = gsl_vector_alloc(N), *x_k_b = gsl_vector_alloc(N);
    gsl_vector *d_k = gsl_vector_alloc(N), *h_k = gsl_vector_alloc(N), *c_k_h = gsl_vector_alloc(N);
    gsl_permutation *per = gsl_permutation_alloc(N);
    double c_k, vareps_k, omega_k, m_k;
    int m_idx;
    double L_c_k, L_c_k_beta, ddot;
    const double gamma = 0.5, r = 2, beta = 0.5, sigma = 0.25;
    
    for (i = 0; i < N; i++)
	for (j = 0; j < N; j++)
	    gsl_matrix_set(D2_xx_L_c_k, i, j, 10086);
    for (i = 0; i < M; i++)
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
    gsl_blas_dsyrk(CblasUpper, CblasNoTrans, c_k, N_k, 1.0, H_k);
    /* Inverse Matrix*/ 
    gsl_linalg_LU_decomp(H_k, per, &sign);
    gsl_linalg_LU_invert(H_k, per, H_c_NN_k_inv);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, N_k, H_c_NN_k_inv, 0.0, M_tmp_M_N);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, M_tmp_M_N,  N_k, 0.0, M_tmp_M_M);
    gsl_blas_dgemv(CblasNoTrans, -1.0, M_tmp_M_N, D_f, 0.0, V_tmp_M);
    gsl_vector_add(V_tmp_M, h_k);
    gsl_linalg_LU_decomp(M_tmp_M_M, per, &sign);
    gsl_linalg_LU_solve(M_tmp_M_M, per, V_tmp_M, lambda_k_b);
    gsl_vector_sub(lambda_k_b, c_k_h);
    
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
    gsl_matrix_free(H_c_NN_k_inv);
    gsl_matrix_free(M_tmp_M_N);
    gsl_matrix_free(M_tmp_M_M);
    gsl_vector_free(V_tmp_M);
    gsl_vector_free(D_x_L_c_k);
    gsl_vector_free(D_L);
    gsl_vector_free(D_f);
    gsl_vector_free(lambda_k); gsl_vector_free(lambda_k_b);
    gsl_vector_free(x_k); gsl_vector_free(x_k_b);
    gsl_vector_free(d_k); gsl_vector_free(h_k); gsl_vector_free(c_k_h);
    gsl_permutation_free(per);
}
