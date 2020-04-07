#include <stdio.h>
#include <stdlib.h>
#include <math.h>


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
