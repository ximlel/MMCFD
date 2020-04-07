#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/var_struc.h"



#define CV_COPY(var)  cv->var[i] = cv->var[pc[i]]
#define FV_COPY(var)  FV->var[i] = FV->var[pc[i]]

void period_ghost(struct cell_var * cv, struct mesh_var mv, struct flu_var * FV, double t)
{
	const int dim = (int)config[0];
	const int order = (int)config[9];
	const int num_cell = (int)config[3];
	const int num_cell_ghost = mv.num_ghost + num_cell;
	const int *pc = mv.peri_cell;	
	
	for(int i = num_cell; i < num_cell_ghost; i++)
		{
			CV_COPY(U_rho);			
			CV_COPY(U_gamma);
			CV_COPY(U_e);
			CV_COPY(U_u);
			FV_COPY(RHO);
			FV_COPY(gamma);
			FV_COPY(P);
			FV_COPY(U);
			if (order > 1)
				{
					CV_COPY(gradx_rho);
					CV_COPY(gradx_e);
					CV_COPY(gradx_u);
				}
			if (dim > 1)
				{
					CV_COPY(U_v);
					FV_COPY(V);
					if (order > 1)
						{
							CV_COPY(grady_rho);
							CV_COPY(grady_e);
							CV_COPY(grady_u);
							CV_COPY(grady_v);
							CV_COPY(gradx_v);
						}
				}
			if (dim > 2)
				{
					CV_COPY(U_w);
					FV_COPY(W);
					if (order > 1)
						{
							CV_COPY(gradz_rho);
							CV_COPY(gradz_e);
							CV_COPY(gradz_u);
							CV_COPY(gradz_v);
							CV_COPY(gradz_w);
							CV_COPY(grady_w);
							CV_COPY(gradx_w);
						}
				}
			if ((int)config[2] == 2)
				{
					CV_COPY(U_phi);
					FV_COPY(PHI);
					if (order > 1)
						{
							CV_COPY(gradx_phi);
							if (dim > 2)
								CV_COPY(grady_phi);
							if (dim > 3)
								CV_COPY(gradz_phi);
						}
					CV_COPY(U_e_a);
					FV_COPY(Z_a);
					if (order > 1)
						{
							CV_COPY(gradx_z_a);
							if (dim > 2)
								CV_COPY(grady_z_a);
							if (dim > 3)
								CV_COPY(gradz_z_a);
						}
				}
			CV_COPY(U_e_a);
		}
}



void period_cell_modi(struct mesh_var * mv)
{
	const int num_cell = mv->num_ghost + (int)config[3];
	int *pc = mv->peri_cell;
	
	int per_num[num_cell], per_n = 0;
	int i, j;
	for (i = 0; i < num_cell; i++)
		{
			if (pc[i] >= 0)
				per_n++;
			per_num[i] = per_n;
		}

	for (i = 0; i < num_cell; i++)
		if (pc[i] >= 0)
			{
				pc[i] -= per_num[pc[i]];
				if (pc[i] < 0)
					pc[i] = 0;
			}

	int *cc_tmp, pc_tmp; 
	for (i = 1; i < num_cell; i++)
		{
			j = i;
			for(j = i; pc[j-1] >= 0 && j >= 1; j--)
				{
					if(pc[j] < 0)
						{
							pc_tmp = pc[j-1];
							pc[j-1] = pc[j];
							pc[j] = pc_tmp;
							cc_tmp = mv->cell_pt[j];
							mv->cell_pt[j] = mv->cell_pt[j-1];
							mv->cell_pt[j-1] = cc_tmp;
						}
					else
						break;
				}	 				
		}
	
	mv->bc = period_ghost;
}
