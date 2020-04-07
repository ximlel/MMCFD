#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/var_struc.h"
#include "../include/meshing.h"

static int line_mesh(struct mesh_var * mv, int n_add)
{
	if(isinf(config[10]) || config[10] < config[4])
		{
			fprintf(stderr, "Without a proper spatial grid size!\n");
			exit(2);
		}
	const int num_cell = (int)config[3] + n_add;
	if (num_cell > (int)config[3])
		printf("There are %d ghost cell!\n", mv->num_ghost = num_cell - (int)config[3]);
	
	mv->num_pt = num_cell + 1;
	mv->X = malloc(mv->num_pt * sizeof(double));
	if(mv->X == NULL)
		{
			printf("Not enough memory in 1-D mesh constructed!\n");
			goto return_0;
		}
	int k;
	for(k = 0; k < mv->num_pt; k++)
		mv->X[k] = k * config[10];		
	
	mv->cell_pt = malloc(num_cell * sizeof(void *));
	if(mv->cell_pt == NULL)
		{
			fprintf(stderr, "Not enough memory in 1-D mesh constructed!\n");
			goto return_0;
		}
	for(k = 0; k < num_cell; k++)
		{
			mv->cell_pt[k] = malloc(3 * sizeof(int));
			if(mv->cell_pt[k] == NULL)
				{
					for(int i = 0; i < k; ++i)
						{
							free(mv->cell_pt[i]);
							mv->cell_pt[i] = NULL;
						}
					fprintf(stderr, "Not enough memory in CELL_POINT[%d]!\n", k);
					goto return_0;
				}
			
			mv->cell_pt[k][0] = 2;
			mv->cell_pt[k][1] = k;
			mv->cell_pt[k][2] = k + 1;
		}

	mv->num_border[0] = 1;	
	mv->num_border[1] = 2;	
	mv->border_pt = malloc(2 * sizeof(int));	
	mv->border_cond = malloc(2 * sizeof(int));
	if(mv->border_pt == NULL || mv->border_cond == NULL)
		{
			printf("Not enough memory in square mesh constructed!\n");
			goto return_0;
		}
	
	mv->border_pt[0] = 0;
	mv->border_pt[1] = mv->num_pt - 1;

	return 1;

 return_0:
	free(mv->X);
	mv->X = NULL;	
	free(mv->border_pt);
	mv->border_pt = NULL;
	free(mv->border_cond);
	mv->border_cond = NULL;	
	if (mv->cell_pt != NULL)
		{
			for(int i = 0; i < num_cell; ++i)
				{
					free(mv->cell_pt[i]);
					mv->cell_pt[i] = NULL;
				}
			free(mv->cell_pt);
			mv->cell_pt = NULL;
		}
	return 0;	
}

static void line_border_cond(struct mesh_var * mv, int n_add, int left, int right)
{
	if (left >= 0 || right >= 0)
		{
			fprintf(stderr, "Input wrong boundary condition in quadrilateral mesh!\n");
			exit(2);
		}
	else if ((left == -7) - (right == -7) != 0)
		{
			fprintf(stderr, "Periodic boundary condition error!\n");
			exit(2);
		}

	const int num_cell = (int)config[3] + n_add;

	mv->border_cond[0] = left;
	mv->border_cond[1] = right;
	
	mv->peri_cell = malloc(num_cell * sizeof(int));
	if(mv->peri_cell == NULL)
		{
			printf("Not enough memory in quad periodic boundary constructed!\n");
			exit(5);
		}
	int k;
	for(k = 0; k < num_cell; k++)
		mv->peri_cell[k] = -1;
	
	if (right == -7)				
		mv->peri_cell[0] = num_cell - 1;
	if (left == -7)
		mv->peri_cell[1] = 0;
			
	if(mv->num_ghost > 0)
		period_cell_modi(mv);
}

void free_1D_mesh(struct mesh_var * mv)
{
	const int n_a = 0;
	if (line_mesh(mv, n_a) == 0)
		exit(5);

	line_border_cond(mv, n_a, -3, -3);
}

void inflow_1D_mesh(struct mesh_var * mv)
{
	const int n_a = 0;
	if (line_mesh(mv, n_a) == 0)
		exit(5);

	line_border_cond(mv, n_a, -1, -3);
}

void periodic_1D_mesh(struct mesh_var * mv)
{
	const int n_a = 2;
	if (line_mesh(mv, n_a) == 0)
		exit(5);

	line_border_cond(mv, n_a, -7, -7);
}
