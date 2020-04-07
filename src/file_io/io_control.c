#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <dirent.h>

#include "../include/var_struc.h"
#include "../include/tools.h"
#include "../include/file_io.h"



static void config_check()
{
	const int dim = (int)config[0];

	if(isinf(config[1]) && isinf(config[5]))
		{
			fprintf(stderr, "The total time or the maximum number of time steps must be setted!\n");
			exit(2);
		}

	// species
	config[2] = isinf(config[2]) ? 1 : config[2];	

	if(isinf(config[4]))
		config[4] = EPS;
	else if(config[4] < 0.0 || config[4] > 0.01)
		{
			fprintf(stderr, "eps(%lf) should in (0, 0.01)!\n", config[4]);
			exit(2);
		}
	
	if(isinf(config[6]))
		config[6] = 1.4;	
	else if(config[6] < (1.0 + config[4]))
		{
			fprintf(stderr, "The constant of the perfect gas(%lf) should be larger than 1.0!\n", config[6]);
			exit(2);
		}

	if (isinf(config[7]))
		{
			if (dim == 1)
				config[7] = 0.9;
			else if (dim == 2)
 				config[7] = 0.45;
		}

	// EUL/LAG
	config[8] = isinf(config[8]) ? 0 : config[8];
	// reconstruction(prim_var/cons_var)
	config[31] = isinf(config[31]) ? 0 : config[31];
	// v_fix
	config[61] = isinf(config[61]) ? 0 : config[61];
}



void example_io(const char *example, char *add_mkdir, const int i_or_o)
{
	const int dim = (int)config[0];
	const int el = (int)config[8];
	const int order = (int)config[9];

	char tmp[15];
	
	if (i_or_o == 0)
		{
			switch (dim)
				{
				case 1 :
					strcpy(add_mkdir, "../data_out/one-dim/");
					break;
				case 2 :
					strcpy(add_mkdir, "../data_out/two-dim/");
					break;
				case 3 :
					strcpy(add_mkdir, "../data_out/three-dim/");
					break;
				default :
					fprintf(stderr, "Strange computational dimension!\n");
					exit(2);
					break;
				}
			switch (el)
				{
				case 0 :
					strcat(add_mkdir, "EUL_");
					break;
				case 1 :
					strcat(add_mkdir, "LAG_");
					break;
				case 2 :
					strcat(add_mkdir, "ALE_");
					break;
				default :
					fprintf(stderr, "Strange description method of fluid motion!\n");
					exit(2);
					break;
				}
			sprintf(tmp, "%d_order/", order);
			strcat(add_mkdir, tmp);
		}
	else
		{		
			switch (dim)
				{
				case 1 :
					strcpy(add_mkdir, "../data_in/one-dim/");
					break;
				case 2 :
					strcpy(add_mkdir, "../data_in/two-dim/");
					break;
				case 3 :
					strcpy(add_mkdir, "../data_in/three-dim/");
					break;
				default :
					fprintf(stderr, "Strange computational dimension!\n");
					exit(2);
					break;
				}
		}
	strcat(add_mkdir, example);
	
	DIR * dir_test = NULL;
	dir_test = opendir(add_mkdir);
	if (dir_test == NULL)
		{
			if (i_or_o == 0)
				{
					if(CreateDir(add_mkdir))
						{
							fprintf(stderr, "Output directory '%s' construction failed.\n", add_mkdir);
							exit(1);
						}
					else
						printf("Output directory '%s' constructed.\n", add_mkdir);				
				}
			else
				{
					fprintf(stderr, "Input directory is not exist!\n");
					exit(1);
				}
		}
	closedir(dir_test);
	strcat(add_mkdir, "/");
}



#define STR_FLU_INI(sfv)										\
	do {														\
		strcpy(add, add_mkdir);									\
		strcat(add, #sfv ".txt");								\
		FV.sfv = malloc((int)config[3] * sizeof(double));		\
		if(FV.sfv == NULL)										\
			{													\
				fprintf(stderr, "Initialize memory fail!\n");	\
				exit(5);										\
			}													\
		r = flu_var_init(add, FV.sfv, 1) ? r : 0;				\
	} while(0)

struct flu_var flu_conf_load(const char *example)
{
	const int dim = (int)config[0];
	struct flu_var FV = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

	int r = 1;
	
	// Find the directory of input data.
	char add_mkdir[FILENAME_MAX];
	example_io(example, add_mkdir, 1);

	// We read the initial data file. 
	char add[FILENAME_MAX];
	strcpy(add, add_mkdir);
	strcat(add, "config.txt");
	configurate(add);

	printf("%s is configurated, dimension = %d .\n", example, dim);

	strcpy(add, add_mkdir);								
	strcat(add, "P.txt");	
	flu_var_init(add, config, 0); // pre read.
		
	STR_FLU_INI(P);
	STR_FLU_INI(RHO);
	STR_FLU_INI(U);	
	if (dim > 1)
		STR_FLU_INI(V);
	else if(!isinf(config[30]))
		{		
			FV.V = calloc((int)config[3], sizeof(double));
			if(FV.V == NULL)									
				{												
					fprintf(stderr, "Initialize memory fail!\n");		
					exit(5);							
				}
		}
	if (dim > 2)		
		STR_FLU_INI(W);
	else if(!isinf(config[30]))
		{		
			FV.W = calloc((int)config[3], sizeof(double));
			if(FV.W == NULL)									
				{												
					fprintf(stderr, "Initialize memory fail!\n");		
					exit(5);							
				}
		}
	if ((int)config[2] == 2)
		{		 							
			STR_FLU_INI(PHI);
			STR_FLU_INI(Z_a);
		}
	if (r == 0)
		{
			printf("Initial date of fluid field failed to read!\n");
			exit(1);
		}

	STR_FLU_INI(gamma);
	if (r == 0)
		for(int i = 0; i < (int)config[3]; i++)
			FV.gamma[i] = 1.0+1.0/(FV.Z_a[i]/(config[6]-1.0)+(1.0-FV.Z_a[i])/(config[106]-1.0));
	
	printf("%s is initialized, grid number = %d .\n", example, (int)config[3]);

	config_check();
	return FV;
}
