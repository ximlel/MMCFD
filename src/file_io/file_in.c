/*!\file file_in.c
 * \brief This file is a collection of functions used to read files
 and check whether they are qualified. 
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "../include/var_struc.h"



/*!\brief This function counts how many numbers are there in the initial data file. 
 * \param[in] fp The pointer of the file to read in.
 * \return The number of the numbers in the initial data file.
 */
static int flu_var_count(FILE * fp, const char * add)
{	
	int num = 0; // Data number.
	
	// "flg" helps us to count.	
	// - flase: when we read a number-using character (0, 1, 2, ..., e, E, minus sign and dot).
	// -  true: when we read a non-number-using character.
	int flg = 0;

	int ch;


	while ((ch = getc(fp)) != EOF)
		{		
			// Count the data number.
			if (ch == 45 || ch == 46 || ch == 69 || ch == 101 || isdigit(ch))
				flg = 1;
			else if (!isspace(ch))
				{
					fprintf(stderr, "Input contains illegal character(ASCII=%d) in the file '%s'!\n", ch, add);
					return 0;
				}
			else if (flg)
				{
					num++;
					flg = 0;
				}				
		}
	
	rewind(fp);

	return num;
}



#define RANGE_WRONG														\
	do {																\
		fprintf(stderr,"Data range is irregular in the file '%s'!\n", add); \
		return 0;														\
	} while(0)

/*!\brief This function counts how many numbers are there in the initial data file. 
 * \param[in] fp The pointer of the file to read in.
 */
static int flu_var_read(FILE * fp, const char * add, double * F)
{
	const int num_cell = (int)config[3];
	const double n_x = config[13], n_y = config[14], n_z = config[15];

	static int D = 0;
	
	if (D == 0 && !isinf(n_x))
		{
			if (n_x < 1.0 || n_y < 1.0 || n_z < 1.0)
				{
					fprintf(stderr, "For the structural mesh, data number in some dimension < 1!\n");
					return 0;
				}
			// Test whether the total number of data is matched.
			else if (num_cell != (int)n_x * (isinf(n_y) ? 1 : (int)n_y) * (isinf(n_z) ? 1 : (int)n_z))
				{
					fprintf(stderr, "Data number is't matched to the structural mesh!\n");
					return 0;
				}
			D = isinf(n_y) ? 1 : (isinf(n_z) ? 2 : 3);
			if (D != (int)config[0])
				{
					fprintf(stderr, "Dimension is't matched to the structural mesh!\n");
					return 0;
				}
		}		

	int r_count = 0, c_count = 0;
	div_t r_div, c_div;

	
	int num = 0;

	int idx = 0;
	char flo_num[50]; // String to store floating numbers.
	char *endptr;
	
	int flg = 0;

	int ch;


	while ((ch = getc(fp)) != EOF)
		{			
			// Read the data number.
			if (ch == 45 || ch == 46 || ch == 69 || ch == 101 || isdigit(ch))
				{							
					flo_num[idx++] = (char)ch;
					flg = 1;
				}
			else if (!isspace(ch))
				{
					fprintf(stderr, "Input contains illegal character(ASCII=%d) in the file '%s'!\n", ch, add);
					return 0;
				}
			else if (flg)
				{
					flo_num[idx] = '\0';
					idx = 0;
					if (num >= num_cell)
						{
							num++;
							break;
						}											
					F[num++] = strtod(flo_num, &endptr);
					if (strlen(endptr) != 0)
						{
							fprintf(stderr,"Reading Sth. that isn't a floating number in the file '%s'!\n", add);
							return 0;
						}
					flg = 0;
				}

			// Test whether the data range is regular and matched to the structual mesh.
			if (ch == '\n' && D > 1)
				{
					r_div = div(num, (int)n_x);
					c_div = div(num, (int)n_x * (int)n_y);
					if (r_div.rem != 0)
						RANGE_WRONG;
					else if (r_div.quot == (r_count+1))								
						r_count = r_div.quot;
					else if (r_div.quot != r_count)
						RANGE_WRONG;
					else if (D == 2)
						; 
					else if (c_div.rem == 0 && (c_div.quot == c_count || c_div.quot == (c_count+1)))
						c_count = c_div.quot;
					else						
						RANGE_WRONG;
				}		
		}
	
	// Test whether the total number of data is matched.
	if (num != num_cell)
		{
			fprintf(stderr, "Data number isn't equal to the given total number in the file '%s'!\n", add);
			return 0;
		}

	return 1;
}



// r_or_c: 0 only count, 1 count and read.
int flu_var_init(const char * add, double * F, const int r_or_c)
{
	FILE * fp;

	// Open the initial data file.
	if ((fp = fopen(add, "r")) == NULL)
		{
			perror(add);
			return 0;
		}
	
	if (isinf(config[3]))		
		config[3] = (double)flu_var_count(fp, add);
	else if (config[3] < 1.0)
		{					
			fprintf(stderr, "Error in the 3-th value of the configuration!\n");
			fclose(fp);
			exit(2);
		}
	else if (F == NULL)
		{					
			fprintf(stderr, "Warning: the incoming pointer is NULL!\n");
			fclose(fp);
			return 0;
		}
	else if (r_or_c)				   		
		if (flu_var_read(fp, add, F) == 0)
			{
				fclose(fp);
				exit(2);							
			}

	fclose(fp);

	return 1;
}



static int config_read(FILE * fp)
{	
	char one_line[200]; // String to store one line.
	char *endptr;

	int i; // Index of config[*].
	
	double tmp;

	while (fgets(one_line, sizeof(one_line), fp) != NULL)
		{
			// A line that doesn't begin with digits is a comment.
			i =strtol(one_line, &endptr, 10);
			for ( ;isspace(*endptr) ; endptr++) ;

			// If the value of config[i] doesn't exit, it is 0 by default.
			if (0 < i && i < N_CONF)
				{
					tmp = strtod(endptr, NULL);
					if(fabs(config[i] - tmp) > EPS)
						CONF_INI(i,tmp);
					config[i] = tmp;
				}
			else if (i != 0 || (*endptr != '#'&& *endptr != '\0'))			   
				fprintf(stderr, "Warning: unknown row occurrs in configuration file!\n");
		}
	if (ferror(fp))
		{
			fprintf(stderr, "Read error occurrs in configuration file!\n");
			return 0;
		}
	return 1;
}



int configurate(char * add)
{
	FILE * fp;

	//open the configuration data file.
	if((fp = fopen(add, "r")) == NULL)
		{
			perror(add);
			return 0;
		}
	
	if(config_read(fp) == 0)
		{
			fclose(fp);
			exit(2);
		}
	fclose(fp);
	return 1;
}
