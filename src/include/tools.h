#ifndef TOOLS_H
#define TOOLS_H

/*!\file tools.h
 * \brief  Some mathematical algorithm functions.
 * \author Lei Xin
 */

void init_mem(double * p[], int n, int * cell_pt[]);
void init_mem_int(int * p[], int n, int * cell_pt[]);
void DispPro(double pro, int step);
int CreateDir(const char* pPath);

/*!\brief A function to caculate the inverse of the input square matrix.
 * \param a The pointer of the input square matrix.
 * \param[in] n The order of the input square matrix.
 */

void mat_mul(double A[], double B[], double C[], int m, int p, int n);
void mat_add(double A[], double B[], double C[], int m, int n);
void mat_sub(double A[], double B[], double C[], int m, int n);
int rinv(double a[], int n);
void Gauss_elimination(int n, double (*a)[n+1], double *x);
double rnd( double *r);

#endif
