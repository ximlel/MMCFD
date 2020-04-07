/*!memory_management
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include <unistd.h>
#include <sys/stat.h>


void DispPro(double pro, int step)
{
        for (int k=0; k<77; k++)  
                putchar('\b');//将当前行全部清空，用以显示最新的进度条状态
        for (int j=0; j<lround(pro/2); j++)  
                putchar('+');//打印进度条上已经完成的部分，用‘+’表示  
        for (int j=1; j<=50-lround(pro/2); j++)  
                putchar('-');//打印进度条上还有多少没有完成的  
        fprintf(stdout, "  %6.2f%%   STEP=%-8d", pro, step);  
        fflush(stdout);
}

int CreateDir(const char * pPath)
{
	if(-1 != access(pPath,0))
		return -1;

	char tmpPath[FILENAME_MAX];
	const char* pCur = pPath;

	memset(tmpPath,0,sizeof(tmpPath));
    
	int pos=0;
    while(*pCur++!='\0')
		{
			tmpPath[pos++] = *(pCur-1);

			if(*pCur=='/' || *pCur=='\0')
				{
					if(0!=access(tmpPath,0)&&strlen(tmpPath)>0)
						{
							mkdir(tmpPath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
						}
				}
		}
	if(!access(pPath,0))
		return 0;
	else
		return 1;
}


void init_mem(double * p[], int n, int * cell_pt[])
{
	for(int k = 0; k < n; ++k)
		{
			p[k] = calloc(cell_pt[k][0], sizeof(double));
			if(p[k] == NULL)
				{
					printf("Initialize memory fail!\n");
					for(int j = 0; j < k; j++)
						{
							free(p[k]);
							p[k] = NULL;
						}
					exit(5);
				}
		}
}

void init_mem_int(int * p[], int n, int * cell_pt[])
{
	for(int k = 0; k < n; ++k)
		{
			p[k] = malloc(cell_pt[k][0] * sizeof(int));
			if(p[k] == NULL)
				{
					printf("Initialize memory fail!\n");
					for(int j = 0; j < k; j++)
						{
							free(p[k]);
							p[k] = NULL;
						}
					exit(5);
				}
		}
}
