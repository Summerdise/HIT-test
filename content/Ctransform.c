/**************************************************
 Mex C Code for computing the sigmoid transformation
***************************************************/
# include <stdio.h>
# include <stdlib.h>
# include "mex.h"        /* the algorithm is connect to matlab */
# include "math.h"
# define PI 3.1415926
# define ABS(x) ((x)>0? (x):(-(x)))
# define MAX(x, y) ((x)>(y)? (x):(y))
# define MIN(x, y) ((x)<(y)? (x):(y))

/* getting the index of matlab image, (x,y) location, (sx, sy) sizes */
int px(int x, int y, int bx, int by)    
{            
   return (x + (y-1)*bx - 1); 
 }

 /* variables */
 int n;                      /* number of images */
 int N;                      /* number of orientations */
 double **fI;                /* filtered images */
 double *Sx1, *Sy1;          /* sizes of image */  
 double Upperbound; 

 
 
void Cmp()
{
   
   int x, y, here, ind, i, ii; 
   int sx,sy;
   for (i=0; i<n; i++)
   {
	   sx = Sx1[i]; sy = Sy1[i];
		for (x=1; x<=sx; x++)
			for (y=1; y<=sy; y++)
			{
				here = px(x, y, sx, sy); 
				for (ind=0; ind<N; ind++)
				{    
					ii = ind*n+i; 
					fI[ii][here] = Upperbound*(2./(1.+exp(-2.*fI[ii][here]/Upperbound))-1.); 
				}
            }
   }
}
   


/* mex function is used to pass on the pointers and scalars from matlab, 
   so that heavy computation can be done by C, which puts the results into 
   some of the pointers. After that, matlab can then use these results. 
   
   So matlab is very much like a managing platform for organizing the 
   experiments, and mex C is like a work enginee for fast computation. */

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])                
{
 int ind, i, c; 
 mxArray *f;  
 
 c = 0; /* counter for input pointers and scalars */
 n = floor(mxGetScalar(prhs[c++])+.5);  /* number of images */
 N = floor(mxGetScalar(prhs[c++])+.5);  /* number of orientations */
 fI = mxCalloc(n*N, sizeof(double*));   /* fitered images */
 for (i=0; i<n; i++)
   {
     for (ind=0; ind<N; ind++)
      {  
       f = mxGetCell(prhs[c], ind*n+i); 
       fI[ind*n+i] = mxGetPr(f);    /* get pointers to filtered images */       
      }
    }
 c++; 
 
 Sx1 = mxGetPr(prhs[c++]);
 Sy1 = mxGetPr(prhs[c++]);   /* size of images */
 
 Upperbound = mxGetScalar(prhs[c++]);  
 
 Cmp();
}

     



 

                    