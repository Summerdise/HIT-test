/**************************************************
 Mex C Code for reproducing experiment 1
***************************************************/
# include <stdio.h>
# include <stdlib.h>
# include "mex.h"        /* the algorithm is connect to matlab */
# include "math.h"
# define PI 3.1415926
# define ABS(x) ((x)>0? (x):(-(x)))
# define MAX(x, y) ((x)>(y)? (x):(y))
# define MIN(x, y) ((x)<(y)? (x):(y))

double *double_vector(int n)
{
    double *v; 
    
    v = (double*) mxCalloc (n, sizeof(double));
    
    return v; 
}

int *int_vector(int n)
{
    int *v; 
    
    v = (int*) mxCalloc (n, sizeof(int));
    
    return v; 
}

double **double_matrix(int m, int n)
{
    double **mat; 
    int i; 
    
    mat = (double**) mxCalloc(m, sizeof(double*)); 
    for (i=0; i<m; i++)
        mat[i] = double_vector(n); 
    
    return mat; 
}

int **int_matrix(int m, int n)
{
    int **mat; 
    int i; 
    mat = (int**) mxCalloc(m, sizeof(int*)); 
    for (i=0; i<m; i++)
        mat[i] = int_vector(n); 
    return mat; 
}


/* getting the index of matlab image, (x,y) location, (sx, sy) sizes */
int px(int x, int y, int bx, int by)
{
   return (x + (y-1)*bx - 1);
}

 /* variables */
 int n;                      /* number of images */
 int N;                      /* number of orientations */
 double **fI;                /* filtered images */
 int h;                      /* halfsizes of filters */
 int sx, sy;                 /* sizes of image */ 
 int binnum;                  /* number of bins to keep the pooled responses */
 double binsize;              /* length of each bin */
 double *histog;              /* histogram of natural images */
 double Upperbound;
 
void Ctransform()
{
   int x, y, here, ind, i, ii; 
   
   for (x=1; x<=sx; x++)
       for (y=1; y<=sy; y++)
       {
        here = px(x, y, sx, sy); 
        for (ind=0; ind<N; ind++)
           {    
             for (i=0; i<n; i++)
             {
                ii = ind*n+i; 
                fI[ii][here] = Upperbound*(2./(1.+exp(-2.*fI[ii][here]/Upperbound))-1.); 
             }
            }
       }
}

void Cmp()
{
   int ind, i, x, y, here, b, tot; 
   
   tot = 0; 
   for (x=h+1; x<sx-h; x++)
      for (y=h+1; y<sy-h; y++)
       {
        here = px(x, y, sx, sy); 
        for (ind=0; ind<N; ind++)
           {    
             for (i=0; i<n; i++)
             {
                 b = MIN(floor(fI[ind*n+i][here]/binsize), binnum-1); 
                 histog[b] += 1.;
                 tot ++;
             }
        }
      }
   for (b=0; b<binnum; b++)
       histog[b] /= (binsize*tot);
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
 
 h = floor(mxGetScalar(prhs[c++])+.5);    /* half size of filters */
 
 sx = floor(mxGetScalar(prhs[c++])+.5);
 sy = floor(mxGetScalar(prhs[c++])+.5);   /* size of images */
 
 binsize = mxGetScalar(prhs[c++]); 
 binnum = floor(mxGetScalar(prhs[c++])+.5); 
 histog = mxGetPr(prhs[c++]);   
 Upperbound = mxGetScalar(prhs[c++]); 
 
 Ctransform();
 Cmp();
}

