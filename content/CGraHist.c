# include <stdio.h>
# include <stdlib.h>
# include "mex.h"        /* the algorithm is connect to matlab */
# include "math.h"
# define PI 3.1415926
# define BIG_NUMBER 1000
# define ABS(x) ((x)>0? (x):(-(x)))
# define MAX(x, y) ((x)>(y)? (x):(y))
# define MIN(x, y) ((x)<(y)? (x):(y))

/*
* Local Maximization.
* Input: Response map.
* Output: Local-maxed response map.
*/

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
 double **hI;				 /* histogram map, same dimention as fI */
 double *varMap;
 int sx, sy;                 /* current size of image */
 int ax, ay;                 /* half size of bounding box of learned model */ 
 double Upperbound;          /* satuation level for response */
 double *Sx1, *Sy1;           /* sizes of images */
 int h;
 int sub;					/* step of sub-sampling */

/* Function graHist() calculates gradient histogram at rectangular area 
   specified by (top,bottom,left,right).
 Indices "top,bottom,left,right" all start from 1.
 "gh" is a vector of length N, as output.*/
double graHist(int i, int top, int bottom, int left, int right, double* gh)
{
	int ind,x,y,here;
	double sum;
	for(ind = 0; ind < N; ind++)
	{
		for(x=top;x<=bottom;x++)
			for(y=left;y<=right;y++)
			{
				here = px(x,y,sx,sy);
				gh[ind] += fI[ind*n+i][here];
			}
	}
	/* normalize the gradient histogram */
	sum = 0;
	for(ind = 0; ind < N; ind++)
		sum += gh[ind];
	if(sum>0)
		for(ind = 0; ind < N; ind++)
			gh[ind] /= sum;
}

void calcGraHistMap(i)
{
	int ind,x,y,top,bottom,left,right;
	double* v = double_vector(N);
	for(x=1+h+ax;x<sx-ax-h;x+=sub)
		for(y=1+ay+h;y<sy-ax-h;y+=sub)
		{
			top = x-ax;
			bottom = x+ax;
			left = y-ay;
			right = y+ay;
			graHist(i,top,bottom,left,right,v);
			for(ind = 0; ind < N; ind++)
				hI[ind*n+i][px(x,y,sx,sy)] = v[ind];
		}
}

void calcVarMap()
{
	int x,y,i,iOri;
	double sum;
	double** hM; /* n gradient histograms, for temperary storage */
	double* hV; /* temperarily store the mean gradient histogram */
	hM = double_matrix(n,N);
	hV = double_vector(N);

	for(x=1+h+ax;x<sx-ax-h;x+=sub)
		for(y=1+ay+h;y<sy-ax-h;y+=sub)
		{
			for(i=1;i<n;i++)
				for(iOri=0;iOri<N;iOri++)
					hM[i][iOri] = 0;
			for(iOri=0;iOri<N;iOri++)
				hV[iOri] = 0;
			/* Compute mean histogram. */
			for(i=0;i<n;i++)
			{
				for(iOri=0;iOri<N;iOri++)
				{
					hM[i][iOri] = hI[iOri*n+i][px(x,y,sx,sy)];
					hV[iOri] += hM[i][iOri];
				}
			}
			for(iOri=0;iOri<N;iOri++)
			{
				hV[iOri] /= n;
			}
			/* Compute variance. */
			sum = 0;
			for(i=0;i<n;i++)
			{
				/* Euclidean distance */
				for(iOri=0;iOri<N;iOri++)
				{
					sum += (hV[iOri]-hM[i][iOri])*(hV[iOri]-hM[i][iOri]);
				}
			}
			varMap[px(x,y,sx,sy)] = sum/n;
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
 int i, ind, j, c; 
 mxArray *f;  
 
 c = 0; /* counter for input pointers and scalars */
 n = floor(mxGetScalar(prhs[c++])+.5);
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

 hI = mxCalloc(n*N, sizeof(double*)); 
 for (i=0; i<n; i++)
   {
     for (ind=0; ind<N; ind++)
      {  
       f = mxGetCell(prhs[c], ind*n+i); 
       hI[ind*n+i] = mxGetPr(f);    /* get pointers to filtered images */       
      }
    }
 c++;

 varMap = mxGetPr(prhs[c++]);
 
 Upperbound = mxGetScalar(prhs[c++]);             /* satuation level */
 
 sx = floor(mxGetScalar(prhs[c++])+.5);
 sy = floor(mxGetScalar(prhs[c++])+.5);
 ax = floor(mxGetScalar(prhs[c++])+.5);
 ay = floor(mxGetScalar(prhs[c++])+.5);
 h = floor(mxGetScalar(prhs[c++])+.5);
 sub = floor(mxGetScalar(prhs[c++])+.5);
 
 for (i=0; i<n; i++)
 {
	calcGraHistMap(i);
 }
 calcVarMap();
}

