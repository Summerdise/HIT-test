/**************************************************
Shared Sketch by Active Bases --- Testing Phase
***************************************************/
# include <stdio.h>
# include <stdlib.h>
# include "mex.h"        /* the algorithm is connect to matlab */
# include "math.h"
# define PI 3.1415926
# define ABS(x) ((x)>0? (x):(-(x)))
# define MAX(x, y) ((x)>(y)? (x):(y))
# define MIN(x, y) ((x)<(y)? (x):(y))
# define NEGMAX -1e10

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

void free_matrix(void **a, int nrow, int ncol)
{
	int i;
	for(i = 0; i < nrow; i++)
		mxFree(a[i]);
	mxFree(a);
}



/* getting the index of matlab image, (x,y) location, (sx, sy) sizes */
int px(int x, int y, int bx, int by)    
{            
	return (x + (y-1)*bx - 1); 
}



int n;                      /* number of images */
int N;                      /* number of orientations */
double **fI;                /* filtered images */
double **C;                 /* inhibition coefficients */
int h;                      /* halfsizes of filters */
int sx, sy;                 /* sizes of image */ 
double *Mi, *Mx, *My, *Mm, *Mm1;  /* storing selected bases */
double **allsymbol;         /* symbol of filters */
double *sym, **Asym;        /* symbol of selected Gabors */
int L, ore;                 /* allowed ranges of shifting in location and orientation */
int sub;                    /* subsampling */
double *tss;                /* log-likelihood score for each testing image */
double Upperbound;          /* satuation level for response */
int Totalsketch;            /* total number of bases */
int **sinsh, **cossh;       /* store the shift to avoid repeated sin and cosin computation */
int M;                      
double *lam, *e, *lz;  

double *bkp, *lbkp;				/* background pooling*/


/* store the shift values so that we do not need to repeat the sin and cos computation */ 
void storeshift()
{
	int ind, l;
	double theta; 

	sinsh = int_matrix(N, L+L+1); 
	cossh = int_matrix(N, L+L+1); 
	for (ind=0; ind<N; ind++)        
	{
		theta = PI*ind/N; 
		for (l=-L; l<=L; l++)
		{
			sinsh[ind][l+L] = floor(l*sub*sin(theta)+.5); 
			cossh[ind][l+L] = floor(l*sub*cos(theta)+.5); 
		}   
	}
}

/* for Gabor(x, y, orientation = ind), find local maximum in image i */
double shiftmax(int i, int ind, int x, int y, int *rm) 
{
	double m;
	int x1, y1, l, ml, mx, my, de, d, d1, here, mo, md; 

	m = NEGMAX; 
	for (l=-L; l<=L; l++)   
	{
		x1 = x + cossh[ind][l+L]; 
		y1 = y + sinsh[ind][l+L];   /* shifting in normal direction */
		if ((x1>=1)&&(x1<=sx)&&(y1>=1)&&(y1<=sy))
		{
			here = px(x1, y1, sx, sy);
			for (de=-ore; de<=ore; de++)
			{
				d = de+ind;
				d1 = d;
				if (d<0)
					d1 = d+N;
				if (d>=N)
					d1 = d-N;   /* shifting in orientation */
				if (fI[d1*n+i][here]>m)
				{
					m = fI[d1*n+i][here];   /* local maximization */
					ml = l; mx = x1; my = y1; mo = d1; md = de; 
				}
			}
		}
	}
	rm[0] = ml; rm[1] = mx; rm[2] = my; rm[3] = mo; rm[4] = md; 
	return(MIN(Upperbound, m)); 
}



/* draw the symbol of Gabor(x0, y0, orientation = ind) with intensity w */
void draw(double *sym, int x0, int y0, int ind, double w)
{
	int x, y; 
	double a; 

	for (x=x0-h; x<=x0+h; x++)
	{
		for (y=y0-h; y<=y0+h; y++)
		{
			a = allsymbol[ind][px(x-x0+h+1, y-y0+h+1, 2*h+1, 2*h+1)]*w; 
			if (sym[px(x, y, sx, sy)]<a)
				sym[px(x, y, sx, sy)] = a;
		}
	}
}


/* the shared sketch algorithm */
void Cmp()
{
	int i, mi, mx, my, t, rm[5]; 
	double mm; 
	int j; j=0; 

	t = 0;
	do
	{
		mi = Mi[t]; mx = Mx[t]; my = My[t];   
		for (i=0; i<n; i++)
		{            
			mm = shiftmax(i, mi, mx, my, rm);  
			tss[i+t*n] += Mm[t]*mm - Mm1[t]; 
			/* tss[i+t*n] -= bkp[n*mi+i]*mm - lbkp[n*mi+i]; */
			if (mm>0.)
			{  
				draw(Asym[i], rm[1], rm[2], rm[3], sqrt(mm)); 
			}
		} 
		t++; 
	}
	while (t<Totalsketch);   /* stopping criterion */
}

void background_pooling()
{
	int ii, jj, hh;

	int j;

	int x, y;

	double ehat, ov;

	for(ii=0; ii<n; ++ii)
	{
		for(jj=0; jj<N; ++jj)
		{
			ehat = 0;

			for(x=1; x<=sx; ++x)
			{
				for(y=1; y<=sy; ++y)
				{
					ehat += fI[jj*n+ii][px(x, y, sx, sy)];
					/*shiftmax(ii, jj, hh%sx+1, hh/sx+1, rm);*/
				}
			}
			ehat /= (sx*sy);
						
			
			j = M-1; 
			while (e[j]> ehat)
			{
				j--;
				if(0==j)	break;
			}

			if (j==M-1)
			{
				bkp[jj*n+ii] = lam[j]; 
				lbkp[jj*n+ii] = lz[j]; 
			}
			else 
			{
				ov = (ehat-e[j])/(e[j+1]-e[j]); 
				bkp[jj*n+ii] = lam[j]+(lam[j+1]-lam[j])*ov; 
				lbkp[jj*n+ii] = lz[j]+(lz[j+1]-lz[j])*ov; 

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
	int ind, i, j, c; 
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
	
	


	C = mxCalloc(N*N, sizeof(double*));    /* C: correlation/inhibition between filters */
	for (ind=0; ind<N; ind++)
	{  
		for (j=0; j<N; j++)
		{
			f = mxGetCell(prhs[c], j*N+ind); 
			C[j*N+ind] = mxGetPr(f);         /* get correlation/inhibition */
		}   
	}
	c++; 

	h = floor(mxGetScalar(prhs[c++])+.5);    /* half size of filters */

	sx = floor(mxGetScalar(prhs[c++])+.5);
	sy = floor(mxGetScalar(prhs[c++])+.5);   /* size of images */

	Mi = mxGetPr(prhs[c++]);   /* orientation of selected Gabor */               
	Mx = mxGetPr(prhs[c++]);         
	My = mxGetPr(prhs[c++]);   /* (x, y) position of selected Gabor */
	Mm = mxGetPr(prhs[c++]);   /* lambda */    
	Mm1 = mxGetPr(prhs[c++]);  /* logZ */ 
	
	

	allsymbol = mxCalloc(N, sizeof(double*));    
	for (ind=0; ind<N; ind++)
	{  
		f = mxGetCell(prhs[c], ind); 
		allsymbol[ind] = mxGetPr(f);  /* symbols of filters */         
	}
	c++; 

	Asym = mxCalloc(n, sizeof(double*));   /* symbols of active Gabors for an individual image */
	for (i=0; i<n; i++)
	{
		f = mxGetCell(prhs[c], i);
		Asym[i] = mxGetPr(f);  
	}
	c++; 

	L = floor(mxGetScalar(prhs[c++])+.5);     /* range of shifting along normal direction of Gabor */
	ore = floor(mxGetScalar(prhs[c++])+.5);   /* range of shifting in orientation */
	sub = floor(mxGetScalar(prhs[c++])+.5);   /* sub-sampling of Gabor filters to improve speed */
	Totalsketch = floor(mxGetScalar(prhs[c++]));  /* total number of bases */
	Upperbound = mxGetScalar(prhs[c++]);             /* satuation level */
	tss = mxGetPr(prhs[c++]);                        /* map of log-likelihood score */
	M = floor(mxGetScalar(prhs[c++])+.5);  
	lam = mxGetPr(prhs[c++]);   
	e = mxGetPr(prhs[c++]);   
	lz = mxGetPr(prhs[c++]);   

	bkp = mxGetPr(prhs[c++]);
	/*bkp = mxCalloc(n*N, sizeof(double));*/   /* storage for background pooling */
	lbkp = mxCalloc(n*N, sizeof(double));
	
	storeshift();
	
	/* background_pooling(); */
	mexPrintf( "before cmp()\n" );
	mexEvalString( "drawnow" );
	Cmp();

	/*mxFree(bkp);*/
	/*mxFree(lbkp);*/
}
