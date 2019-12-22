/**************************************************
 Active Basis Learning
***************************************************/

/*
Note: modified to include normalized data weights.
*/

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
 double **mI;                /* locally maximized images */
 double **pI;                /* pooled maximized image */
 double **C;                 /* inhibition coefficients */
 int h;                      /* halfsizes of filters */
 int sx, sy;                 /* sizes of image */ 
 double T;                   /* threshold of selected bases */
 double *Mi, *Mx, *My, *Mm, *Mm1;  /* storing selected bases */
 double **allsymbol;         /* symbol of filters */
 double *sym, **Asym;        /* symbol of selected Gabors */
 int L, ore;                 /* allowed ranges of shifting in location and orientation */
 int sub;                    /* subsampling */
 int Totalsketch;            /* total number of Gabors */
 double Upperbound;          /* If response is greater than Upperbound, then set response = Upperbound */
 double *tss;                /* score for each example */
 double *dataWeight;		 /* normalized data weights of the data instances */
 double *gain;               /* coding gain for each selected base per example */
 double SHUTUP;              /* thresholding for shutting down a candidate basis element */
 int **sinsh, **cossh;       /* store the shift to avoid repeated sin and cosin computation */
 int **xI, **yI, **indI;      /* keep track the shifted base for local maximum pooling */ 
 int binnum;                  /* number of bins to keep the pooled responses */
 double binsize;              /* length of each bin */
 int *headx, *heady, *headind; /* initial position and orientation in each bin */
 int **nxI, **nyI, **nindI;    /* next position and location */
 int *tailx, *taily, *tailind; /* last position and orientation in each bin */
 int **bxI, **byI, **bindI;    /* precedent position and orientation */
 int **binI;                   /* the bin that a certain position orientation belongs to */
 int topbin;                   /* the top non-empty bin */
 int M;   
 double *lam, *e, *lz;
 
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
 
   m = 0.; 
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


/* the Gabor(mx, my, orientation = mi) inhibits overlapping Gabors for image i */
void inhibit(int i, int mi, int mx, int my) 
{
   int x0, y0, x, y, ind, here, rm[5], x1, y1, ii, ind1; 
   double *f, *fc, mm, mm1;
   
   /* set neighboring basis elements to be 0 */
   for (ind=0; ind<N; ind++)   
     {
       f = fI[ind*n+i];   
       fc = C[mi+ind*N];   /* inhibition coefficients, zero = inhibit */
       for (x=MAX(1, mx-2*h); x<=MIN(sx, mx+2*h); x++)
         for (y=MAX(1, my-2*h); y<=MIN(sy, my+2*h); y++)
         {
          f[px(x, y, sx, sy)] *= fc[px(x-mx+2*h+1, y-my+2*h+1, 4*h+1, 4*h+1)];
         }
      }
   
   /* adjust the local maximum caused by inhibition */
   for (ind=0; ind<N; ind++)   
     {
      for (x0 = MAX( floor((h+1)/sub)+1+L, floor((mx-h*2-L*sub)/sub)+1); x0 <= MIN (floor((sx-h)/sub)-1-L, floor((mx+2*h+L*sub)/sub)) ; x0++) /* modified !! */
         for (y0=MAX(floor((h+1)/sub)+1+L, floor((my-2*h-L*sub)/sub)+1); y0<=MIN(floor((sy-h)/sub)-1-L, floor((my+2*h+L*sub)/sub)); y0++) /* modified */
         {
           x = x0*sub; y = y0*sub; 
           ii = ind*n+i;
           here = px(x, y, sx, sy); 
           if (pI[ind][here]>0.)   /* if the orientation/location is still in the game */
           {
             mm = mI[ii][here];
             x1 = xI[ii][here];
             y1 = yI[ii][here];
             ind1 = indI[ii][here];   /* find the local maximum position and orientation */  
              /* if the local maximum is within the range of inhibition */
             if ((x1-mx>=-2*h)&&(x1-mx<=2*h)&&(y1-my>=-2*h)&&(y1-my<=2*h)) 
             {
                 /* if the local maximum is indeed inhibited */
                 if(C[mi+ind1*N][px(x1-mx+2*h+1,y1-my+2*h+1,4*h+1,4*h+1)]==0.)
                 {
                     mm1 = shiftmax(i, ind, x, y, rm);  /* re-locate the local maximum */
                     xI[ii][here] = rm[1];
                     yI[ii][here] = rm[2];
                     indI[ii][here] = rm[3];    /* update re-located local maximum */
                     mI[ii][here] = mm1;
                     /* update the pooled sum */
                     pI[ind][here] += (mm1-mm)*dataWeight[i];
                 }
             }
           }
         }
   }
}


/* draw the symbol of Gabor(x0, y0, orientation = ind) with intensity w */
void draw(double *sym, int x0, int y0, int ind, double w)
{
  int x, y; 
  double a; 
          
  for (x=x0-h; x<=x0+h; x++)
     for (y=y0-h; y<=y0+h; y++)
       {
         a = allsymbol[ind][px(x-x0+h+1, y-y0+h+1, 2*h+1, 2*h+1)]*w; 
         if (sym[px(x, y, sx, sy)]<a)
             sym[px(x, y, sx, sy)] = a;
       }
}

/* local maximization */
void Cmm()
{
   int i, ind, x, y, rm[5], here, ii;  
   
   /* we keep track the position and orientation of local maximum explicitly */
   pI = double_matrix(N, sx*sy);
   mI = double_matrix(n*N, sx*sy);
   xI = int_matrix(n*N, sx*sy);
   yI = int_matrix(n*N, sx*sy);
   indI = int_matrix(n*N, sx*sy);   /* keep track the position and orientation of the local maximum */
      
   for (x=floor((h+1)/sub)+1+L; x<=floor((sx-h)/sub)-1-L; x++)
      for (y=floor((h+1)/sub)+1+L; y<=floor((sy-h)/sub)-1-L; y++)
       {
        here = px(x*sub, y*sub, sx, sy); 
        for (ind=0; ind<N; ind++)
           {
             pI[ind][here] = 0.;
             for (i=0; i<n; i++)
             {
                ii = ind*n+i;
                mI[ii][here] = shiftmax(i, ind, x*sub, y*sub, rm);   /* local maximum of active Gabor */  
                
                pI[ind][here] += mI[ii][here]*dataWeight[i];
                
                xI[ii][here] = rm[1];      /* record the position of local maximum */
                yI[ii][here] = rm[2]; 
                indI[ii][here] = rm[3];    /* record the orientation of local maximum */
              }
             if (pI[ind][here]<SHUTUP) 
                 pI[ind][here] = 0.;
            }
       }
}

void putfront(int b, int x, int y, int ind)
{
    int hx0, hy0, hi0, here0, here;
    
    here = px(x, y, sx, sy);
    if (headx[b] == sx+1) 
	  {
		headx[b] = x; heady[b] = y; headind[b] = ind;
		nxI[ind][here] = sx+1; nyI[ind][here] = sy+1; nindI[ind][here] = N+1;
		bxI[ind][here] = -1; byI[ind][here] = -1;  bindI[ind][here] = -1;
	  }
    else
      {
        hx0 = headx[b]; hy0 = heady[b]; hi0 = headind[b];
        here0 = px(hx0, hy0, sx, sy);
        headx[b] = x; heady[b] = y; headind[b] = ind;
        nxI[ind][here] = hx0; nyI[ind][here] = hy0; nindI[ind][here] = hi0;
        bxI[ind][here] = -1; byI[ind][here] = -1;  bindI[ind][here] = -1;
        bxI[hi0][here0] = x; byI[hi0][here0] = y; bindI[hi0][here0] = ind;
      }
}

void Cbinchain()
{
  int b, ind, x, y, here, x0, y0;
   /* the following divides the pooled responses into bins, and chain them in each bin
      the purpose is to avoid searching all the positions and orientations 
      we only need to go through the top non-empty bin in each iteration of matching pursuit */
   headx = int_vector(binnum);
   heady = int_vector(binnum);
   headind = int_vector(binnum);   /* initial position and orientation of each bin */
   nxI = int_matrix(N, sx*sy);   
   nyI = int_matrix(N, sx*sy); 
   nindI = int_matrix(N, sx*sy);   /* next position and orientation in the chain for each bin */
   
   bxI = int_matrix(N, sx*sy);    
   byI = int_matrix(N, sx*sy); 
   bindI = int_matrix(N, sx*sy);   /* precedent position and orientation in the chain for each bin */
   
   binI = int_matrix(N, sx*sy); 
   
   for (b=0; b<binnum; b++)
   {
       headx[b] = sx+1; 
       heady[b] = sy+1; 
       headind[b] = N+1; 
   }
      
   for (x0=floor((h+1)/sub)+1+L; x0<=floor((sx-h)/sub)-1-L; x0++)
      for (y0=floor((h+1)/sub)+1+L; y0<=floor((sy-h)/sub)-1-L; y0++)
       {
        x = x0*sub; 
        y = y0*sub; 
        here = px(x, y, sx, sy); 
        for (ind=0; ind<N; ind++)
           {
             if(pI[ind][here] > 0.)
             {
              b = MIN(floor((pI[ind][here]-SHUTUP)/binsize), binnum-1);   /* find the bin number */
              binI[ind][here] = b;
              /* add (x, y, ind) to the front of the chain in bin b */
              putfront(b, x, y, ind);
             }
             else
             {
                 binI[ind][here] = -1; 
             }
        }
      }
   
   topbin = binnum - 1;
   while (headx[topbin] == sx+1)
   {
          topbin --;
   }
}

/* the shared sketch algorithm */
void Cmp()
{
   int i, ind, x0, y0, x, y, mi, mx, my, t, here, xi, yi, indi, ii, b, b0, j; 
   int nx, ny, nind, nhere, bx, by, bind, bhere;
   double m, mm, ehat, ov;
   int count;
   t = 0; 
   do
   {
     m = 0.;
     x = headx[topbin]; y = heady[topbin]; ind = headind[topbin]; 
     count = 0;
     do
     {
     here = px(x, y, sx, sy);
     if (m<pI[ind][here])
     {
         m = pI[ind][here]; 
         mi = ind; mx = x; my = y;
     }
     x = nxI[ind][here]; y = nyI[ind][here]; ind = nindI[ind][here]; 
     count++;
     }
     while (x!=sx+1 && count < 10000);
     
     Mi[t] = mi; Mx[t] = mx; My[t] = my;
     ehat = m;
     j = M-1;
     while (e[j]>ehat)
		j--;
 
     if (j==M-1)
      {
        Mm[t] = lam[j];
        Mm1[t] = lz[j];
       }
     else
     {
        ov = (ehat-e[j])/(e[j+1]-e[j]);
        Mm[t] = lam[j]+(lam[j+1]-lam[j])*ov;
        Mm1[t] = lz[j]+(lz[j+1]-lz[j])*ov;
     }
     
      draw(sym, mx, my, mi, sqrt(ehat));            
        
     here = px(mx, my, sx, sy);
     for (i=0; i<n; i++)
        {
          ii = mi*n+i; 
          mm = mI[ii][here];
          tss[i] += Mm[t]*mm - Mm1[t];
          xi = xI[ii][here];
          yi = yI[ii][here];
          indi = indI[ii][here];
           
          if (mm>0.)
          {  
               draw(Asym[i], xi, yi, indi, sqrt(mm));
               inhibit(i, indi, xi, yi);
          }
        }

     
     /* go through the affected pixels and re-bin and re-chain them */     
      for (x0=floor(MAX(1, mx-h*2-L*sub*2)/sub)+1; x0<=floor(MIN(sx, mx+2*h+L*sub*2)/sub); x0++)
         for (y0=floor(MAX(1, my-2*h-L*sub*2)/sub)+1; y0<=floor(MIN(sy, my+2*h+L*sub*2)/sub); y0++)
         {
           x = x0*sub; y = y0*sub; 
           here = px(x, y, sx, sy); 
           for (ind=0; ind<N; ind++)   
           {
           if (pI[ind][here]<SHUTUP)
             {
               pI[ind][here] = 0.;
               b = -1;
             }
           else 
              {
                 b = MIN(floor((pI[ind][here]-SHUTUP)/binsize), binnum-1);
               }
           b0 = binI[ind][here];
           if (b != b0)
                 {
                     binI[ind][here] = b; 
                     
                     /* remove (x, y, ind) from the chain in bin b0 */
                     nx = nxI[ind][here]; ny = nyI[ind][here]; nind = nindI[ind][here];         
                     bx = bxI[ind][here]; by = byI[ind][here]; bind = bindI[ind][here];         
                     if(bx != -1)
                     {
                         if (bx != 0) /* IF branch added */
						 { 
						  bhere = px(bx, by, sx, sy); 
						  nxI[bind][bhere] = nx; nyI[bind][bhere] = ny; nindI[bind][bhere] = nind; 
						 }
                     }
                     else 
                     {
                         headx[b0] = nx; heady[b0] = ny; headind[b0] = nind; 
                     }
                     
                     if(nx != sx+1)
                     {
                         if(nx != 0) /* IF branch added */
						 {
						  nhere = px(nx, ny, sx, sy); 
						  bxI[nind][nhere] = bx; byI[nind][nhere] = by; bindI[nind][nhere] = bind; 
						 }
                     }
                   
                  if (b != -1)
                     {
                     /* add (x, y, ind) to the front of the chain in bin b */
						putfront(b, x, y, ind); 
                     }
                }
           }
        }
        
         while (headx[topbin] == sx+1)
         {
          topbin --;
		  /* added: */
		  if(topbin==0)
		  {
			  printf("notice!!\n");
		  }
          }
     gain[t] = Mm[t]*ehat-Mm1[t];
      t++;
   }
  while (t<Totalsketch);
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
 dataWeight = mxGetPr(prhs[c++]);      /* normalized data weights of data items */
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
 T = mxGetScalar(prhs[c++]);              /* threshold */  
     
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
 
 sym = mxGetPr(prhs[c++]);   /* symbols of selected Gabors for all the images */          
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

 Totalsketch = floor(mxGetScalar(prhs[c++])+.5);  /* total number of bases to be selected */
 Upperbound = mxGetScalar(prhs[c++]);       /* saturation level */
 tss = mxGetPr(prhs[c++]);                  /* log-likelihood scores for examples */
 gain = mxGetPr(prhs[c++]);                 /* coding gain for selected base */
 
 SHUTUP = mxGetScalar(prhs[c++]); 
 binnum = floor(mxGetScalar(prhs[c++])+.5); 
 
 M = floor(mxGetScalar(prhs[c++])+.5);  
 lam = mxGetPr(prhs[c++]);   
 e = mxGetPr(prhs[c++]);   
 lz = mxGetPr(prhs[c++]);   

 binsize = (Upperbound-SHUTUP)/binnum; 
 
 storeshift();
 Cmm(); 
 Cbinchain();
 Cmp();
}
