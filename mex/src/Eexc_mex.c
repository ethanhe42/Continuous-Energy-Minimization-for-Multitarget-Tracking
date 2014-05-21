/* (C) Anton Andriyenko, 2012 */

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {


    /* /Declarations */
    const mxArray *Xmtx, *Ymtx, *gridStepdata, *itToInddata, *xldata;

    double *X, *Y, gridStep;
	double *itToInd;

    int i,j,t, xl;
    int F,N;
	int F1,N1;

    double *f;
    double *df;

	double gs_squared;
	int ind, indj, xind, yind;

    /* //Copy input pointer x */
    Xmtx = prhs[0];
	Ymtx = prhs[1];
    gridStepdata = prhs[2];
    itToInddata = prhs[3];
	xldata = prhs[4];

    /* //Get matrix x */
    X = mxGetPr(Xmtx);
	Y = mxGetPr(Ymtx);
    N = mxGetN(Xmtx);
    F = mxGetM(Xmtx);
    N1 = mxGetN(itToInddata);
    F1 = mxGetM(itToInddata);

    gridStep = (double)(mxGetScalar(gridStepdata));
    itToInd = mxGetPr(itToInddata);
	xl = (int)(mxGetScalar(xldata));

    /* Allocate memory and assign output pointer */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(xl, 1, mxREAL);

    /* Get a pointer to the data space in our newly allocated memory */
    f = mxGetPr(plhs[0]);
    df = mxGetPr(plhs[1]);


	gs_squared=gridStep*gridStep;
	xind=0;
	yind=0;



/* for each frame */
for(t=0;t<F;t++)
{
	/* for each existing target */
	for (i = 0; i < N; i++)
    {
		double a;
	  ind=(i*F)+t;
      a = X[ind];
      if (a != 0)
      {
        for (j = 0; j < N; j++)
        {
		  double c;
		  indj=(j*F)+t;
          c = X[indj];
          if (i != j && c != 0)
          {
			const double b = Y[ind];
            const double d = Y[indj];
            double tmp, dda, ddb;

            double denom = a*a - 2*a*c + b*b - 2*b*d + c*c + d*d;
            *f += gs_squared / denom;

			if (nlhs>1) {
				
				denom *= denom;
				tmp = gs_squared / denom;
				dda = -(2*a-2*c)*tmp;
				ddb = -(2*b-2*d)*tmp;
				
				xind=(int)itToInd[ind]-1;
				df[xind] += dda + dda;
				df[xind+1] += ddb + ddb;
				
			}


          }
        }
      }
    }
}
}