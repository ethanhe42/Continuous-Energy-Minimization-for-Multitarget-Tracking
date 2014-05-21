/* (C) Anton Andriyenko, 2012 */

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {


	/* Declarations */
	const mxArray *Xmtx, *Ymtx, *targetSizedata,  *tEdata, *xldata;

	double *X, *Y, targetSize, *targetsExist, xl;

	int i,t;
	int F,N;

	double *fx;
	double *df;

	int ind, indp, indn, xind, yind;

	/* //Copy input pointer x */
	Xmtx = prhs[0];
	Ymtx = prhs[1];
	targetSizedata = prhs[2];
	tEdata = prhs[3];
	xldata = prhs[4];

	/* //Get matrix x */
	X = mxGetPr(Xmtx);
	Y = mxGetPr(Ymtx);
	N = mxGetN(Xmtx);
	F = mxGetM(Xmtx);

	targetSize = (double)(mxGetScalar(targetSizedata));
	targetsExist = mxGetPr(tEdata);
	xl = (int)(mxGetScalar(xldata));


	/* // Allocate memory and assign output pointer */
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(xl, 1, mxREAL);

	/* //Get a pointer to the data space in our newly allocated memory  */
	fx = mxGetPr(plhs[0]);
	df = mxGetPr(plhs[1]);

	xind=0;	yind=0;

	/* for each target */
	for (i = 0; i < N; i++)
	{

	int from = (int) targetsExist[i];
	int to = (int) targetsExist[i+N];
	int tlength=to-from+1;
	xind += 2;
	yind = xind+1;

	for (t = from; t < to-1; t++)
	{
		double a,b,c,d,e,f;
		double diff;
		indp=(i*F)+t-1;
		ind=(i*F)+t;
		indn=(i*F)+t+1;
		a = X[indp];      b = Y[indp];
		c = X[ind];      d = Y[ind];
		e = X[indn];      f = Y[indn];

		diff = a*a + a*(2*e - 4*c) + b*b + b*(2*f - 4*d) + 4*c*c  - 4*c*e + 4*d*d  - 4*d*f + e*e  + f*f;
		*fx += diff;

		if (nlhs>1)
		{
			df[xind-2] += (2*a-4*c+2*e);
			df[yind-2] +=  (2*b-4*d+2*f);
			df[xind] += (-4*a+8*c-4*e);
			df[yind] += (-4*b+8*d-4*f);
			df[xind+2] += (2*a-4*c+2*e);
			df[yind+2] += (2*b-4*d+2*f);
		}

		xind += 2;
		yind = xind+1;
		}
		if (tlength>1) {
			xind += 2;
		}
	}

	/* normalize with target Size */
	*fx = *fx / targetSize;
	for (i = 0; i < xl; ++i)
	{
		df[i] = df[i] / targetSize;
	}
}