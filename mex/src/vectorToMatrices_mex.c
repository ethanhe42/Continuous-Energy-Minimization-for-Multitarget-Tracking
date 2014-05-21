/* (C) Anton Andriyenko, 2012 */

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* Declarations */ 
    const mxArray *xdata, *tiToIndData, *Fdata, *Ndata;
    double *x, *tiToInd;
    
    double *X, *Y;
    int i,t,ind, F,N, tiToIndInt;
    
    
	/* //Copy input pointer x */
	xdata = prhs[0];
	tiToIndData = prhs[1];
    Fdata = prhs[2];
    Ndata = prhs[3];
    
    /* //Get matrix x */
    x = mxGetPr(xdata);
    tiToInd = mxGetPr(tiToIndData);
    
    F = (int)(mxGetScalar(Fdata));
    N = (int)(mxGetScalar(Ndata));
    
	/* // Allocate memory and assign output pointer */
	plhs[0] = mxCreateDoubleMatrix(F, N, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(F, N, mxREAL);
    
	/*
	//Get a pointer to the data space in our newly allocated memory
	*/
	X = mxGetPr(plhs[0]);
	Y = mxGetPr(plhs[1]);
    
    for (i=0; i<N; i++) {
        for (t=0; t<F; t++)            
            /*             mexPrintf("%d %d %d %d %d\n",i,t,F,N,(int)tiToInd[i*F+t]);*/
        {
            ind=i*F+t;
            tiToIndInt=(int)tiToInd[ind];
            if (tiToIndInt) {
                X[ind]=x[tiToIndInt-1];
                Y[ind]=x[tiToIndInt];
            }
        }
    } 
            
}
