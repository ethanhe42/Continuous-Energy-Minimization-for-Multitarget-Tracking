#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    /* /Declarations */
    const mxArray *Dw2x, *Dw2y, *hh, *binsdata;
    double *xValues, *yValues;
    const char *hhValues;
    int i, j;
    int rowLen, colLen;

    int nbins;
    int ind;
    
    double *resx, *resy;
    
    if (!(mxIsUint8(prhs[2]))) {
        mexErrMsgTxt("Input argument 3 must be of type uint8.");
    }
    
    /* //Copy input pointer x */
    Dw2x = prhs[0];
    Dw2y = prhs[1];
    hh = prhs[2];
    binsdata = prhs[3];
    
    /* //Get matrix x */
    xValues = mxGetPr(Dw2x);
    rowLen = mxGetN(Dw2x);
    colLen = mxGetM(Dw2x);
    
    yValues = mxGetPr(Dw2y);
    hhValues = mxGetPr(hh);
    
    /* //Get the Integer */
    nbins = (int)(mxGetScalar(binsdata));
   
    plhs[0]=mxCreateDoubleMatrix(1, nbins, mxREAL);
    plhs[1]=mxCreateDoubleMatrix(1, nbins, mxREAL);
    resx = mxGetPr(plhs[0]);
    resy = mxGetPr(plhs[1]);
    
    /* //Print the integer avg of each col to matlab console */
        for(i=0;i<rowLen;i++) {
            for(j=0;j<colLen;j++) {
                ind=(j*rowLen)+i;
                resx[hhValues[ind]-1]-=xValues[ind];
                resy[hhValues[ind]-1]-=yValues[ind];

                        /*             res[ind] = xValues[ind]*yValues[ind]; */
            }
            
            
        }
    
    
    
    return;
    
}