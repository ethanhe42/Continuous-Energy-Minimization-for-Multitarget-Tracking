#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    /* /Declarations */
    mxArray *xData, *yData, *zData;
    double *xValues;
    char *yValues;
    int i, j;
    int rowLen, colLen;
    double avg;
    int nbins;
    int ind;

    double *res;

    if (!(mxIsUint8(prhs[1]))) {
        mexErrMsgTxt("Input argument must be of type uint8...");
    }

    /* //Copy input pointer x */
    xData = prhs[0];
    yData = prhs[1];
    zData = prhs[2];


    /* //Get matrix x */
    xValues = mxGetPr(xData);
    rowLen = mxGetN(xData);
    colLen = mxGetM(xData);

    yValues = mxGetPr(yData);

    /* //Get the Integer */
    nbins = (int)(mxGetScalar(zData));

    plhs[0]=mxCreateDoubleMatrix(1, nbins, mxREAL);
    res = mxGetPr(plhs[0]);

    /* //Print the integer avg of each col to matlab console */
        for(i=0;i<rowLen;i++) {
            for(j=0;j<colLen;j++) {
                ind=(j*rowLen)+i;
                res[yValues[ind]-1]+=xValues[ind];

                        /*             res[ind] = xValues[ind]*yValues[ind]; */
            }


        }



    return;

}