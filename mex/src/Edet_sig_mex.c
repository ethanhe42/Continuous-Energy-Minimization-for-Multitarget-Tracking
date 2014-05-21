/* (C) Anton Andriyenko, 2012 */
#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {


	/* /Declarations */
	const mxArray *Xmtx, *Ymtx, *Xdmtx, *Ydmtx, *Sdmtx, *gridStepdata, *lambdadata, *tEdata, *xldata, *poccdata;

	double *X, *Y, *Xd, *Yd, *Sd, *occ, *tvis, *occx, *occy, *ddvix, *ddviy, *targetsExist, gridStep, lambda, pocc;


	int i,j,t;
	int F,N,Ndet,Fdet,FN;

	double *f;
	double *df;

	double d, ddfx, ddfy;

	double csig;
	int ind, indj, dfxind, xl, icnt, jcnt;
	int dims[3];
	
	double xdet, ydet, sdet, xxd, yyd, xxdsq, yydsq;
	double xxyycsig, xxyycsigsq;
	double xpos, ypos;
	double gg;
	double siga, sigb, epsilon, tmp1, tmp2;

	/* //Copy input pointer x */
	Xmtx = prhs[0];
	Ymtx = prhs[1];
	Xdmtx = prhs[2];
	Ydmtx = prhs[3];
	Sdmtx = prhs[4];
	gridStepdata = prhs[5];
	lambdadata = prhs[6];
	tEdata = prhs[7];
	xldata = prhs[8];
	poccdata = prhs[9];

	/* //Get matrix x */
	X = mxGetPr(Xmtx);
	Y = mxGetPr(Ymtx);
	Xd = mxGetPr(Xdmtx);
	Yd = mxGetPr(Ydmtx);
	Sd = mxGetPr(Sdmtx);
	targetsExist = mxGetPr(tEdata);
	N = mxGetN(Xmtx);
	F = mxGetM(Xmtx);
	FN=F*N;
	Ndet = mxGetN(Xdmtx);
	Fdet = mxGetM(Xdmtx);

	gridStep = (double)(mxGetScalar(gridStepdata));
	lambda = (double)(mxGetScalar(lambdadata));
	xl = (int)(mxGetScalar(xldata));
	pocc = (int)(mxGetScalar(poccdata));



	dims[0]=F;dims[1]=N;dims[2]=N;

	/* // Allocate memory and assign output pointer */
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(xl, 1, mxREAL);
	
	/*
// // 		plhs[2] = mxCreateDoubleMatrix(F, N, mxREAL);
// 		plhs[3] = mxCreateDoubleMatrix(F, N, mxREAL);
// 		plhs[4] = mxCreateDoubleMatrix(F, N, mxREAL);
// 		plhs[5] = mxCreateNumericArray(3, dims,mxDOUBLE_CLASS,   mxREAL);
// 		plhs[6] = mxCreateNumericArray(3, dims,mxDOUBLE_CLASS,   mxREAL);
// */


	/*
	//Get a pointer to the data space in our newly allocated memory
	*/
	f = mxGetPr(plhs[0]);
	df = mxGetPr(plhs[1]);
	
/* 	occ=mxGetPr(plhs[2]);
	occx=mxGetPr(plhs[3]);
	occy=mxGetPr(plhs[4]);
	ddvix=mxGetPr(plhs[5]);
	ddviy=mxGetPr(plhs[6]);  */
	*f=0;
	for (i=0; i<xl; i++) {
		df[i]=0;
	}




	/* if occlusions are on, call computeOcclusions2 */
	if (pocc!=0) {
		mxArray *out[6];
		mxArray *in[2];
		in[0] = mxDuplicateArray(prhs[0]);
		in[1] = mxDuplicateArray(prhs[1]);
		mexCallMATLAB(6,out,2,in,"computeOcclusions2");

		occ=mxGetPr(out[0]);
		tvis=mxGetPr(out[1]);
		occx=mxGetPr(out[2]);
		occy=mxGetPr(out[3]);
		ddvix=mxGetPr(out[4]);
		ddviy=mxGetPr(out[5]);


		/* copy occlusion computations */
		for (i=0; i<=5; i++) {
			plhs[i+2]=mxDuplicateArray(out[i]);
		}

	}
	else {
		double *p1, *p2, *p3, *p4, *p5, *p6;
		
		/* make sure to pass it to mex output */
		plhs[2] = mxCreateDoubleMatrix(F, N, mxREAL);
		plhs[3] = mxCreateDoubleMatrix(F, N, mxREAL);
		plhs[4] = mxCreateDoubleMatrix(F, N, mxREAL);
		plhs[5] = mxCreateDoubleMatrix(F, N, mxREAL);
		plhs[6] = mxCreateNumericArray(3, dims,mxDOUBLE_CLASS,   mxREAL);
		plhs[7] = mxCreateNumericArray(3, dims,mxDOUBLE_CLASS,   mxREAL);

/*
// 		p1 = mxGetPr(plhs[2]); memcpy(p1,occ, F*N * sizeof (double));
// 		p2 = mxGetPr(plhs[3]); memcpy(p2,tvis, F*N * sizeof (double));
// 		p3 = mxGetPr(plhs[4]); memcpy(p3,occx, F*N * sizeof (double));
// 		p4 = mxGetPr(plhs[5]); memcpy(p4,occy, F*N * sizeof (double));
// 		p5 = mxGetPr(plhs[6]); memcpy(p5,ddvix, F*N*N * sizeof (double));
// 		p6 = mxGetPr(plhs[7]); memcpy(p6,ddviy, F*N*N * sizeof (double));
	*/	
		p1 = mxGetPr(plhs[2]);
		p2 = mxGetPr(plhs[3]);
		p3 = mxGetPr(plhs[4]);
		p4 = mxGetPr(plhs[5]);
		p5 = mxGetPr(plhs[6]);
		p6 = mxGetPr(plhs[7]);
		
		for (i=0; i<F*N; i++) {
			p1[i]=1; p2[i]=1;
		}
		occ=p1;
		occx=p3;
		occy=p4;
		

	}


	/* now compute Edet	 */
	
	
	csig=gridStep*gridStep;
	dfxind=0;
	siga = 0.1;
	sigb = gridStep/2 * siga;
	epsilon = 1;
	
	/* in each target */
	for(i=0;i<N;i++)
	{

		int from = (int) targetsExist[i];
		int to = (int) targetsExist[i+N];
		for (t = from-1; t < to; t++)
		{
			gg=0;
			ind=(i*F)+t;
			xpos = X[ind];
			ypos = Y[ind];

			d=lambda*occ[ind];
			ddfx=lambda*occx[ind];
			ddfy=lambda*occy[ind];
			for (j = 0; j < Ndet; j++)
			{
				indj=(j*Fdet)+t;
				xdet = Xd[indj];

				if (xdet != 0)
				{
					ydet = Yd[indj];
					xxd = xdet - xpos;
					yyd = ydet - ypos;
					xxdsq = xxd*xxd;
					yydsq = yyd*yyd;
					sdet = Sd[indj];
					/* sdet=1;					*/
					


					 // gg -= sdet*csig / (xxdsq + yydsq + csig);
					gg += sdet * (1. / (exp(-siga * sqrt(xxdsq + yydsq + epsilon) + sigb)+1) - 1);
					// gg += xpos;
					// printf("%d %d %d %f\n",t,i,j,xpos);
					// printf("%f %f %f %f %f %f %f %f\n",sdet,siga,sigb,xxdsq,yydsq,epsilon,sqrt(xxdsq + yydsq + epsilon),gg);
					// df[dfxind]   += 1;
					if (nlhs>1) { 
						xxyycsig  = xxdsq + yydsq + csig;
						xxyycsigsq = xxyycsig * xxyycsig;
						
						 // df[dfxind] -= 2*csig*sdet * xxd / xxyycsigsq;						 
						 // df[dfxind+1] -= 2*csig*sdet * yyd / xxyycsigsq;
						 
						tmp2= exp(sigb - siga*sqrt(epsilon + xxdsq + yydsq));
						tmp1 = (tmp2 + 1);
						tmp1 = tmp1*tmp1;
						
						df[dfxind]   += -(sdet*siga*tmp2*(2*xdet - 2*xpos))/(2*tmp1*sqrt(epsilon + xxdsq + yydsq));
						df[dfxind+1] += -(sdet*siga*tmp2*(2*ydet - 2*ypos))/(2*tmp1*sqrt(epsilon + xxdsq + yydsq));


						// df[dfxind]   += -(sdet*siga*tmp2*(2*xdet - 2*xpos))/(2*tmp1*sqrt(epsilon + xxdsq + yydsq));
						// df[dfxind+1] += -(sdet*siga*tmp2*(2*ydet - 2*ypos))/(2*tmp1*sqrt(epsilon + xxdsq + yydsq));
						
						// df[dfxind]   += 1;
						// df[dfxind+1] += 1;


					}


				}
			}
			// printf("%f %f\n",gg,df[dfxind]);

			*f+=(d+gg);
			// *f += gg;

			df[dfxind]+=ddfx;
			df[dfxind+1]+=ddfy;
			dfxind+=2;

		}
	}
	


}