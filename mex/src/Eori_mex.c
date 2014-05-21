/* (C) Anton Andriyenko, 2012 */

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {


	/* Declarations */
	const mxArray *Xmtx, *Ymtx, *targetSizedata,  *tEdata, *xldata;
	const mxArray *Xdmtx, *Ydmtx, *Dxdmtx, *Dydmty;
	
	double *X, *Y, *Xd, *Yd, *Dx, *Dy, targetSize, *targetsExist, xl;

	int i,t, j,k,Ndet, Fdet;
	int indj, indk;
	int F,N;

	double *fx;
	double *df;
	double d1[2], d2[2], o1[2], o2[2];
	double xdet, ydet;
	double d1norm, d2norm, v1norm, v2norm;
	double so, epsilon, csig;
	double m1,m2,ff1,ff2,f1,f2,f1acos,f2acos, dist1, dist2, g1, g2;
	double f1b,f2b,ff1b,ff2b;
	double d1x,d1y,d2x,d2y;
	double det1x, det1y, det2x, det2y;
	double diff1, diff2, diff3, diff4, diff5, diff6;
	double amc, bmd, cme, dmf, amcs, bmds, cmes, dmfs;
	double cma, dmb, emc, fmd;
	double exp1,exp2,exp3,exp4,exp5,exp6,exp7,exp8;
	double exp9,exp10,exp11,exp12,exp13,exp14,exp15,exp16;
	double exp17,exp18,exp19,exp20,exp21,exp22,exp23,exp24;
	double exp25,exp26,exp27,exp28,exp29,exp30,exp31,exp32;
	double tmp1, tmp2, tmp3, tmp4;
	
	double a2c2, b2d2, c2e2, d2f2;
	double cutoffthr;

	int ind, indp, indn, xind, yind;

	double pauseTime = 1.0;
	double wt=0.8; 
	double kappa=.5;
/* 	double intv=13.3865972477803;*/
	double intv=13.3837563612694;
	double distrconst=0.318309886183791;
    mxArray *prhsp = mxCreateDoubleScalar(pauseTime);
    
	
	/* //Copy input pointer x */
	Xmtx = prhs[0];
	Ymtx = prhs[1];
	Xdmtx = prhs[2];
	Ydmtx = prhs[3];
	Dxdmtx = prhs[4];
	Dydmty = prhs[5];
	
	targetSizedata = prhs[6];
	tEdata = prhs[7];
	xldata = prhs[8];

	/* //Get matrix x */
	X = mxGetPr(Xmtx);
	Y = mxGetPr(Ymtx);
	Xd = mxGetPr(Xdmtx);
	Yd = mxGetPr(Ydmtx);
	Dx = mxGetPr(Dxdmtx);
	Dy = mxGetPr(Dydmty);
	N = mxGetN(Xmtx);
	F = mxGetM(Xmtx);
	Ndet = mxGetN(Xdmtx);
	Fdet = mxGetM(Xdmtx);
	
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

	so=10;
	epsilon=1;
	csig=targetSize;
	cutoffthr=50;

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

		amc=(a-c );            bmd=(b-d );
		cme=(c-e );            dmf=(d-f );
		cma=(c-a );            dmb=(d-b );
		emc=(e-c );            fmd=(f-d );
		amcs=amc*amc;            bmds=bmd*bmd;
		cmes=cme*cme;            dmfs=dmf*dmf;

		m1=sqrt(amcs+bmds+epsilon);
		m2=sqrt(cmes+dmfs+epsilon);

		/* auxiliary to speed up */
		
		o1[0]=cma/m1; o1[1]=dmb/m1;
		o2[0]=emc/m2; o2[1]=fmd/m2;
		
		if (nlhs>1) {

			a2c2 = (a*2.0-c*2.0);
			b2d2 = (b*2.0-d*2.0);
			c2e2 = (c*2.0-e*2.0);
			d2f2 = (d*2.0-f*2.0);
			exp1 = epsilon+amcs+bmds;
			exp2 = epsilon+cmes+dmfs;
			exp3 = sqrt(exp1);
			exp4 = sqrt(exp2);
			exp5 = pow(exp1,3.0/2.0);
			exp6 = pow(exp2,3.0/2.0);
			exp7 = 1.0/exp3;
			exp8 = 1.0/exp4;
			exp9 = 1.0/exp5;
			exp10 = 1.0/exp6;

		}
		
		ff2=1./(1+exp(-m2 + so));
		ff1=1./(1+exp(-m1 + so));
		
		ff1b=1-ff1;
		ff2b=1-ff2;
		for (j = 0; j < Ndet; j++)
		{
			indj=(j*Fdet)+t;
			det2x = Xd[indj];

			if (det2x != 0)
			{
				det2y = Yd[indj];					
				dist2 = sqrt((c-det2x)*(c-det2x) + (d-det2y)*(d-det2y));
				
				
					if (dist2 < csig*cutoffthr) {
						d2x = Dx[indj];d2y = Dy[indj];
					
						tmp2=csig / (dist2*dist2 + csig);
						f2 = ff2*tmp2;
						f2b=ff2b * (csig - tmp2);
						/*g2 = acos( (emc*d2x + fmd*d2y) / m2);
						g2 = (-(g2-3.141592653589793/2)*(g2-3.141592653589793/2) + (3.141592653589793*3.141592653589793)/4);
						g2 = (wt*exp(kappa*(o2[0]*d2x + o2[1]*d2y)) + (1-wt)*exp(kappa*(-o2[0]*d2x - o2[1]*d2y)));
						*/
						g2 = exp(kappa*(o2[0]*d2x + o2[1]*d2y));
						g2=-log(g2/intv);
					
						/*f2acos = f2 * g2;*/
						/*printf("a %i %i %i %i %f %f %f %f\n",i,t,j,k,f2,exp(-m2+so),(dist2*dist2+csig),f2acos);*/
						if (nlhs>1) {
							

						}
					}
					for (k = 0; k < Ndet; k++) {
						indk=(k*Fdet)+t-1;
						det1x = Xd[indk];

						if (det1x != 0)
						{
							det1y = Yd[indk];			
							
							dist1 = sqrt((a-det1x)*(a-det1x) + (b-det1y)*(b-det1y));
							
							if (dist2 < csig*cutoffthr || dist1 < csig*cutoffthr) {
								d1x = Dx[indk];d1y = Dy[indk];
								
							
								tmp1=csig/(dist1*dist1+csig);
								f1 = ff1*tmp1;  
								f1b = ff1b * (csig-tmp1);
								/*g1=acos((cma*d1x + dmb*d1y) / m1);
 								g1 = (-(g1-3.141592653589793/2)*(g1-3.141592653589793/2) + (3.141592653589793*3.141592653589793)/4);
 								g1 = (wt*exp(kappa*(o1[0]*d1x + o1[1]*d1y)) + (1-wt)*exp(kappa*(-o1[0]*d1x - o1[1]*d1y)));
								
								*/
								g1 = exp(kappa*(o1[0]*d1x + o1[1]*d1y));
								g1=-log(g1/intv);
								

								/*f1acos=f1 * g1;*/
								
								
	/* 							printf("b %i %i %i %i %f %f %f %f %f %f %f %f\n",i,t,j,k,f1acos,f2acos,cma,d1x,dmb,d1y,m1,f1); */
								/*mexCallMATLAB(0, NULL, 1, &prhsp, "pause");*/
								diff = f1*g1 + f2*g2 + f1b*distrconst + f2b*distrconst;
								
								/*
								if (t==2 && i == 1) {
								printf("%f %f %f %f  %f %f  %f %f\n",o1[0],o1[1],o2[0],o2[1],d1x,d1y,d2x,d2y);
								printf("%f %f %f %f %f\n",f1,f2,g1,g2,diff);
								}
								*/
								
                                /*
								if (isnan(diff)) printf("AAA"); 
								if (isinf(diff)) printf("AAA"); 
								if (diff<0) printf("AHA!");
                                 **/
								*fx += diff;
								if (nlhs>1)
								{
									

									exp11 = pow(a-det1x,2.0);
									exp12 = pow(b-det1y,2.0);
									exp13 = pow(c-det2x,2.0);
									exp14 = pow(d-det2y,2.0);
									exp15 = exp(so-exp3);
									exp16 = exp(so-exp4);
									exp17 = d1x*amc;
									exp18 = d1y*bmd;
									exp19 = d2x*cme;
									exp20 = d2y*dmf;
									exp21 = (exp17*exp7+exp18*exp7);
									exp22 = (exp19*exp8+exp20*exp8);
									exp23 = 1.0/pow(csig+exp11+exp12,2.0);
									exp24 = 1.0/pow(exp16+1.0,2.0);
									exp25 = exp(-kappa*exp21);
									exp26 = exp(kappa*exp21);
									exp27 = exp(-kappa*exp22);
									exp28 = exp(kappa*exp22);
									exp29 = (wt-1.0);
									exp30 = 1.0/pow(exp15+1.0,2.0);
			exp31=log(-(exp26*exp29-wt*exp25)/intv);
			exp32=log(-(exp28*exp29-wt*exp27)/intv);
									/*
									exp13=d1x*amc;								exp14=d1y*bmd;
									exp19=1.0/sqrt(-pow(exp13+exp14,2.0)/(epsilon+exp1+exp2)+1.0);
									exp21=acos((exp13+exp14)*1.0/exp5);
									
									
									diff1 = (csig*(d1x*1.0/exp5-(a*2.0-c*2.0)*(exp13+exp14)*1.0/exp6*(1.0/2.0))*exp19)/((exp15+1.0)*(csig+exp7+exp9))-(csig*(3.141592653589793-exp21)*(a*2.0-det1x*2.0)*1.0/pow(csig+exp7+exp9,2.0))/(exp15+1.0)+(csig*exp15*(3.141592653589793-exp21)*exp23*(a*2.0-c*2.0)*1.0/exp5*(1.0/2.0))/(csig+exp7+exp9);
									diff2 = (csig*(d1y*1.0/exp5-(b*2.0-d*2.0)*(exp13+exp14)*1.0/exp6*(1.0/2.0))*exp19)/((exp15+1.0)*(csig+exp7+exp9))-(csig*(3.141592653589793-exp21)*(b*2.0-det1y*2.0)*1.0/pow(csig+exp7+exp9,2.0))/(exp15+1.0)+(csig*exp15*(3.141592653589793-exp21)*exp23*(b*2.0-d*2.0)*1.0/exp5*(1.0/2.0))/(csig+exp7+exp9);
									diff3= -(csig*(d1x*1.0/exp5-(a*2.0-c*2.0)*(exp13+exp14)*1.0/exp6*(1.0/2.0))*exp19)/((exp15+1.0)*(csig+exp7+exp9))-(csig*(3.141592653589793-exp22)*(c*2.0-det2x*2.0)*1.0/pow(csig+exp8+exp10,2.0))/(exp16+1.0)+(csig*exp17*(d2x*exp18-(c*2.0-e*2.0)*(exp11+exp12)*1.0/exp20*(1.0/2.0)))/((exp16+1.0)*(csig+exp8+exp10))+(csig*exp16*(3.141592653589793-exp22)*(c*2.0-e*2.0)*exp24*exp18*(1.0/2.0))/(csig+exp8+exp10)-(csig*exp15*(3.141592653589793-exp21)*exp23*(a*2.0-c*2.0)*1.0/exp5*(1.0/2.0))/(csig+exp7+exp9);
									diff4= -(csig*(d1y*1.0/exp5-(b*2.0-d*2.0)*(exp13+exp14)*1.0/exp6*(1.0/2.0))*exp19)/((exp15+1.0)*(csig+exp7+exp9))-(csig*(3.141592653589793-exp22)*(d*2.0-det2y*2.0)*1.0/pow(csig+exp8+exp10,2.0))/(exp16+1.0)+(csig*exp17*(d2y*exp18-(d*2.0-f*2.0)*(exp11+exp12)*1.0/exp20*(1.0/2.0)))/((exp16+1.0)*(csig+exp8+exp10))+(csig*exp16*(3.141592653589793-exp22)*(d*2.0-f*2.0)*exp24*exp18*(1.0/2.0))/(csig+exp8+exp10)-(csig*exp15*(3.141592653589793-exp21)*exp23*(b*2.0-d*2.0)*1.0/exp5*(1.0/2.0))/(csig+exp7+exp9);
									diff5= -(csig*exp17*(d2x*exp18-(c*2.0-e*2.0)*(exp11+exp12)*1.0/exp20*(1.0/2.0)))/((exp16+1.0)*(csig+exp8+exp10))-(csig*exp16*(3.141592653589793-exp22)*(c*2.0-e*2.0)*exp24*exp18*(1.0/2.0))/(csig+exp8+exp10);
									diff6= -(csig*exp17*(d2y*exp18-(d*2.0-f*2.0)*(exp11+exp12)*1.0/exp20*(1.0/2.0)))/((exp16+1.0)*(csig+exp8+exp10))-(csig*exp16*(3.141592653589793-exp22)*(d*2.0-f*2.0)*exp24*exp18*(1.0/2.0))/(csig+exp8+exp10);
									*/
	
									/*
diff1 = (csig*(kappa*wt*exp25*(-d1x*exp7+exp17*a2c2*exp9*(1.0/2.0)+exp18*a2c2*exp9*(1.0/2.0))+kappa*exp26*exp29*(-d1x*exp7+exp17*a2c2*exp9*(1.0/2.0)+exp18*a2c2*exp9*(1.0/2.0))))/((exp15+1.0)*(csig+exp11+exp12))+(csig*(exp26*exp29-wt*exp25)*(a*2.0-det1x*2.0)*exp23)/(exp15+1.0)-(csig*exp15*(exp26*exp29-wt*exp25)*exp30*a2c2*exp7*(1.0/2.0))/(csig+exp11+exp12);;
diff2 = (csig*(kappa*wt*exp25*(-d1y*exp7+exp17*b2d2*exp9*(1.0/2.0)+exp18*b2d2*exp9*(1.0/2.0))+kappa*exp26*exp29*(-d1y*exp7+exp17*b2d2*exp9*(1.0/2.0)+exp18*b2d2*exp9*(1.0/2.0))))/((exp15+1.0)*(csig+exp11+exp12))+(csig*(exp26*exp29-wt*exp25)*(b*2.0-det1y*2.0)*exp23)/(exp15+1.0)-(csig*exp15*(exp26*exp29-wt*exp25)*exp30*b2d2*exp7*(1.0/2.0))/(csig+exp11+exp12);;
diff3 = (csig*(kappa*exp28*exp29*(-d2x*exp8+exp19*c2e2*exp10*(1.0/2.0)+exp20*c2e2*exp10*(1.0/2.0))+kappa*wt*exp27*(-d2x*exp8+exp19*c2e2*exp10*(1.0/2.0)+exp20*c2e2*exp10*(1.0/2.0))))/((exp16+1.0)*(csig+exp13+exp14))-(csig*(kappa*wt*exp25*(-d1x*exp7+exp17*a2c2*exp9*(1.0/2.0)+exp18*a2c2*exp9*(1.0/2.0))+kappa*exp26*exp29*(-d1x*exp7+exp17*a2c2*exp9*(1.0/2.0)+exp18*a2c2*exp9*(1.0/2.0))))/((exp15+1.0)*(csig+exp11+exp12))+(csig*(exp28*exp29-wt*exp27)*(c*2.0-det2x*2.0)*1.0/pow(csig+exp13+exp14,2.0))/(exp16+1.0)+(csig*exp15*(exp26*exp29-wt*exp25)*exp30*a2c2*exp7*(1.0/2.0))/(csig+exp11+exp12)-(csig*exp16*(exp28*exp29-wt*exp27)*c2e2*exp24*exp8*(1.0/2.0))/(csig+exp13+exp14);;
diff4 = (csig*(kappa*exp28*exp29*(-d2y*exp8+exp19*d2f2*exp10*(1.0/2.0)+exp20*d2f2*exp10*(1.0/2.0))+kappa*wt*exp27*(-d2y*exp8+exp19*d2f2*exp10*(1.0/2.0)+exp20*d2f2*exp10*(1.0/2.0))))/((exp16+1.0)*(csig+exp13+exp14))-(csig*(kappa*wt*exp25*(-d1y*exp7+exp17*b2d2*exp9*(1.0/2.0)+exp18*b2d2*exp9*(1.0/2.0))+kappa*exp26*exp29*(-d1y*exp7+exp17*b2d2*exp9*(1.0/2.0)+exp18*b2d2*exp9*(1.0/2.0))))/((exp15+1.0)*(csig+exp11+exp12))+(csig*(exp28*exp29-wt*exp27)*(d*2.0-det2y*2.0)*1.0/pow(csig+exp13+exp14,2.0))/(exp16+1.0)+(csig*exp15*(exp26*exp29-wt*exp25)*exp30*b2d2*exp7*(1.0/2.0))/(csig+exp11+exp12)-(csig*exp16*(exp28*exp29-wt*exp27)*d2f2*exp24*exp8*(1.0/2.0))/(csig+exp13+exp14);;
diff5 = -(csig*(kappa*exp28*exp29*(-d2x*exp8+exp19*c2e2*exp10*(1.0/2.0)+exp20*c2e2*exp10*(1.0/2.0))+kappa*wt*exp27*(-d2x*exp8+exp19*c2e2*exp10*(1.0/2.0)+exp20*c2e2*exp10*(1.0/2.0))))/((exp16+1.0)*(csig+exp13+exp14))+(csig*exp16*(exp28*exp29-wt*exp27)*c2e2*exp24*exp8*(1.0/2.0))/(csig+exp13+exp14);;
diff6 = -(csig*(kappa*exp28*exp29*(-d2y*exp8+exp19*d2f2*exp10*(1.0/2.0)+exp20*d2f2*exp10*(1.0/2.0))+kappa*wt*exp27*(-d2y*exp8+exp19*d2f2*exp10*(1.0/2.0)+exp20*d2f2*exp10*(1.0/2.0))))/((exp16+1.0)*(csig+exp13+exp14))+(csig*exp16*(exp28*exp29-wt*exp27)*d2f2*exp24*exp8*(1.0/2.0))/(csig+exp13+exp14);;

diff1 = (csig*(kappa*wt*exp25*(-d1x*exp7+exp17*a2c2*exp9*(1.0/2.0)+exp18*a2c2*exp9*(1.0/2.0))+kappa*exp26*exp29*(-d1x*exp7+exp17*a2c2*exp9*(1.0/2.0)+exp18*a2c2*exp9*(1.0/2.0))))/((exp26*exp29-wt*exp25)*(exp15+1.0)*(csig+exp11+exp12))+(csig*exp31*(a*2.0-det1x*2.0)*exp23)/(exp15+1.0)-(csig*exp15*exp31*exp30*a2c2*exp7*(1.0/2.0))/(csig+exp11+exp12);
diff2 = (csig*(kappa*wt*exp25*(-d1y*exp7+exp17*b2d2*exp9*(1.0/2.0)+exp18*b2d2*exp9*(1.0/2.0))+kappa*exp26*exp29*(-d1y*exp7+exp17*b2d2*exp9*(1.0/2.0)+exp18*b2d2*exp9*(1.0/2.0))))/((exp26*exp29-wt*exp25)*(exp15+1.0)*(csig+exp11+exp12))+(csig*exp31*(b*2.0-det1y*2.0)*exp23)/(exp15+1.0)-(csig*exp15*exp31*exp30*b2d2*exp7*(1.0/2.0))/(csig+exp11+exp12);
diff3 = -(csig*(kappa*wt*exp25*(-d1x*exp7+exp17*a2c2*exp9*(1.0/2.0)+exp18*a2c2*exp9*(1.0/2.0))+kappa*exp26*exp29*(-d1x*exp7+exp17*a2c2*exp9*(1.0/2.0)+exp18*a2c2*exp9*(1.0/2.0))))/((exp26*exp29-wt*exp25)*(exp15+1.0)*(csig+exp11+exp12))+(csig*exp32*(c*2.0-det2x*2.0)*1.0/pow(csig+exp13+exp14,2.0))/(exp16+1.0)+(csig*(kappa*exp28*exp29*(-d2x*exp8+exp19*c2e2*exp10*(1.0/2.0)+exp20*c2e2*exp10*(1.0/2.0))+kappa*wt*exp27*(-d2x*exp8+exp19*c2e2*exp10*(1.0/2.0)+exp20*c2e2*exp10*(1.0/2.0))))/((exp28*exp29-wt*exp27)*(exp16+1.0)*(csig+exp13+exp14))+(csig*exp15*exp31*exp30*a2c2*exp7*(1.0/2.0))/(csig+exp11+exp12)-(csig*exp32*exp16*c2e2*exp24*exp8*(1.0/2.0))/(csig+exp13+exp14);
diff4 = (csig*(kappa*exp28*exp29*(-d2y*exp8+exp19*d2f2*exp10*(1.0/2.0)+exp20*d2f2*exp10*(1.0/2.0))+kappa*wt*exp27*(-d2y*exp8+exp19*d2f2*exp10*(1.0/2.0)+exp20*d2f2*exp10*(1.0/2.0))))/((exp28*exp29-wt*exp27)*(exp16+1.0)*(csig+exp13+exp14))-(csig*(kappa*wt*exp25*(-d1y*exp7+exp17*b2d2*exp9*(1.0/2.0)+exp18*b2d2*exp9*(1.0/2.0))+kappa*exp26*exp29*(-d1y*exp7+exp17*b2d2*exp9*(1.0/2.0)+exp18*b2d2*exp9*(1.0/2.0))))/((exp26*exp29-wt*exp25)*(exp15+1.0)*(csig+exp11+exp12))+(csig*exp32*(d*2.0-det2y*2.0)*1.0/pow(csig+exp13+exp14,2.0))/(exp16+1.0)+(csig*exp15*exp31*exp30*b2d2*exp7*(1.0/2.0))/(csig+exp11+exp12)-(csig*exp32*exp16*d2f2*exp24*exp8*(1.0/2.0))/(csig+exp13+exp14);
diff5 = -(csig*(kappa*exp28*exp29*(-d2x*exp8+exp19*c2e2*exp10*(1.0/2.0)+exp20*c2e2*exp10*(1.0/2.0))+kappa*wt*exp27*(-d2x*exp8+exp19*c2e2*exp10*(1.0/2.0)+exp20*c2e2*exp10*(1.0/2.0))))/((exp28*exp29-wt*exp27)*(exp16+1.0)*(csig+exp13+exp14))+(csig*exp32*exp16*c2e2*exp24*exp8*(1.0/2.0))/(csig+exp13+exp14);
diff6 = -(csig*(kappa*exp28*exp29*(-d2y*exp8+exp19*d2f2*exp10*(1.0/2.0)+exp20*d2f2*exp10*(1.0/2.0))+kappa*wt*exp27*(-d2y*exp8+exp19*d2f2*exp10*(1.0/2.0)+exp20*d2f2*exp10*(1.0/2.0))))/((exp28*exp29-wt*exp27)*(exp16+1.0)*(csig+exp13+exp14))+(csig*exp32*exp16*d2f2*exp24*exp8*(1.0/2.0))/(csig+exp13+exp14);
*/									
diff1 = (csig*log(exp25/intv)*(a*2.0-det1x*2.0)*exp23)/(exp15+1.0)-(csig*kappa*(-d1x*exp7+exp17*a2c2*exp9*(1.0/2.0)+exp18*a2c2*exp9*(1.0/2.0)))/((exp15+1.0)*(csig+exp11+exp12))-csig*distrconst*(1.0/(exp15+1.0)-1.0)*(a*2.0-det1x*2.0)*exp23-distrconst*exp15*(csig-csig/(csig+exp11+exp12))*exp30*a2c2*exp7*(1.0/2.0)-(csig*exp15*log(exp25/intv)*exp30*a2c2*exp7*(1.0/2.0))/(csig+exp11+exp12);
diff2 = (csig*log(exp25/intv)*(b*2.0-det1y*2.0)*exp23)/(exp15+1.0)-(csig*kappa*(-d1y*exp7+exp17*b2d2*exp9*(1.0/2.0)+exp18*b2d2*exp9*(1.0/2.0)))/((exp15+1.0)*(csig+exp11+exp12))-csig*distrconst*(1.0/(exp15+1.0)-1.0)*(b*2.0-det1y*2.0)*exp23-distrconst*exp15*(csig-csig/(csig+exp11+exp12))*exp30*b2d2*exp7*(1.0/2.0)-(csig*exp15*log(exp25/intv)*exp30*b2d2*exp7*(1.0/2.0))/(csig+exp11+exp12);
diff3 = (csig*log(exp27/intv)*(c*2.0-det2x*2.0)*1.0/pow(csig+exp13+exp14,2.0))/(exp16+1.0)+(csig*kappa*(-d1x*exp7+exp17*a2c2*exp9*(1.0/2.0)+exp18*a2c2*exp9*(1.0/2.0)))/((exp15+1.0)*(csig+exp11+exp12))-csig*distrconst*(1.0/(exp16+1.0)-1.0)*(c*2.0-det2x*2.0)*1.0/pow(csig+exp13+exp14,2.0)-(csig*kappa*(-d2x*exp8+exp19*c2e2*exp10*(1.0/2.0)+exp20*c2e2*exp10*(1.0/2.0)))/((exp16+1.0)*(csig+exp13+exp14))+distrconst*exp15*(csig-csig/(csig+exp11+exp12))*exp30*a2c2*exp7*(1.0/2.0)-distrconst*exp16*(csig-csig/(csig+exp13+exp14))*c2e2*exp24*exp8*(1.0/2.0)-(csig*exp16*log(exp27/intv)*c2e2*exp24*exp8*(1.0/2.0))/(csig+exp13+exp14)+(csig*exp15*log(exp25/intv)*exp30*a2c2*exp7*(1.0/2.0))/(csig+exp11+exp12);
diff4 = (csig*log(exp27/intv)*(d*2.0-det2y*2.0)*1.0/pow(csig+exp13+exp14,2.0))/(exp16+1.0)+(csig*kappa*(-d1y*exp7+exp17*b2d2*exp9*(1.0/2.0)+exp18*b2d2*exp9*(1.0/2.0)))/((exp15+1.0)*(csig+exp11+exp12))-csig*distrconst*(1.0/(exp16+1.0)-1.0)*(d*2.0-det2y*2.0)*1.0/pow(csig+exp13+exp14,2.0)-(csig*kappa*(-d2y*exp8+exp19*d2f2*exp10*(1.0/2.0)+exp20*d2f2*exp10*(1.0/2.0)))/((exp16+1.0)*(csig+exp13+exp14))+distrconst*exp15*(csig-csig/(csig+exp11+exp12))*exp30*b2d2*exp7*(1.0/2.0)-distrconst*exp16*(csig-csig/(csig+exp13+exp14))*d2f2*exp24*exp8*(1.0/2.0)-(csig*exp16*log(exp27/intv)*d2f2*exp24*exp8*(1.0/2.0))/(csig+exp13+exp14)+(csig*exp15*log(exp25/intv)*exp30*b2d2*exp7*(1.0/2.0))/(csig+exp11+exp12);
diff5 = (csig*kappa*(-d2x*exp8+exp19*c2e2*exp10*(1.0/2.0)+exp20*c2e2*exp10*(1.0/2.0)))/((exp16+1.0)*(csig+exp13+exp14))+distrconst*exp16*(csig-csig/(csig+exp13+exp14))*c2e2*exp24*exp8*(1.0/2.0)+(csig*exp16*log(exp27/intv)*c2e2*exp24*exp8*(1.0/2.0))/(csig+exp13+exp14);
diff6 = (csig*kappa*(-d2y*exp8+exp19*d2f2*exp10*(1.0/2.0)+exp20*d2f2*exp10*(1.0/2.0)))/((exp16+1.0)*(csig+exp13+exp14))+distrconst*exp16*(csig-csig/(csig+exp13+exp14))*d2f2*exp24*exp8*(1.0/2.0)+(csig*exp16*log(exp27/intv)*d2f2*exp24*exp8*(1.0/2.0))/(csig+exp13+exp14);


									/*
diff1 = (csig*log(exp25/intv)*(a*2.0-det1x*2.0)*exp23)/(exp15+1.0)-(csig*kappa*(-d1x*exp7+exp17*a2c2*exp9*(1.0/2.0)+exp18*a2c2*exp9*(1.0/2.0)))/((exp15+1.0)*(csig+exp11+exp12))-(csig*exp15*log(exp25/intv)*exp30*a2c2*exp7*(1.0/2.0))/(csig+exp11+exp12);
diff2 = (csig*log(exp25/intv)*(b*2.0-det1y*2.0)*exp23)/(exp15+1.0)-(csig*kappa*(-d1y*exp7+exp17*b2d2*exp9*(1.0/2.0)+exp18*b2d2*exp9*(1.0/2.0)))/((exp15+1.0)*(csig+exp11+exp12))-(csig*exp15*log(exp25/intv)*exp30*b2d2*exp7*(1.0/2.0))/(csig+exp11+exp12);
diff3 = (csig*log(exp27/intv)*(c*2.0-det2x*2.0)*1.0/pow(csig+exp13+exp14,2.0))/(exp16+1.0)+(csig*kappa*(-d1x*exp7+exp17*a2c2*exp9*(1.0/2.0)+exp18*a2c2*exp9*(1.0/2.0)))/((exp15+1.0)*(csig+exp11+exp12))-(csig*kappa*(-d2x*exp8+exp19*c2e2*exp10*(1.0/2.0)+exp20*c2e2*exp10*(1.0/2.0)))/((exp16+1.0)*(csig+exp13+exp14))-(csig*exp16*log(exp27/intv)*c2e2*exp24*exp8*(1.0/2.0))/(csig+exp13+exp14)+(csig*exp15*log(exp25/intv)*exp30*a2c2*exp7*(1.0/2.0))/(csig+exp11+exp12);
diff4 = (csig*log(exp27/intv)*(d*2.0-det2y*2.0)*1.0/pow(csig+exp13+exp14,2.0))/(exp16+1.0)+(csig*kappa*(-d1y*exp7+exp17*b2d2*exp9*(1.0/2.0)+exp18*b2d2*exp9*(1.0/2.0)))/((exp15+1.0)*(csig+exp11+exp12))-(csig*kappa*(-d2y*exp8+exp19*d2f2*exp10*(1.0/2.0)+exp20*d2f2*exp10*(1.0/2.0)))/((exp16+1.0)*(csig+exp13+exp14))-(csig*exp16*log(exp27/intv)*d2f2*exp24*exp8*(1.0/2.0))/(csig+exp13+exp14)+(csig*exp15*log(exp25/intv)*exp30*b2d2*exp7*(1.0/2.0))/(csig+exp11+exp12);
diff5 = (csig*kappa*(-d2x*exp8+exp19*c2e2*exp10*(1.0/2.0)+exp20*c2e2*exp10*(1.0/2.0)))/((exp16+1.0)*(csig+exp13+exp14))+(csig*exp16*log(exp27/intv)*c2e2*exp24*exp8*(1.0/2.0))/(csig+exp13+exp14);
diff6 = (csig*kappa*(-d2y*exp8+exp19*d2f2*exp10*(1.0/2.0)+exp20*d2f2*exp10*(1.0/2.0)))/((exp16+1.0)*(csig+exp13+exp14))+(csig*exp16*log(exp27/intv)*d2f2*exp24*exp8*(1.0/2.0))/(csig+exp13+exp14);
*/
/*
									if (isnan(diff1) || isinf(diff1)) printf("AAA"); 
									if (isnan(diff2) || isinf(diff2)) printf("AAA"); 
									if (isnan(diff3) || isinf(diff3)) printf("AAA"); 
									if (isnan(diff4) || isinf(diff4)) printf("AAA"); 
									if (isnan(diff5) || isinf(diff5)) printf("AAA"); 
									if (isnan(diff6) || isinf(diff6)) printf("AAA"); 
 */
									
									df[xind-2] += diff1;
									df[yind-2] += diff2;
									df[xind] += diff3;
									df[yind] += diff4;
									df[xind+2] += diff5;
									df[yind+2] += diff6;
									
								}
							}


							/* det exists */
							
							
						}
					}
				
			}
		}
		xind += 2;
		yind = xind+1;
		}
		if (tlength>1) {
			xind += 2;
		}	
	}
}
