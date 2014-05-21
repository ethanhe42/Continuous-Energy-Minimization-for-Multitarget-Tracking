/* (C) Anton Andriyenko, 2012 */

#include "mex.h"
#include <math.h>

/* // WRONG  GRADIENT !!!  if (ci==1 || ci==3) { si*=-1;}*/
void min_distances(const double *areaLimits, double x, double y, double* dm, double* ci, int* si)
{
  const double dl = fabs(areaLimits[0] - x);
  const double dr = fabs(areaLimits[1] - x);
  const double du = fabs(areaLimits[2] - y);
  const double dd = fabs(areaLimits[3] - y);

  double dm1 = (dl < dr) ? dl : dr;
  double dm2 = (du < dd) ? du : dd;
  double ci1 = (dl < dr) ? 0 : 1;
  double ci2 = (du < dd) ? 2 : 3;

  *dm = (dm1 < dm2) ? dm1 : dm2;
  *ci = (dm1 < dm2) ? ci1 : ci2;

	*si=1;
	if (x<areaLimits[0] || x>areaLimits[1] || y<areaLimits[2] || y>areaLimits[3]) {
		*si=-1;
	}

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {


    /* /Declarations */
    const mxArray *Xmtx, *Ymtx, *areaLimitsdata, *gridStepdata,  *tEdata, *xldata;

    double *X, *Y, *areaLimits, gridStep, *targetsExist, xl;

    int i,t;
    int F,N;

    double *fx;
    double *dfx;
    double scale;

	int st, en, xind, yind, ind;

    /* //Copy input pointer x */
    Xmtx = prhs[0];
	Ymtx = prhs[1];
	areaLimitsdata = prhs[2];
    gridStepdata = prhs[3];
    tEdata = prhs[4];
	xldata = prhs[5];

    /* //Get matrix x */
    X = mxGetPr(Xmtx);
	Y = mxGetPr(Ymtx);
    N = mxGetN(Xmtx);
    F = mxGetM(Xmtx);

	areaLimits = mxGetPr(areaLimitsdata);
    gridStep = (double)(mxGetScalar(gridStepdata));
	targetsExist = mxGetPr(tEdata);
	xl = (int)(mxGetScalar(xldata));


    /* Allocate memory and assign output pointer */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(xl, 1, mxREAL);

    /* Get a pointer to the data space in our newly allocated memory */
    fx = mxGetPr(plhs[0]);
    dfx = mxGetPr(plhs[1]);
	*fx = 0;
	scale = 1.0 / gridStep;

	xind=0;
	yind=0;

/* for each target */
  for (i = 0; i < N; i++)
  {
	int st = (int) targetsExist[i]-1;
    int en = (int) targetsExist[i+N]-1;

	if (st > 0)
    {
		double xp,yp,dm,ci;
		double d,tmp,val,tmp2;
		int si;
	  ind=(i*F)+st;
      xp = X[ind];
      yp = Y[ind];
      dm = 0;
      ci = 0;
	  si = 1;
      min_distances(areaLimits, xp, yp, &dm, &ci, &si);
	  if (ci==1 || ci==3) { si*=-1;}
      d = 1.0 / (1 + exp(2 - scale*dm)) + 1.0 / (1 + exp(2 + scale*dm)) ;
      *fx = *fx + d;
      tmp = exp(2 - scale * dm) + 1;      tmp *= tmp;
	  tmp2 = exp(2 + scale * dm) + 1;      tmp2 *= tmp2;
	  
      val = si*scale * exp(2 - scale * dm) / tmp + si*scale * exp(2 + scale * dm) / tmp2;

      if (ci < 2)
      {
        dfx[xind] = val;
      }
      else
      {
        dfx[xind+1] = val;
      }
    }
    xind+=(en-st)*2;

	if (en < F-1)
    {
		double xp,yp,dm,ci;
		double d,tmp,val,tmp2;
		int si;

	  ind=(i*F)+en;
      xp = X[ind];
      yp = Y[ind];
      dm = 0;
      ci = 0;
	  si = 1;
      min_distances(areaLimits, xp, yp, &dm, &ci,&si);
	  if (ci==1 || ci==3) { si*=-1;}
      d = 1.0 / (1 + exp(2 - scale*dm)) + 1.0 / (1 + exp(2 + scale*dm));
      *fx = *fx + d;
      tmp = exp(2 - scale * dm) + 1;      tmp *= tmp;
      tmp2 = exp(2 + scale * dm) + 1;      tmp2 *= tmp2;
	  
      val = si*scale * exp(2 - scale * dm) / tmp + si*scale * exp(2 + scale * dm) / tmp2;
      if (ci < 2)
      {
        dfx[xind] = val;
      }
      else
      {
        dfx[xind+1] = val;
      }
    }
    xind += 2;
  }

}