/* (C) Anton Andriyenko, 2012

 The code may be used free of charge for non-commercial and
 educational purposes, the only requirement is that this text is
 preserved within the derivative work. For any other purpose you
 must contact the authors for permission. This code may not be
 redistributed without written permission from the authors.
*/

#include "mex.h"
#include <math.h>

void dfnddALL(double *prhs, double *plhs) {


    /* /Declarations */
    const mxArray *Dpxdata, *Dpydata, *Sxdata, *Xwadata, *Xwbdata, *Ywadata, *Ywbdata, *Zwadata, *Zwbdata, *focaldata, \
            *mR11data, *mR12data, *mR13data, *mR21data, *mR22data, *mR23data, *mR31data, *mR32data, *mR33data, \
            *txdata, *tydata, *tzdata, *txidata, *tyidata, *tzidata;


    double Dpx, Dpy, Sx, Xwa, Xwb, Ywa, Ywb, Zwa, Zwb, focal, \
            mR11, mR12, mR13, mR21, mR22, mR23, mR31, mR32, mR33, \
            tx, ty, tz, txi, tyi, tzi;

    double Xa, Ya, Xb, Yb;

double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10,\
 t11, t12, t13, t14, t15, t16, t17, t18, t19, t20,\
 t21, t22, t23, t24, t25, t26, t27, t28, t29, t30,\
 t31, t32, t33, t34, t35, t36, t37, t38, t39, t40,\
 t41, t42, t43, t44, t45, t46, t47, t48, t49, t50,\
 t51, t52, t53, t54, t55, t56, t57, t58, t59, t60,\
 t61, t62, t63, t64, t65, t66, t67, t68, t69, t70,\
 t71, t72, t73, t74, t75, t76, t77, t78, t79, t80,\
 t81, t82, t83, t84, t85, t86, t87, t88, t89, t90,\
 t91, t92, t93, t94, t95, t96, t97, t98, t99, t100,\
 t101, t102, t103, t104, t105, t106, t107, t108, t109, t110,\
 t111, t112, t113, t114, t115, t116, t117, t118, t119, t120,\
 t121, t122, t123, t124, t125, t126, t127, t128, t129, t130,\
 t131, t132, t133, t134, t135, t136, t137, t138, t139, t140,\
 t141, t142, t143, t144, t145, t146, t147, t148, t149, t150,\
 t151, t152, t153, t154, t155, t156, t157, t158, t159, t160,\
 t161, t162, t163, t164, t165, t166, t167, t168, t169, t170,\
 t171, t172, t173, t174, t175, t176, t177, t178, t179, t180,\
 t181;

    /* //Copy input pointer x */
    Dpx= prhs[0];
    Dpy= prhs[1];
    Sx= prhs[2];
    Xwa= prhs[3];
    Xwb= prhs[4];
    Ywa= prhs[5];
    Ywb= prhs[6];
    Zwa= prhs[7];
    Zwb= prhs[8];
    focal= prhs[9];
    mR11= prhs[10];
    mR12= prhs[11];
    mR13= prhs[12];
    mR21= prhs[13];
    mR22= prhs[14];
    mR23= prhs[15];
    mR31= prhs[16];
    mR32= prhs[17];
    mR33= prhs[18];
    tx= prhs[19];
    ty= prhs[20];
    tz= prhs[21];
    txi = prhs[22];
    tyi = prhs[23];
    tzi = prhs[24];




t2 = 1/Dpx;
t15 = Xwa*mR31;
t16 = Ywa*mR32;
t17 = Zwa*mR33;
t18 = t15+t16+t17+tz;
t19 = 1/t18;
t21 = Xwb*mR31;
t22 = Ywb*mR32;
t23 = Zwb*mR33;
t24 = t21+t22+t23+tz;
t25 = 1/t24;
t83 = Xwa*mR11;
t84 = Ywa*mR12;
t85 = Zwa*mR13;
t86 = t83+t84+t85+tx;
t87 = Sx*focal*t19*t2*t86;
t88 = Xwb*mR11;
t89 = Ywb*mR12;
t90 = Zwb*mR13;
t91 = t88+t89+t90+tx;
t92 = Sx*focal*t2*t25*t91;
t3 = t87-t92;
t4 = pow(focal,2);
t5 = pow(txi,2);
t6 = pow(Dpy,2);
t7 = pow(tyi,2);
t8 = pow(tzi,2);
t9 = pow(Xwa,2);
t10 = pow(Xwb,2);
t11 = pow(Ywb,2);
t12 = pow(Ywa,2);
t13 = pow(Zwb,2);
t14 = pow(Zwa,2);
t20 = 1/Dpy;
t73 = Xwa*mR21;
t74 = Ywa*mR22;
t75 = Zwa*mR23;
t76 = t73+t74+t75+ty;
t77 = focal*t19*t20*t76;
t78 = Xwb*mR21;
t79 = Ywb*mR22;
t80 = Zwb*mR23;
t81 = t78+t79+t80+ty;
t82 = focal*t20*t25*t81;
t26 = t77-t82;
t27 = pow(t5,2);
t28 = t27*t6;
t29 = pow(t7,2);
t30 = t29*t6;
t31 = pow(t8,2);
t32 = t31*t6;
t33 = t10*t6*t9;
t34 = t11*t6*t9;
t35 = t10*t12*t6;
t36 = t13*t6*t9;
t37 = t10*t14*t6;
t38 = t11*t12*t6;
t39 = t12*t13*t6;
t40 = t11*t14*t6;
t41 = t13*t14*t6;
t42 = t5*t6*t9;
t43 = t10*t5*t6;
t44 = t6*t7*t9;
t45 = t10*t6*t7;
t46 = t6*t8*t9;
t47 = t10*t6*t8;
t48 = t12*t5*t6;
t49 = t11*t5*t6;
t50 = t12*t6*t7;
t51 = t11*t6*t7;
t52 = t12*t6*t8;
t53 = t11*t6*t8;
t54 = t14*t5*t6;
t55 = t13*t5*t6;
t56 = t14*t6*t7;
t57 = t13*t6*t7;
t58 = t14*t6*t8;
t59 = t13*t6*t8;
t60 = t5*t6*t7*2.0;
t61 = t5*t6*t8*2.0;
t62 = t6*t7*t8*2.0;
t63 = Xwa*Xwb*t5*t6*4.0;
t64 = Ywa*Ywb*t6*t7*4.0;
t65 = Zwa*Zwb*t6*t8*4.0;
t66 = Xwa*Ywb*t6*txi*tyi*4.0;
t67 = Xwb*Ywa*t6*txi*tyi*4.0;
t68 = Xwa*Zwb*t6*txi*tzi*4.0;
t69 = Xwb*Zwa*t6*txi*tzi*4.0;
t70 = Ywa*Zwb*t6*tyi*tzi*4.0;
t71 = Ywb*Zwa*t6*tyi*tzi*4.0;
t133 = Xwa*t5*t6*txi*2.0;
t134 = Xwb*t5*t6*txi*2.0;
t135 = Ywa*t6*t7*tyi*2.0;
t136 = Ywb*t6*t7*tyi*2.0;
t137 = Zwa*t6*t8*tzi*2.0;
t138 = Zwb*t6*t8*tzi*2.0;
t139 = Xwa*t10*t6*txi*2.0;
t140 = Xwb*t6*t9*txi*2.0;
t141 = Xwa*t11*t6*txi*2.0;
t142 = Xwb*t12*t6*txi*2.0;
t143 = Ywb*t6*t9*tyi*2.0;
t144 = Ywa*t10*t6*tyi*2.0;
t145 = Xwa*t13*t6*txi*2.0;
t146 = Xwb*t14*t6*txi*2.0;
t147 = Ywa*t11*t6*tyi*2.0;
t148 = Ywb*t12*t6*tyi*2.0;
t149 = Zwb*t6*t9*tzi*2.0;
t150 = Zwa*t10*t6*tzi*2.0;
t151 = Ywa*t13*t6*tyi*2.0;
t152 = Ywb*t14*t6*tyi*2.0;
t153 = Zwb*t12*t6*tzi*2.0;
t154 = Zwa*t11*t6*tzi*2.0;
t155 = Zwa*t13*t6*tzi*2.0;
t156 = Zwb*t14*t6*tzi*2.0;
t157 = Xwa*t6*t7*txi*2.0;
t158 = Xwb*t6*t7*txi*2.0;
t159 = Xwa*t6*t8*txi*2.0;
t160 = Xwb*t6*t8*txi*2.0;
t161 = Ywa*t5*t6*tyi*2.0;
t162 = Ywb*t5*t6*tyi*2.0;
t163 = Ywa*t6*t8*tyi*2.0;
t164 = Ywb*t6*t8*tyi*2.0;
t165 = Zwa*t5*t6*tzi*2.0;
t166 = Zwb*t5*t6*tzi*2.0;
t167 = Zwa*t6*t7*tzi*2.0;
t168 = Zwb*t6*t7*tzi*2.0;
t72 = -t133-t134-t135-t136-t137-t138-t139-t140-t141-t142-t143-t144-t145-t146-t147-t148-t149-t150-t151-t152-t153-t154-t155-t156-t157-t158-t159-t160-t161-t162-t163-t164-t165-t166-t167-t168+t28+t30+t32+t33+t34+t35+t36+t37+t38+t39+t40+t41+t42+t43+t44+t45+t46+t47+t48+t49+t50+t51+t52+t53+t54+t55+t56+t57+t58+t59+t60+t61+t62+t63+t64+t65+t66+t67+t68+t69+t70+t71;
t93 = pow(t3,2);
t94 = t4*t9*2.025e5;
t95 = t10*t4*2.025e5;
t96 = t12*t4*2.025e5;
t97 = t11*t4*2.025e5;
t98 = t14*t4*2.025e5;
t99 = t13*t4*2.025e5;
t100 = t4*t5*4.05e5;
t101 = t4*t7*4.05e5;
t102 = t4*t8*4.05e5;
t127 = Xwa*t4*txi*4.05e5;
t128 = Xwb*t4*txi*4.05e5;
t129 = Ywa*t4*tyi*4.05e5;
t130 = Ywb*t4*tyi*4.05e5;
t131 = Zwa*t4*tzi*4.05e5;
t132 = Zwb*t4*tzi*4.05e5;
t103 = t100+t101+t102-t127-t128-t129-t130-t131-t132+t94+t95+t96+t97+t98+t99;
t104 = 1/t103;
t105 = pow(t26,2);
t106 = t4*t9*8.1e5;
t107 = t10*t4*8.1e5;
t108 = t12*t4*8.1e5;
t109 = t11*t4*8.1e5;
t110 = t14*t4*8.1e5;
t111 = t13*t4*8.1e5;
t112 = t4*t5*1.62e6;
t113 = t4*t7*1.62e6;
t114 = t4*t8*1.62e6;
t170 = Xwa*t4*txi*1.62e6;
t171 = Xwb*t4*txi*1.62e6;
t172 = Ywa*t4*tyi*1.62e6;
t173 = Ywb*t4*tyi*1.62e6;
t174 = Zwa*t4*tzi*1.62e6;
t175 = Zwb*t4*tzi*1.62e6;
t115 = t106+t107+t108+t109+t110+t111+t112+t113+t114-t170-t171-t172-t173-t174-t175;
t116 = 1/t115;
/* */
/*  Xa */
t117 = Xwa*t10*t6*2.0;
t118 = Xwa*t11*t6*2.0;
t119 = Xwa*t13*t6*2.0;
t120 = Xwa*t5*t6*2.0;
t121 = Xwb*t5*t6*4.0;
t122 = Xwa*t6*t7*2.0;
t123 = Xwa*t6*t8*2.0;
t124 = Ywb*t6*txi*tyi*4.0;
t125 = Zwb*t6*txi*tzi*4.0;
t126 = t117+t118+t119+t120+t121+t122+t123+t124+t125-t10*t6*txi*2.0-t11*t6*txi*2.0-t13*t6*txi*2.0-t5*t6*txi*2.0-t6*t7*txi*2.0-t6*t8*txi*2.0-Xwa*Xwb*t6*txi*4.0-Xwa*Ywb*t6*tyi*4.0-Xwa*Zwb*t6*tzi*4.0;
t169 = 1/pow(t18,2);
t176 = exp(t26);
t177 = t105*t116*t72*(-1.0/2.0)-t104*t72*t93*(1.0/2.0);
t178 = exp(t177);
t179 = focal*mR21*t19*t20;
t180 = t179-focal*mR31*t169*t20*t76;
t181 = t176+1.0;
Xa = -(t178*(t105*t116*t126*(1.0/2.0)+t104*t126*t93*(1.0/2.0)+t116*t180*t26*t72+t104*t3*t72*(Sx*focal*mR11*t19*t2-Sx*focal*mR31*t169*t2*t86)-1/pow(t103,2)*t72*t93*(Xwa*t4*4.05e5-t4*txi*4.05e5)*(1.0/2.0)-t105*1/pow(t115,2)*t72*(Xwa*t4*1.62e6-t4*txi*1.62e6)*(1.0/2.0)))/t181-t176*t178*t180*1/pow(t181,2);

/*  Ya */
t117 = Ywa*t10*t6*2.0;
t118 = Ywa*t11*t6*2.0;
t119 = Ywa*t13*t6*2.0;
t120 = Ywa*t5*t6*2.0;
t121 = Ywa*t6*t7*2.0;
t122 = Ywb*t6*t7*4.0;
t123 = Ywa*t6*t8*2.0;
t124 = Xwb*t6*txi*tyi*4.0;
t125 = Zwb*t6*tyi*tzi*4.0;
t126 = t117+t118+t119+t120+t121+t122+t123+t124+t125-t10*t6*tyi*2.0-t11*t6*tyi*2.0-t13*t6*tyi*2.0-t5*t6*tyi*2.0-t6*t7*tyi*2.0-t6*t8*tyi*2.0-Xwb*Ywa*t6*txi*4.0-Ywa*Ywb*t6*tyi*4.0-Ywa*Zwb*t6*tzi*4.0;
t169 = 1/pow(t18,2);
t176 = exp(t26);
t177 = t105*t116*t72*(-1.0/2.0)-t104*t72*t93*(1.0/2.0);
t178 = exp(t177);
t179 = focal*mR22*t19*t20;
t180 = t179-focal*mR32*t169*t20*t76;
t181 = t176+1.0;
Ya = -(t178*(t105*t116*t126*(1.0/2.0)+t104*t126*t93*(1.0/2.0)+t116*t180*t26*t72+t104*t3*t72*(Sx*focal*mR12*t19*t2-Sx*focal*mR32*t169*t2*t86)-1/pow(t103,2)*t72*t93*(Ywa*t4*4.05e5-t4*tyi*4.05e5)*(1.0/2.0)-t105*1/pow(t115,2)*t72*(Ywa*t4*1.62e6-t4*tyi*1.62e6)*(1.0/2.0)))/t181-t176*t178*t180*1/pow(t181,2);

/*  Xb */
t117 = Xwb*t6*t9*2.0;
t118 = Xwb*t12*t6*2.0;
t119 = Xwb*t14*t6*2.0;
t120 = Xwa*t5*t6*4.0;
t121 = Xwb*t5*t6*2.0;
t122 = Xwb*t6*t7*2.0;
t123 = Xwb*t6*t8*2.0;
t124 = Ywa*t6*txi*tyi*4.0;
t125 = Zwa*t6*txi*tzi*4.0;
t126 = t117+t118+t119+t120+t121+t122+t123+t124+t125-t12*t6*txi*2.0-t14*t6*txi*2.0-t5*t6*txi*2.0-t6*t7*txi*2.0-t6*t8*txi*2.0-t6*t9*txi*2.0-Xwa*Xwb*t6*txi*4.0-Xwb*Ywa*t6*tyi*4.0-Xwb*Zwa*t6*tzi*4.0;
t169 = 1/pow(t24,2);
t176 = exp(t26);
t177 = t105*t116*t72*(-1.0/2.0)-t104*t72*t93*(1.0/2.0);
t178 = exp(t177);
t179 = focal*mR21*t20*t25;
t180 = t179-focal*mR31*t169*t20*t81;
t181 = t176+1.0;
Xb = (t178*(t105*t116*t126*(-1.0/2.0)-t104*t126*t93*(1.0/2.0)+t116*t180*t26*t72+t104*t3*t72*(Sx*focal*mR11*t2*t25-Sx*focal*mR31*t169*t2*t91)+1/pow(t103,2)*t72*t93*(Xwb*t4*4.05e5-t4*txi*4.05e5)*(1.0/2.0)+t105*1/pow(t115,2)*t72*(Xwb*t4*1.62e6-t4*txi*1.62e6)*(1.0/2.0)))/t181+t176*t178*t180*1/pow(t181,2);

/* Yb */
t117 = Ywb*t6*t9*2.0;
t118 = Ywb*t12*t6*2.0;
t119 = Ywb*t14*t6*2.0;
t120 = Ywb*t5*t6*2.0;
t121 = Ywa*t6*t7*4.0;
t122 = Ywb*t6*t7*2.0;
t123 = Ywb*t6*t8*2.0;
t124 = Xwa*t6*txi*tyi*4.0;
t125 = Zwa*t6*tyi*tzi*4.0;
t126 = t117+t118+t119+t120+t121+t122+t123+t124+t125-t12*t6*tyi*2.0-t14*t6*tyi*2.0-t5*t6*tyi*2.0-t6*t7*tyi*2.0-t6*t8*tyi*2.0-t6*t9*tyi*2.0-Xwa*Ywb*t6*txi*4.0-Ywa*Ywb*t6*tyi*4.0-Ywb*Zwa*t6*tzi*4.0;
t169 = 1/pow(t24,2);
t176 = exp(t26);
t177 = t105*t116*t72*(-1.0/2.0)-t104*t72*t93*(1.0/2.0);
t178 = exp(t177);
t179 = focal*mR22*t20*t25;
t180 = t179-focal*mR32*t169*t20*t81;
t181 = t176+1.0;
Yb = (t178*(t105*t116*t126*(-1.0/2.0)-t104*t126*t93*(1.0/2.0)+t116*t180*t26*t72+t104*t3*t72*(Sx*focal*mR12*t2*t25-Sx*focal*mR32*t169*t2*t91)+1/pow(t103,2)*t72*t93*(Ywb*t4*4.05e5-t4*tyi*4.05e5)*(1.0/2.0)+t105*1/pow(t115,2)*t72*(Ywb*t4*1.62e6-t4*tyi*1.62e6)*(1.0/2.0)))/t181+t176*t178*t180*1/pow(t181,2);

    plhs[0]=Xa;
    plhs[1]=Ya;
    plhs[2]=Xb;
    plhs[3]=Yb;


}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {


    /* /Declarations */
    const mxArray *Dpxdata, *Dpydata, *Sxdata, *focaldata, \
            *mR11data, *mR12data, *mR13data, *mR21data, *mR22data, *mR23data, *mR31data, *mR32data, *mR33data, \
            *txdata, *tydata, *tzdata, *Xmtx, *Ymtx, *Ximtx, *Yimtx, *rhodata, *txidata, *tyidata, *tzidata;

    double Dpx, Dpy, Sx, Xwa, Xwb, Ywa, Ywb, Zwa, Zwb, focal, \
            mR11, mR12, mR13, mR21, mR22, mR23, mR31, mR32, mR33, \
            tx, ty, tz, *X, *Y, *Xi, *Yi, rho, txi, tyi, tzi;

		double Za, Zb;
		double sizeOnImagea, sizeOnImageb;
		double tmpdouble;
		double tmp, tmp2, tmp2x, tmp2y, tmpx, tmpy;
		double sigfac;


			int i,j,k,t, ind, ind2, ind3, ind4,icnt,jcnt;
			double B[4], C[4];
double *v;
double *tvis;
double *occx;
double *occy;
double *ddvix;
double *ddviy;

	unsigned int N;
	unsigned int F;
	int dims[3];

			double camPar[25];
				double dX[4];
				double *sf;

			double *occmat;
			double *occmatdX;
			double *occmatdY;
			double *occmatdX2;
			double *occmatdY2;
			double *ol;

    /* //Copy input pointer x */
	Xmtx = prhs[0];	Ymtx = prhs[1];	Ximtx = prhs[2];	Yimtx = prhs[3];
    Dpxdata = prhs[4];    Dpydata = prhs[5];    Sxdata = prhs[6];

    focaldata = prhs[7];
    mR11data = prhs[8];    mR12data = prhs[9];    mR13data = prhs[10];
    mR21data = prhs[11];    mR22data = prhs[12];    mR23data = prhs[13];
    mR31data = prhs[14];    mR32data = prhs[15];    mR33data = prhs[16];
    txdata = prhs[17];    tydata = prhs[18];    tzdata = prhs[19];
	txidata = prhs[20];    tyidata = prhs[21];    tzidata = prhs[22];
	rhodata = prhs[23];


    Dpx = (double)(mxGetScalar(Dpxdata));    Dpy = (double)(mxGetScalar(Dpydata));
    Sx = (double)(mxGetScalar(Sxdata));
    /*
//     Xwa = (double)(mxGetScalar(Xwadata));    Xwb = (double)(mxGetScalar(Xwbdata));
//     Ywa = (double)(mxGetScalar(Ywadata));    Ywb = (double)(mxGetScalar(Ywbdata));
//     Zwa = (double)(mxGetScalar(Zwadata));    Zwb = (double)(mxGetScalar(Zwbdata));
     */
    focal = (double)(mxGetScalar(focaldata));
    mR11 = (double)(mxGetScalar(mR11data));    mR12 = (double)(mxGetScalar(mR12data));    mR13 = (double)(mxGetScalar(mR13data));
    mR21 = (double)(mxGetScalar(mR21data));    mR22 = (double)(mxGetScalar(mR22data));    mR23 = (double)(mxGetScalar(mR23data));
    mR31 = (double)(mxGetScalar(mR31data));    mR32 = (double)(mxGetScalar(mR32data));    mR33 = (double)(mxGetScalar(mR33data));
    tx = (double)(mxGetScalar(txdata));    ty = (double)(mxGetScalar(tydata));    tz = (double)(mxGetScalar(tzdata));
	txi = (double)(mxGetScalar(txidata));    tyi = (double)(mxGetScalar(tyidata));    tzi = (double)(mxGetScalar(tzidata));
	X = mxGetPr(Xmtx);	Y = mxGetPr(Ymtx);	Xi = mxGetPr(Ximtx);	Yi = mxGetPr(Yimtx); rho=(double)(mxGetScalar(rhodata));





	N = mxGetN(Xmtx);
	F = mxGetM(Xmtx);
/*// 	mexPrintf("%i %i\n",F,N);*/


	dims[0]=F;dims[1]=N;dims[2]=N;

	plhs[0] = mxCreateDoubleMatrix(F, N, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(F, N, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(F, N, mxREAL);
	plhs[3] = mxCreateDoubleMatrix(F, N, mxREAL);
	plhs[4] = mxCreateNumericArray(3, dims,mxDOUBLE_CLASS,   mxREAL);
	plhs[5] = mxCreateNumericArray(3, dims,mxDOUBLE_CLASS,   mxREAL);

	v = mxGetPr(plhs[0]);
	tvis = mxGetPr(plhs[1]);
	occx = mxGetPr(plhs[2]);
	occy = mxGetPr(plhs[3]);
	ddvix = mxGetPr(plhs[4]);
	ddviy = mxGetPr(plhs[5]);

	icnt=0;
/*
// 	for (t=0; t<F; t++) {
// 		for (i=0; i<N; i++) {
// 			for (j=0; j<N; j++) {
// 				ddvix[icnt++]=icnt;
// // 				mexPrintf("%f %i\n",ddvix[icnt],icnt++);
//
// 			}
// 		}
// 	}
//
// icnt=0;
// 	for (t=0; t<F; t++) {
// 		for (i=0; i<N; i++) {
// 			for (j=0; j<N; j++) {
// 				mexPrintf("%i %i %i %f\n",t,i,j,ddvix[icnt++]);
//
// 			}
// 		}
// 	}
*/

									camPar[0]=Dpx;
									camPar[1]=Dpy;
									camPar[2]=Sx;
									camPar[3]=Xwa;
									camPar[4]=Xwb;
									camPar[5]=Ywa;
									camPar[6]=Ywb;
									camPar[7]=Zwa;
									camPar[8]=Zwb;
									camPar[9]=focal;
											camPar[10]=mR11;
											camPar[11]=mR12;
											camPar[12]=mR13;
											camPar[13]=mR21;
											camPar[14]=mR22;
											camPar[15]=mR23;
											camPar[16]=mR31;
											camPar[17]=mR32;
											camPar[18]=mR33;
											camPar[19]=tx;
											camPar[20]=ty;
											camPar[21]=tz;
											camPar[22]=txi;
											camPar[23]=tyi;
											camPar[24]=tzi;
				dX[0]=dX[1]=dX[2]=dX[3]=0;
				sf = (double*)malloc(N*N * sizeof (double));

			occmat = (double*)malloc(N*N * sizeof (double));
			occmatdX = (double*)malloc(N*N * sizeof (double));
			occmatdY = (double*)malloc(N*N * sizeof (double));
			occmatdX2 = (double*)malloc(N*N * sizeof (double));
			occmatdY2 = (double*)malloc(N*N * sizeof (double));
			ol = (double*)malloc(N*N * sizeof (double));


		for (t=0; t<F; t++) {
			for (i=0; i<N*N; i++) {
				occmat[i]=0.;occmatdX[i]=0.;occmatdY[i]=0.;
				occmatdX2[i]=0.; occmatdY2[i]=0.; ol[i]=0.;
			}

			icnt=0;

			for (i=0; i<N; i++) {
				ind=i*F+t;
				Xwa=X[ind];
				if (Xwa != 0) {
					double a[2];
					double A[2];
					icnt++;
					Ywa=Y[ind];
					Zwa=900; Zwb=900;

/*// 					Za=sqrt((tx-Xwa)*(tx-Xwa) + (ty-Ywa)*(ty-Ywa) +(tz-Zwa)*(tz-Zwa));*/
					Za=sqrt((txi-Xwa)*(txi-Xwa) + (tyi-Ywa)*(tyi-Ywa) +(tzi-Zwa)*(tzi-Zwa));
					sizeOnImagea=1800*focal/Za/Dpy;
					tmpdouble=(sizeOnImagea/2)*(sizeOnImagea/2);
					A[0]=tmpdouble/4;
					A[1]=tmpdouble;


					for (j=0; j<N; j++) {
						ind2=j*F+t;
						Xwb=X[ind2];

						if (Xwb != 0) {
							if (j>i) {
								double b[2];
								double B[2], C[2];
								double alessb[2];
								ind3=i*N+j;
								Ywb=Y[ind2];

/*// 								Zb=sqrt((tx-Xwb)*(tx-Xwb) + (ty-Ywb)*(ty-Ywb) +(tz-Zwb)*(tz-Zwb));*/
								Zb=sqrt((txi-Xwb)*(txi-Xwb) + (tyi-Ywb)*(tyi-Ywb) +(tzi-Zwb)*(tzi-Zwb));
								sizeOnImageb=1800*focal/Zb/Dpy;
								tmpdouble=(sizeOnImageb/2)*(sizeOnImageb/2);
								B[0]=tmpdouble/4; B[1]=tmpdouble;
								C[0]=1/(A[0]+B[0]); C[1]=1/(A[1]+B[1]);

								a[0]=Xi[ind]; a[1]=Yi[ind];
								b[0]=Xi[ind2]; b[1]=Yi[ind2];
								alessb[0]=a[0]-b[0]; alessb[1]=a[1]-b[1];
								tmpdouble=alessb[0]*alessb[0]*C[0]+alessb[1]*alessb[1]*C[1];

								ol[ind3]=exp(-.5 * tmpdouble);
/*// 								mexPrintf("%i %i %f %f %f %f %f\n",i,j,alessb[0],alessb[1],C[0],C[1],ol[i][j]);*/
							}
						}
					}


					for (j=0; j<N; j++) {
						ind2=j*F+t;
						Xwb=X[ind2];
						if (Xwb != 0) {
							ind3=i*N+j;
							if (j<i) {
								ol[ind3]=ol[j*N+i];
							}

							if (i!=j) {

								Ywb=Y[ind2];
								sigfac=1./(1+exp(Yi[ind]-Yi[ind2]));
								sf[ind3]=sigfac;
								occmat[ind3]=ol[ind3]*sigfac;

								if (ol[ind3] > 0.0001){

									camPar[3]=Xwa;camPar[4]=Xwb;
									camPar[5]=Ywa;camPar[6]=Ywb;
									camPar[7]=Zwa;camPar[8]=Zwb;
									dX[0]=0;dX[1]=0;dX[2]=0;dX[3]=0;

									dfnddALL(camPar,dX);

									occmatdX[ind3]=dX[0];
									occmatdY[ind3]=dX[1];
									occmatdX2[ind3]=dX[2];
									occmatdY2[ind3]=dX[3];

								}

							}
						}
					}
				}
			}
            /*
// 			if (t==0) {
// 				for (k=0; k<N; k++) {
// 					for (j=0; j<N; j++) {
// 						if (isnan(occmatdX[k][j])) {occmatdX[k][j]=0;}
// 						if (isnan(occmatdY[k][j])) {occmatdY[k][j]=0;}
// 						if (isnan(occmatdX2[k][j])) {occmatdX2[k][j]=0;}
// 						if (isnan(occmatdY2[k][j])) {occmatdY2[k][j]=0;}
// 					}
// // 					mexPrintf("\n");
// 				}
// // 				mexPrintf("\n");
// 			}
			ind3=0;
			if (t==0) {
				for (k=0; k<N; k++) {
					for (j=0; j<N; j++) {
						mexPrintf("%.15f  ",occmatdX[ind3++]*1000);
					}
					mexPrintf("\n");
				}
				mexPrintf("\n");
			}		
			if (t==0) {
				for (k=0; k<N; k++) {
					for (j=0; j<N; j++) {
						mexPrintf("%.15f  ",occmatdX2[k][j]);
					}
					mexPrintf("\n");
				}
				mexPrintf("\n");
			}		
//    if nargout > 1
            */
			for (i=0; i<N; i++) {
				ind=i*F+t;
				Xwa=X[ind];
				if (Xwa != 0) {

					double sumdx=0;
					double sumdy=0;
					for (j=0; j<N; j++) {
						ind2=j*F+t;

						Xwb=X[ind2];
						if (Xwb != 0) {

							ind4=j*F*N+i*F+t;

							if (j!=i) {
								tmp=0;
								for (k=0; k<N; k++){
									ind3=j*N+k;
									tmp+=occmat[ind3];
								}
								ind3=j*N+i;
								tmpx=-rho*exp(-rho*tmp)*occmatdX2[ind3];
								sumdx+= tmpx;
								ddvix[ind4] = tmpx;

								tmpy=-rho*exp(-rho*tmp)*occmatdY2[ind3];
								sumdy+= tmpy;
								ddviy[ind4] = tmpy;
							}
							else {
								tmp=0; tmp2x=0; tmp2y=0;
								for (k=0; k<N; k++){
									ind3=i*N+k;
									tmp+=occmat[ind3];
									tmp2x+=occmatdX[ind3];
									tmp2y+=occmatdY[ind3];
	/*									if (t==0)	{
											mexPrintf("%i %i %f %f %f\n",i+1,k+1,tmp,tmpx,occmatdX[i][k]);
										}	*/

								}
								tmpx=-rho*exp(-rho*tmp)*tmp2x;
	/*							if (t==0)	{
									mexPrintf("%i %i %f %f %f\n",t+1,i+1,tmp,tmpx,tmp2x);
								}							*/

								sumdx+= tmpx;
								ddvix[ind4] = tmpx;

								tmpy=-rho*exp(-rho*tmp)*tmp2y;
								sumdy+= tmpy;
								ddviy[ind4] = tmpy;
							}
						}
					}
/*				if (t==0)	{
					mexPrintf("%i %i %f %f\n",t+1,i+1,sumdx,sumdy);
				}*/
				occx[ind]=sumdx;
				occy[ind]=sumdy;
				}
			}

		for (i=0; i<N; i++) {
			ind=i*F+t;
			tmp=0;
			for (j=0; j<N; j++) {
				ind3=i*N+j;
				tmp+=occmat[ind3];
			}
			v[ind]=exp(-rho*tmp);
			tvis[ind]=exp(-tmp);
		}
	}

	icnt=0;
	jcnt=0;
    /*
// 	for (t=0; t<F; t++) {
// 				for (i=0; i<N; i++) {
// 					if (isnan(v[icnt])) {mexPrintf("%i %i v\n",t,i); }
// 					if (isnan(tvis[icnt])) {mexPrintf("%i %i tv\n",t,i); }
// 					if (isnan(occx[icnt])) {mexPrintf("%i %i ox\n",t,i); }
// 					if (isnan(occy[icnt])) {mexPrintf("%i %i oy\n",t,i); }
// 					for (j=0; j<N; j++) {
// 						if (isnan(ddvix[jcnt])) {mexPrintf("%i %i ddx\n",t,i); }
// 						if (isnan(ddviy[jcnt])) {mexPrintf("%i %i ddy\n",t,i); }
// 						jcnt++;
// 					}
// 					icnt++;
// // 					mexPrintf("\n");
// 				}
// 	}
				for (i=0; i<N; i++) {
					for (j=0; j<N; j++) {
						mexPrintf("%.15f ",sf[i][j]);
					}
					mexPrintf("\n");
				}*/
			free(sf);
			free(occmat);
			free(occmatdX);
			free(occmatdY);
			free(occmatdX2);
			free(occmatdY2);
			free(ol);
			}
