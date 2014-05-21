#include "mex.h"
#include <assert.h>
#include <math.h>
#include <string.h>

#define HUGE 1000000

const short kInvalidIndex = -1;

/*squared 2d norm (so in the end its just x^2+y^2, which is faster than also taking the squareroot and enough for comparisons)*/
double norm2dSqr(double x, double y) {
	return x * x + y * y;
}

/*2d norm*/
double norm2d(double x, double y) {
	return sqrt(x * x + y * y);
}

/*counts the amount of nonzero elements in a*/
int countNonZero(int length, const double* const a) {
	int i, nz;
	nz = 0;

	for (i = 0; i < length; i++) {
		if (a[i] != 0) {
			nz++;
		}
	}
	return nz;

}

/*finds all values in a != valueToIgnore and returns there indices in res*/
void find(int length, const short * const a, int* const foundCount, short * const res, short valueToIgnore){
	int i;
	*foundCount = 0;

	for (i = 0; i < length; i++) {
		if (a[i] != valueToIgnore) {
			res[(*foundCount)++] = i;
		}
	}
}

/*finds all values in a != 0 and returns there indices in res*/
void findD(int length, const double * const a, int* const foundCount, short * const res){
	int i;
	*foundCount = 0;

	for (i = 0; i < length; i++) {
		if (a[i] != 0) {
			res[(*foundCount)++] = i;
		}
	}
}

/*executes find(~M(t,:) & Xgt(t,:)); in c*/
void findNotMbutXgt(double* const Xgt, short* const M, int Ngt, short* const GTsNotMapped, int* const GTsNotMappedCount){
	int i;
	*GTsNotMappedCount = 0;
	for(i = 0; i < Ngt; i++){
		if((M[i] == kInvalidIndex) && Xgt[i])
			GTsNotMapped[(*GTsNotMappedCount)++] = i;
	}
}

/*executes setdiff(find(X(t,:)),M(t,:)); in c (well the "(t," needs to be already factored into X/M)*/
void setDiffFind(double* const X, int N, short* const M, int Ngt, short* const EsNotMapped, int* const EsNotMappedCount){
	int i;
	char* temp = (char*) calloc(N, sizeof(char));
	*EsNotMappedCount = 0;
	/*first get a boolean map of X, true == x at i was in M, otherwise false*/
	for(i = 0; i < Ngt; i++){
		if(M[i] != kInvalidIndex)
			temp[M[i]] = true;
	}
	/*now we only need to check if there is a X[i] for a false temp[i]
	well this isn't really correct because the matlab code would also find stuff that is in M but not in X, however I just assume that never happens :)*/
	*EsNotMappedCount = 0;
	for(i = 0; i < N; i++){
		if(X[i] && !temp[i]){
			EsNotMapped[(*EsNotMappedCount)++] = i;
		}
	}
	free(temp);
}

double sum(short* a, int len){
	double sum = 0.0;
	int i;
	for ( i = 0; i < len; i++){
		sum += a[i];
	}
	return sum;
}

double sumD(double* a, int len){
	double sum = 0.0;
	int i;
	for ( i = 0; i < len; i++){
		sum += a[i];
	}
	return sum;
}

void setDiff(const short* const a, const int aLength, const short* const b, const int bLength, short* const out, int* const outLength){
	int i,j;
	bool contains;
	*outLength = 0;
	for(i = 0; i < aLength; i++){
		contains = false;
		for(j = 0; j < bLength; j++){
			if(a[i] == b[j]) {
				contains = true;
				break;
			}
		}
		if(!contains)
			out[(*outLength)++] = a[i];
	}
}

/*intersect for two sorted lists*/
void intersectSorted(const short* const a, const int aLength, const short* const b, const int bLength, short* const out, int* const outLength){
	int i,j;
	i = j = 0;
	*outLength = 0;
	while(i < aLength){
		if(a[i] == b[j]) {
			out[(*outLength)++] = a[i];
			/*only increase if we found it in b*/
			i++;
		}
		else {
			/*only increase j if we couldn't find it at the current b*/
			j++;
		}
	}
}

/*intersect for two lists */
void intersect(const short* const a, const int aLength, const short* const b, const int bLength, short* const out, int* const outLength){
	int i,j;
	*outLength = 0;
	for(i = 0; i < aLength; i++){
		for(j = 0; j < bLength; j++){
			if(a[i] == b[j]) {
				out[(*outLength)++] = a[i];
				break;
			}
		}
	}
}

/*
computes stuff
input:
0 Xgt
1 Ygt
2 X
3 Y
4 td

output:
0 MOTA
1 MOTP
2 ma
3 fpa
4 mmea
5 idsw
6 missed
7 falsepositives
8 idswitches
9 alltracked
10 allfalsepos
11 MT
12 PT
13 ML
14 recall
15 precision
16 fafrm
17 FRA
18 MOTAL
19 alld
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    /* /Declarations */ 
	/* in vars*/
    const mxArray *XgtMtx, *YgtMtx, *XMtx, *YMtx, *tdata;
    double *Xgt, *Ygt, *X, *Y, td;
    int F, N, Fgt, Ngt;
	
	/* out vars*/
	double *MOTA, *MOTP, *ma, *fpa, *mmea, *idsw, *missed, *falsepositives, *idswitches, *alltracked, *allfalsepos, *MT, *PT, *ML, *recall, *precision, *fafrm, *FRA, *MOTAL, *alld;

	/* temp*/
	int i, t, foundCount, map, GTsNotMappedCount, EsNotMappedCount, j, lastNotEmpty;
	int minDistIndexGT, minDistIndexE,trlengtha;
	short *mme, *allCurTrackedCount;
	int tNgtOffset, tMinus1NgtOffset, tNOffset, NgtCurMapOffset, NLastMapOffset, XFoundCount, gtOffset, allTrackersCount, mappedTrackersCount, curMFoundCount, tempCount, curTrackedCount, gtFramesCount;
    short *g, *GTsNotMapped, *EsNotMapped, *curTracked, *M, *mappings, *XFound, *mappedTrackers, *curMFound, *temp, *fp, *m, *gtFrames;
	short eid;
	const short *allTrackers;
	double mindist, dist, res;
	double *alldist, *allTrackedDist;

    /*Copy input pointer x */
    XgtMtx = prhs[0];
	YgtMtx = prhs[1];
    XMtx = prhs[2];
	YMtx = prhs[3];
	tdata = prhs[4];

    /*Get input matrices*/
	Xgt = mxGetPr(XgtMtx);
	Ygt = mxGetPr(YgtMtx);
	/*F and N are switched because we expect the matrices transposed*/
    Fgt = mxGetN(XgtMtx);
    Ngt = mxGetM(XgtMtx);

    X = mxGetPr(XMtx);
	Y = mxGetPr(YMtx);
    F = mxGetN(XMtx);
    N = mxGetM(XMtx);

    assert(F==Fgt);

    td = mxGetScalar(tdata);

	/*init scalar outputs (everything but alld, allfalsepos, alltracked)*/
	for(i = 0; i <= 8; i++) {
		plhs[i] = mxCreateDoubleMatrix(1, 1, mxREAL);
	}
	for(i = 11; i <= 18; i++) {
		plhs[i] = mxCreateDoubleMatrix(1, 1, mxREAL);
	}

	/*init multi dim output*/
	plhs[9] = mxCreateDoubleMatrix(Ngt, F, mxREAL);
	plhs[10] = mxCreateDoubleMatrix(N, F, mxREAL);
	plhs[19] = mxCreateDoubleMatrix(Ngt, F, mxREAL);

	/*get output pointer*/
	MOTA = mxGetPr(plhs[0]);
	MOTP = mxGetPr(plhs[1]);
	ma = mxGetPr(plhs[2]);
	fpa = mxGetPr(plhs[3]);
	mmea = mxGetPr(plhs[4]);
	idsw = mxGetPr(plhs[5]);
	missed = mxGetPr(plhs[6]);
	falsepositives = mxGetPr(plhs[7]);
	idswitches = mxGetPr(plhs[8]);
	alltracked = mxGetPr(plhs[9]);
	allfalsepos = mxGetPr(plhs[10]);
	MT = mxGetPr(plhs[11]);
	PT = mxGetPr(plhs[12]);
	ML = mxGetPr(plhs[13]);
	recall = mxGetPr(plhs[14]);
	precision = mxGetPr(plhs[15]);
	fafrm = mxGetPr(plhs[16]);
	FRA = mxGetPr(plhs[17]);
	MOTAL = mxGetPr(plhs[18]);
	alld = mxGetPr(plhs[19]);

	/*init temp vars*/
	M = (short*) malloc(sizeof(short) * F * Ngt);
	mappings = (short*) malloc(sizeof(short) * Ngt);
	g = (short*) malloc(sizeof(short) * F);
	GTsNotMapped = (short*) malloc(sizeof(short) * Ngt);
	EsNotMapped = (short*) malloc(sizeof(short) * Ngt);
	XFound = (short*) malloc(sizeof(short) * N);
	curMFound = (short*) malloc(sizeof(short) * N);
	temp = (short*) malloc(sizeof(short) * N);
	curTracked = (short*) malloc(sizeof(short) * Ngt);
	mme = (short*) malloc(sizeof(int) * F);
	alldist = (double*) malloc(sizeof(double) * Ngt * N);
	allCurTrackedCount = (short*) malloc(sizeof(short) * F);
	allTrackedDist = (double*) malloc(sizeof(double) *F * Ngt);
	fp = (short*) malloc(sizeof(short) * F);
	m = (short*) malloc(sizeof(short) * F);
	gtFrames = (short*) malloc(sizeof(short) * F);
	mappedTrackers = (short*) malloc(sizeof(short) * Ngt);

	/*init M to -1s (because 0s are people(aka indices) too)*/
	for(i = 0; i < F * Ngt; i++){
		M[i] = kInvalidIndex;
	}

	for(i = 0; i < F * N; i++){
		allfalsepos[i] = 0;
	}

	for(i = 0; i < F* Ngt; i++) {
		allTrackedDist[i] = 0.0;
	}

	/*now to the work part*/
	for(t = 0; t < F; t++){

		/*how many gt targets are there?
		also compute the offset for arrays of F*Ngt for the current t*/
		tNgtOffset = Ngt * t;
		g[t] = countNonZero(Ngt, Xgt + tNgtOffset);

		/*some indices precomputation because we will need it later anyway*/
		tNOffset = N * t;
		tMinus1NgtOffset = Ngt * (t - 1);

		if(t > 0){
			/*do some stuff
			 find all non 0 M values of the last iteration and store them in mappings*/
			find(Ngt, M + tMinus1NgtOffset, &foundCount, mappings, kInvalidIndex);
			
			/*now take the last mappings if they are close enough to the ground truth*/
			for(i = 0; i < foundCount; i++){
				map = mappings[i];

				/*get offsets for the current time plus the current map and
				the current timestep plus the last M*/
				NgtCurMapOffset = tNgtOffset + map;
				NLastMapOffset = tNOffset + M[tMinus1NgtOffset + map];

				/*if both X[bla] and Xgt[bla] are nonZero
				and they are close enough to each other, set the old map as the new one*/
				if( Xgt[NgtCurMapOffset] 
				    && X[NLastMapOffset]
				    && norm2dSqr(Xgt[NgtCurMapOffset] - X[NLastMapOffset],
				              Ygt[NgtCurMapOffset] - Y[NLastMapOffset]) <= td*td) {
					M[NgtCurMapOffset] = M[tMinus1NgtOffset + map];
				}
			}
		}

		/*setup alldist to infinity(== HUGE_VAL)*/
		for(i = 0; i < Ngt * N; i++){
			alldist[i] = HUGE_VAL;
		}

		/*find all gt for which there is no entry in M*/
		findNotMbutXgt(Xgt + tNgtOffset, M + tNgtOffset, Ngt, GTsNotMapped, &GTsNotMappedCount);
		
		/*find all X != 0 (TODO: change setdifffind so it uses this result)*/
		findD(N, X + tNOffset, &XFoundCount, XFound);

		/*find all x for which there is no entry in M*/
		setDiffFind(X + tNOffset, N, M + tNgtOffset, Ngt, EsNotMapped, &EsNotMappedCount);

		mindist = 0;
		minDistIndexGT = kInvalidIndex;
		minDistIndexE = kInvalidIndex;

		/*compute distances*/
		while (mindist < td && GTsNotMappedCount > 0 && EsNotMappedCount > 0) {
			/*compute distances for all not yet mapped GTs and Es and also find the closest x and gt*/
			mindist = HUGE;
			for (i = 0; i < GTsNotMappedCount; i++){ /*Ngt
				gtOffset = i * N;
				for(j = 0; j < EsNotMappedCount; j++){ /*N
					dist = norm2dSqr(Xgt[tNgtOffset + GTsNotMapped[i]] - X[tNOffset + EsNotMapped[j]],
					              Ygt[tNgtOffset +GTsNotMapped[i]] - Y[tNOffset + EsNotMapped[j]]);
					alldist[gtOffset + j] = dist;
					if(dist < mindist){
						mindist = dist;
						minDistIndexE = j;
						minDistIndexGT = i;
					}
				}
			}

			if (mindist <= td*td){
				/*add the mapping to M*/
				M[tNgtOffset + GTsNotMapped[minDistIndexGT]] = EsNotMapped[minDistIndexE];
				/*set all distances for the current gt to Inf*/
				gtOffset = minDistIndexGT * N;
				for(i = 0; i < EsNotMappedCount; i++){
					alldist[gtOffset + i] = HUGE;
				}
				/*remove the mapped gt and e from the notMapped arrays (should be changed later for now this is the simplest impl.)*/
				findNotMbutXgt(Xgt + tNgtOffset, M + tNgtOffset, Ngt, GTsNotMapped, &GTsNotMappedCount);
				setDiffFind(X + tNOffset, N, M + tNgtOffset, Ngt, EsNotMapped, &EsNotMappedCount);
			}
		}

		find(Ngt, M + tNgtOffset, &curTrackedCount, curTracked, kInvalidIndex);

		/*since alltrackers == Xfound we can just write it as follows*/
		allTrackers = XFound;
		allTrackersCount = XFoundCount;
		
		/*find mapped trackers
		TODO: replace find by assigning to curMfound already above, however for now this is good enough
		find currently assigned Xgts*/
		find(Ngt, M + tNgtOffset, &curMFoundCount, curMFound, kInvalidIndex);
		/*change each Xgt position to its respective value*/
		for(i = 0; i < curMFoundCount; i++){
			curMFound[i] = M[tNgtOffset + curMFound[i]];
		}
		intersect(curMFound, curMFoundCount, allTrackers, allTrackersCount, mappedTrackers, &mappedTrackersCount);
		
		/*find falsepositives*/
		setDiff(allTrackers, allTrackersCount, mappedTrackers, mappedTrackersCount, temp, &tempCount);
		/*convert them to double and store them in allfalsepos*/
		for(i = 0; i < tempCount; i++){
			/*temp+1 because -1 == invalid index in c however 0 == invalid index in matlab*/
			allfalsepos[tNOffset + i] = (double)(temp[i] + 1);
		}

		/*compute mismatch errors*/
		mme[t] = 0;
		if (t > 1){
			for(i = 0; i < curTrackedCount; i++){
				/*find last not empty ct*/
				lastNotEmpty = kInvalidIndex;
				j = t;
				while(j --> 0){
					if(M[j * Ngt + curTracked[i]] != kInvalidIndex) {
						lastNotEmpty = j;
						break;
					}
				}

				/*if we don't match increase mismatch*/
				if(Xgt[tMinus1NgtOffset + curTracked[i]] != 0
					&& lastNotEmpty != kInvalidIndex
					&& M[tNgtOffset + curTracked[i]] != M[lastNotEmpty * Ngt + curTracked[i]]){
						mme[t] += 1;
				}
			}
		}

		allCurTrackedCount[t] = curTrackedCount;

		/*compute distances for all currently tracked Xs*/
		for(i = 0; i < curTrackedCount; i++){
			eid = M[tNgtOffset + curTracked[i]];
			allTrackedDist[tNgtOffset + i] = norm2d(
				Xgt[tNgtOffset + curTracked[i]] - X[tNOffset + eid],
				Ygt[tNgtOffset + curTracked[i]] - Y[tNOffset + eid]);
		}

		fp[t] = XFoundCount - curTrackedCount;

		m[t] = g[t] - curTrackedCount;
	}

	/*copy alltrackedDist to alld*/
	for(i = 0; i < F * Ngt; i++){
		alld[i] = allTrackedDist[i];
	}

	/*copy M to alltracked*/
	for(i = 0; i < F * Ngt; i++){
		alltracked[i] = M[i] + 1;
	}

	*missed = sum(m, F);
	*falsepositives = sum(fp, F);
	*idswitches = sum(mme, F);

	/*compute average ditance to [0,1]*/
	*MOTP = 1 - sumD(allTrackedDist, F * Ngt) / sum(allCurTrackedCount, F) / td;

	*MOTAL = 1 - ( ( *missed + *falsepositives + log10(*idswitches + 1) ) / sum(g, F) );
	*MOTA = 1 - ( ( *missed + *falsepositives + *idswitches ) / sum(g, F) );
	*ma = *missed / sum(g, F);
	*fpa = *falsepositives / sum(g, F);
	*mmea = *idswitches / sum(g, F);
	*idsw = *idswitches / Ngt / Fgt;

	res = sum(allCurTrackedCount, F);
	*recall = sum(allCurTrackedCount, F) / sum(g, F);
	*precision = sum(allCurTrackedCount, F) / (*falsepositives + sum(allCurTrackedCount, F));
	*fafrm = *falsepositives / Fgt;

	/* MT PT ML*/
	*MT = *PT = *ML = 0;
	for(i = 0; i < Ngt; i++){
		gtFramesCount = 0;
		/*first get a bitmask of values in the column*/
		for(j = 0; j < F; j++) {
			gtFrames[j] = Xgt[j * Ngt + i] != 0;
		}
		/*Nowget their indices*/
		find(F, gtFrames, &gtFramesCount, gtFrames, 0);

		trlengtha = 0;
		for(j = 0; j < gtFramesCount; j++){
			if(alltracked[gtFrames[j] * Ngt + i] > 0)
				trlengtha += 1;
		}

		res = trlengtha/(double)gtFramesCount;
		if(res < 0.2){
			(*ML)++;
		}
		else if( t >= gtFrames[gtFramesCount - 1] && res <= 0.8){
			(*PT)++;
		}
		else if(res >= 0.8){
			(*MT)++;
		}
	}

	/*count fragments*/
	FRA = 0;
	for(i = 0; i < Ngt; i++){
		
		FRA += 0;
	}

	/*free temp variables*/
	free(mappedTrackers);
	free(gtFrames);
	free(m);
	free(fp);
	free(allCurTrackedCount);
	free(allTrackedDist);
	free(alldist);
	free(mme);
	free(curTracked);
	free(temp);
	free(curMFound);
	free(XFound);
	free(EsNotMapped);
	free(GTsNotMapped);
	free(g);
	free(mappings);
	free(M);
}