#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define UNDEF 1e+30

//#include "genepair_f.h"

double log2d(double x)
{
    return log(x)/log(2);
}

float meanf(float *data, int numSamples)
{
    int curSample;
    float mean = 0;

    for (curSample = 0; curSample < numSamples; curSample++)
    {
        mean += data[curSample];
    }
    return mean / (float) numSamples;
}

float stdf(float *data, int numSamples)
{
    int curSample;
    float m = meanf(data, numSamples);
    float std = 0;

    for (curSample = 0; curSample < numSamples; curSample++)
    {
        std += (data[curSample] - m)*(data[curSample] - m);
    }
    return sqrt(1/(float)(numSamples - 1) * std);
}

/*
Computes max for undefined as well
Use this for undefined!!!
*/
double maxU(const double *data, int numSamples)
{
    int curSample;
    double curMax = -UNDEF;

    for (curSample = 0; curSample < numSamples; curSample++)
    {
        if ( (data[curSample] != UNDEF)&&(data[curSample] > curMax) )
        {
            curMax = data[curSample];
        }
    }
    return curMax;
}

double min(const double *data, int numSamples)
{
    int curSample;
    double curMin = data[0];

    for (curSample = 1; curSample < numSamples; curSample++)
    {
        if (data[curSample] < curMin)
        {
            curMin = data[curSample];
        }
    }
    return curMin;
}

/* "Lifted" from http://local.wasp.uwa.edu.au/~pbourke/curves/spline/ */
/*
  Calculate the blending value, this is done recursively.

  If the numerator and denominator are 0 the expression is 0.
  If the denominator is 0, the expression is 0
  SplineBlend(curBin, splineOrder, knots, z[curSample], numBins);
  http://en.wikipedia.org/wiki/B-spline
  Computes the basis B-spline of degree t and knot u[k] using the Cox-deBoor recursion formula
*/
double SplineBlend(int k,int t, const int *u,double v, int n)
{
    double value = 0.0;
    double d1, d2;

    if (t == 1)
    {
        if (((u[k] <= v) && (v < u[k+1])) ||
                ( fabs(v - u[k+1]) < 1e-10 && (k+1 == n) ) )
            value = 1;
        else
            value = 0;
    }
    else
    {
        d1 = u[k+t-1] - u[k];
        d2 = u[k+t] - u[k+1];

        if ( (d1 == 0) && (d2 == 0) )
            value = 0;
        else if (d1 == 0)
            value = (u[k+t] - v) / (double)d2 * SplineBlend(k+1,t-1,u,v,n);
        else if (d2 == 0)
            value = (v - u[k]) / (double)d1 * SplineBlend(k,t-1,u,v,n);
        else
            value = (v - u[k]) / (double)d1 * SplineBlend(k,t-1,u,v,n) +
                    (u[k+t] - v) / (double)d2 * SplineBlend(k+1,t-1,u,v,n);
    }
    if (value < 0)
    {
        value = 0; /* rounding sometimes makes this < 0, e.g. -0.000000001 */
    }
    return(value);
}

/*
  The positions of the subintervals of v and breakpoints, the position
  on the curve are called knots. Breakpoints can be uniformly defined
  by setting u[j] = j, a more useful series of breakpoints are defined
  by the function below. This set of breakpoints localises changes to
  the vicinity of the control point being modified.
  int *knots = calloc(numBins + splineOrder, sizeof(int));
  SplineKnots(knots, numBins, splineOrder);
  positions on hte curve 0 1 ... 10 11 11 11 ..7 times (10 bins, 7 order spline)
*/
void SplineKnots(int *u,int n,int t)
{
    int j;
    int d=n-1;

    for (j=0;j<=d + t;j++)
    {
        if (j < t)
            u[j] = 0;
        else if (j <= d)
            u[j] = u[j-1] + 1;
        else if (j > d)
            u[j] = u[d] + 1;
    }
}

/*scales the expressions to between [0,1]*(numBins-splineOrder+1)
to work with undefined, the maximum function was changed and z is set to 1e+30 for undefined x
this will have to change findWeights accordingly
*/
void xToZ(const double *fromData, double *toData, int numSamples, int splineOrder, int numBins)
{
    int curSample;
    double xMax = maxU(fromData, numSamples);
    double xMin = min(fromData, numSamples);
    //printf("min: %f% max: %f\n",xMin,xMax);

    for (curSample = 0; curSample < numSamples; curSample++)
    {
        if (fromData[curSample] == UNDEF)
        {
            toData[curSample]=UNDEF;
        }
        else
        {
            toData[curSample] = (fromData[curSample] - xMin) * (numBins - splineOrder + 1) / (double) (xMax - xMin);
        }
    }
}

/*
Computes the smoothed probability that the expression x is in a certain bin
To work with undef, the weight is set to 1e+30 for all undefined conditions across the whole bins
This will have to change histogram computations accordingly (hist1 and hist2)
*/
void findWeights(const double *x, const int *knots, double *weights, int numSamples, int splineOrder, int numBins)
{
    int curSample;
    int curBin;
    double *z = (double*) calloc(numSamples, sizeof(double));

//    for (curSample = 0; curSample < numBins + splineOrder; curSample++) {
//        printf("knot %d: %d\n", curSample, knots[curSample]);/*0 0 1 2 3 3 3*/
//    }

    /*todo: xToZ doesnt work with undefined*/
    xToZ(x, z, numSamples, splineOrder, numBins); /*scales the expressions to between [0,1]*(numBins-splineOrder+1) - WHY? - because its the Daub et al algorithm*/

    //printf("x[0]:\n");snPrintString(x, numSamples);
    for (curSample = 0; curSample < numSamples; curSample++)
    {
        if (z[curSample] == UNDEF)
        {
            for (curBin = 0; curBin < numBins; curBin++)
            {
                weights[curBin * numSamples + curSample] = UNDEF;
                //printf("b:%d,s:%d|%e(%e)\t", curBin, curSample, weights[curBin * numSamples + curSample],z[curSample]);
            }
            //exit(0);
        }
        else
        {
            for (curBin = 0; curBin < numBins; curBin++)
            {
                weights[curBin * numSamples + curSample] = SplineBlend(curBin, splineOrder, knots, z[curSample], numBins);
                //printf("b:%d,s:%d|%e(%e)\t", curBin, curSample, weights[curBin * numSamples + curSample],z[curSample]);
            }
        }
    }
    free(z);
}

/*
Changed to take into account the number of defined variables
The probabilities will be scaled by the numbed of defined samples not by the number of samples
*/
void hist1d(const double *x, const int *knots, double *hist, double *w, int numSamples, int splineOrder, int numBins)
{
    int curSample;
    int curBin;

    /*first get the number of defined samples, compute it from bin curBin=0 since should be the same in every bin*/
    int nbDef = 0;
    for (curSample = 0; curSample < numSamples; curSample++){
        if (w[curSample]!=UNDEF) nbDef += 1;
    }
    for (curBin = 0; curBin < numBins; curBin++)
    {
        for (curSample = 0; curSample < numSamples; curSample++)
        {
            if (w[curBin * numSamples + curSample] != UNDEF)
            {
                hist[curBin] += w[curBin * numSamples + curSample]/(double)nbDef;
            }
        }
        //printf("hist[%d]:%f, ", curBin, hist[curBin]);
    }
    //exit(0);
}

void histU(const double *x, const double *y, const int *knots, const double *wx, const double *wy, double *hist12, double *hist1, double *hist2, int numSamples, int splineOrder, int numBins)
{
    int curSample;
    int curBin, curBinX, curBinY;
    for (curBin = 0; curBin < numBins; curBin++)
    {
        int nbDef = 0;
        for (curSample = 0; curSample < numSamples; curSample++)
            if ((wx[curBin * numSamples + curSample] != UNDEF) && (wy[curBin * numSamples + curSample] != UNDEF)) nbDef += 1;
        hist1[curBin] = 0.0;
        hist2[curBin] = 0.0;
        for (curSample = 0; curSample < numSamples; curSample++)
            if ((wx[curBin * numSamples + curSample] != UNDEF) && (wy[curBin * numSamples + curSample] != UNDEF))
            {
                hist1[curBin] += wx[curBin * numSamples + curSample]/(double)nbDef;
                hist2[curBin] += wy[curBin * numSamples + curSample]/(double)nbDef;
            }
        //printf("hist1[%d]:%f, ", curBin, hist1[curBin]);
        //printf("hist2[%d]:%f, ", curBin, hist2[curBin]);
    }
    for (curBinX = 0; curBinX < numBins; curBinX++)
        for (curBinY = 0; curBinY < numBins; curBinY++)
        {
            int nbDef = 0;
            for (curSample = 0; curSample < numSamples; curSample++)
                if ((wx[curBinX * numSamples + curSample]!=UNDEF) && (wy[curBinY * numSamples + curSample]!=UNDEF)) nbDef += 1;
            hist12[curBinX * numBins + curBinY] = 0.0;
            for (curSample = 0; curSample < numSamples; curSample++)
                if ((wx[curBinX * numSamples + curSample]!=UNDEF) && (wy[curBinY * numSamples + curSample]!=UNDEF))
                    hist12[curBinX * numBins + curBinY] += wx[curBinX * numSamples + curSample] * wy[curBinY * numSamples + curSample]/(double)nbDef;
            //printf("hist12[%d,%d]:%f, ", curBinX, curBinY, hist12[curBinX * numBins + curBinY]);
        }
}

/*
Changed to take into account the number of defined variables
The probabilities will be scaled by the numbed of commonly defined samples not by the number of samples
*/
void hist2d(const double *x, const double *y, const int *knots, const double *wx, const double *wy, double *hist, int numSamples, int splineOrder, int numBins)
{
    int curSample;
    int curBinX, curBinY;

    /*findWeights(x, knots, wx, numSamples, splineOrder, numBins);
      findWeights(y, knots, wy, numSamples, splineOrder, numBins); */

    int nbDef = 0;
    for (curSample = 0; curSample < numSamples; curSample++){
        if ((wx[curSample]!=UNDEF) && (wy[curSample]!=UNDEF)) nbDef += 1;
    }
    for (curBinX = 0; curBinX < numBins; curBinX++)
    {
        for (curBinY = 0; curBinY < numBins; curBinY++)
        {
            for (curSample = 0; curSample < numSamples; curSample++)
            {
                if ((wx[curBinX * numSamples + curSample]!=UNDEF) && (wy[curBinY * numSamples + curSample]!=UNDEF))
                {
                    hist[curBinX * numBins + curBinY] += wx[curBinX * numSamples + curSample] * wy[curBinY * numSamples + curSample]/(double)nbDef;
                }
            }
            //printf("hist[%d,%d]:%f, ", curBinX, curBinY, hist[curBinX * numBins + curBinY]);
        }
        //exit(0);
    }

    /*
      for (curBinX = 0; curBinX < numBinsX; curBinX++) {
      for (curBinY = 0; curBinY < numBinsY; curBinY++) {
      mexPrintf("%f\t", hist[curBinX * numBinsY + curBinY]);
      }
      mexPrintf("\n");
      }
    */
}

double entropy1d(const double *x, const int *knots, double *weights, int numSamples, int splineOrder, int numBins)
{
    int curBin;
    double *hist = (double*) calloc(numBins, sizeof(double));
    double H = 0;

    hist1d(x, knots, hist, weights, numSamples, splineOrder, numBins);
    for (curBin = 0; curBin < numBins; curBin++)
    {
        if (hist[curBin] > 0)
        {
            H -= hist[curBin] * log2d(hist[curBin]);
        }
    }
    free(hist);
    return H;
}

double entropy2d(const double *x, const double *y, const int *knots, const double *wx, const double *wy, int numSamples, int splineOrder, int numBins)
{
    int curBinX, curBinY;
    double *hist = (double*) calloc(numBins * numBins, sizeof(double));
    double H = 0;
    double incr;

    hist2d(x, y, knots, wx, wy, hist, numSamples, splineOrder, numBins);
    for (curBinX = 0; curBinX < numBins; curBinX++)
    {
        for (curBinY = 0; curBinY < numBins; curBinY++)
        {
            incr = hist[curBinX * numBins + curBinY];
            if (incr > 0)
            {
                H -= incr * log2d(incr);
            }
        }
    }
    free(hist);
    return H;
}

double miundefined(const double *x, const double *y, const int *knots, const double *wx, const double *wy, int numSamples, int splineOrder, int numBins)
{
    int curBin, curBinX, curBinY;
    double *h1 = (double*) calloc(numBins, sizeof(double));
    double *h2 = (double*) calloc(numBins, sizeof(double));
    double *h12 = (double*) calloc(numBins * numBins, sizeof(double));
    double e1 = 0;
    double e2 = 0;
    double e12 = 0;
    double p;
    histU(x, y, knots, wx, wy, h12, h1, h2, numSamples, splineOrder, numBins);
    for (curBin = 0; curBin < numBins; curBin++)
    {
        p = h1[curBin];
        if (p > 0) e1 -= p*log2d(p);
        p = h2[curBin];
        if (p > 0) e2 -= p*log2d(p);
    }

    for (curBinX = 0; curBinX < numBins; curBinX++)
        for (curBinY = 0; curBinY < numBins; curBinY++)
        {
            p = h12[curBinX * numBins + curBinY];
            if (p > 0) e12 -= p*log2d(p);
        }
    free(h12);
    free(h1);
    free(h2);
    return e1+e2-e12;
}

void snPrintArrayDouble(const char *fn, const double *x, int size)
{
    FILE *f;
    int i;
    f = fopen(fn,"w");
    for (i = 0; i < size; i++)
    {
        if (i%10==0) fprintf(f,"\n");
        fprintf(f, "%f ", x[i]);
    }
    fclose(f);
}

void snPrintArrayInt(const char *fn, const int *x, int size)
{
    FILE *f;
    int i;
    f = fopen(fn,"w");
    for (i = 0; i < size; i++)
    {
        if (i%10==0) fprintf(f,"\n");
        fprintf(f, "%d ", x[i]);
    }
    fclose(f);
}

double* readMicroarray(const char *nmicrf, int ngenes, int nexp){
    printf("Alocating the microarray\n");
    double *X = (double*) calloc(ngenes*nexp,sizeof(double));
    FILE *mif;
    mif = fopen(nmicrf,"r");
    int i,j;
    for (i = 0; i < ngenes; i++){
        for (j = 0; j < nexp; j++){
            double t;
            char tt[80];
            fscanf(mif,"%s",tt);
            t = atof(tt);
            X[i*nexp + j] = t;
        }
    }
    fclose(mif);
    printf("Done alocating the microarray\n");
    return X;
}

float getMIpair(const double* X, int ngenes, int nexp, int row, int col, int numbins, int splineord){
    int *knots = calloc(numbins + splineord, sizeof(int));
    SplineKnots(knots, numbins, splineord);
    /* allocate room for marginal weights one weight for each bin in each each condition point*/
    double *weights1 = (double*) calloc(nexp * numbins, sizeof(double));
    double *weights2 = (double*) calloc(nexp * numbins, sizeof(double));
    /*this is the core function, it computes the weights using the spline smoothing, not from the expression itself but from an xToz transformaiton of it (step 1a in the article)*/
    findWeights(X + col*nexp, knots, weights1, nexp, splineord, numbins);
    findWeights(X + row*nexp, knots, weights2, nexp, splineord, numbins);
    float e1 = (float)entropy1d(X + col*nexp,knots,weights1,nexp,splineord,numbins);
    float e2 = (float)entropy1d(X + row*nexp,knots,weights2,nexp,splineord,numbins);
    float e12 = (float)entropy2d(X + col*nexp, X + row*nexp, knots, weights1, weights2, nexp, splineord, numbins);
    float mi = e1 + e2 - e12;
//    if (col == 0){
//        snPrintArrayDouble("weight_final.txt",weights2,nexp * numbins);
//        snPrintArrayDouble("micr_final.txt",X + row*nexp,nexp);
//        snPrintArrayInt("knots_final.txt",knots,numbins + splineord);
//        printf("Num bins %d, spline order %d",numbins,splineord);
//    }
    free(knots);
    free(weights1);
    free(weights2);
    return mi;
}

float* getMIline(const char* nmilinef, const char* nmif, int ngenes, int nexp, int numbins, int splineord, int linei){
    int i;
	double *X = readMicroarray(nmif,ngenes,nexp);
    float *row = (float*) calloc(ngenes, sizeof(float));
	FILE *fo;
	fo = fopen(nmilinef,"w");
    for (i = 0; i<ngenes; i++)
    {
        row[i] = getMIpair(X,ngenes,nexp,linei,i,numbins,splineord);
		fprintf(fo,"%f",row[i]);
		if (i<ngenes-1) fprintf(fo," ");
    }
	fprintf(fo,"\n");
	fclose(fo);
	free(X);
    return row;
}

// float getCLRpair(const double* X, int ngenes, int nexp, int g1, int g2, int numbins, int splineord){
    // float *mi1, *mi2;
    // float m1, m2, s1, s2;
    // float z, zg1, zg2;
    // mi1 = getMIline(X, ngenes, nexp, numbins, splineord, g1);
    // mi2 = getMIline(X, ngenes, nexp, numbins, splineord, g2);
    // mi1[g1]=0; mi2[g2]=0;
    // m1 = meanf(mi1,ngenes); m2 = meanf(mi2,ngenes);
    // s1 = stdf(mi1,ngenes); s2 = stdf(mi2,ngenes);
    // zg1 = (mi1[g2]-m1)/s1; zg2 = (mi2[g1]-m2)/s2;
    // //printf("%f %f\n",mi1[g2],mi2[g1]);
    // if (zg1 < 0) zg1 = 0; if (zg2 < 0) zg2 = 0;
    // z = sqrt(zg1*zg1 + zg2*zg2);
    // printf("CLR Z-score: %f",z);
    // free(mi1); free(mi2);
    // return z;
// }

/*
computes the MI scores for a list of gene pairs
*/
void getMIfile(const double* X, int ngenes, int nexp, int numbins, int splineord){
    int i;
    int nlines = 36;
    FILE *fi, *fo;
    fi = fopen("genepairs.txt","r");
    fo = fopen("genepairsmi.txt","w");
    for (i=0; i< nlines; i++){
        int ng1, ng2;
        fscanf(fi,"%d %d",&ng1,&ng2);
        float mi = getMIpair(X,ngenes,nexp,ng1,ng2,numbins,splineord);
        fprintf(fo,"%d %d %f\n",ng1,ng2,mi);
    }
    fclose(fi); fclose(fo);
}

void loadTorgeirUpperMatrix(float* v, int ngenes, const char *nvalf){
    int i,j;
    FILE *fi;
    fi = fopen(nvalf,"r");
    char c[200];

    for (i=0; i<ngenes+1; i++) fscanf(fi,"%s",c);//read the header
    for (i=0; i<ngenes; i++){
        fscanf(fi,"%s",c);//read the gene name
        fscanf(fi,"%s",c);//read all spaces and the 1
        for (j=i+1; j<ngenes; j++) fscanf(fi,"%f",&v[(j-1)*j/2+i]);
    }
    fclose(fi);

//    for (i=0; i<ngenes; i++){
//        for (j=0; j < i; j++)
//            printf("%f ",v[i*(i-1)/2+j]);
//        printf("\n");
//    }
//
//    FILE *ft;
//    ft = fopen("torgeir_genereadtest.txt","a");
//    i = 3;
//    for (j = 0; j<i; j++) fprintf(ft,"%f ",fabs(v[i*(i-1)/2+j]));
//    fprintf(ft,"%f ",0);
//    for (j = i+1; j<ngenes; j++) fprintf(ft,"%f ",fabs(v[j*(j-1)/2+i]));
//    fprintf(ft,"\n");
//    i = 6;
//    for (j = 0; j<i; j++) fprintf(ft,"%f ",fabs(v[i*(i-1)/2+j]));
//    fprintf(ft,"%f ",0);
//    for (j = i+1; j<ngenes; j++) fprintf(ft,"%f ",fabs(v[j*(j-1)/2+i]));
//    fprintf(ft,"\n");
//    fclose(ft);

    printf("Done reading correlations\n");
}

void loadLowerMatrix(float* v, int ngenes, const char *nvalf){
    int i,j;
    FILE *fi;
    fi = fopen(nvalf,"r");
    for (i=1; i<ngenes; i++)
        for (j=0; j<i; j++) fscanf(fi,"%f",&v[(i-1)*i/2+j]);
    fclose(fi);
    printf("v1: %f\n",v[1]);
    printf("Done reading correlations\n");
}

void computeCLR(const char *nvalf, int ngenes, const char *nclrf, int inputftype){
    int i,j;
    float *v = (float*) calloc(ngenes*(ngenes-1)/2,sizeof(float));
    if (inputftype == 0) loadTorgeirUpperMatrix(v, ngenes, nvalf);
    if (inputftype == 1) loadLowerMatrix(v, ngenes, nvalf);
    //printf("v2: %f\n",v[1]);

    float *m = (float*) calloc(ngenes,sizeof(float));
    float *s = (float*) calloc(ngenes,sizeof(float));
    for (i=0; i<ngenes; i++){
        float *tv = (float*) calloc(ngenes,sizeof(float));
        for (j = 0; j<i; j++) tv[j]=fabs(v[i*(i-1)/2+j]);
        tv[i]=0;
        for (j = i+1; j<ngenes; j++) tv[j]=fabs(v[j*(j-1)/2+i]);
        //for (j = 0; j<ngenes; j++) printf("%f ",tv[j]);
        //printf("\n");
        m[i] = meanf(tv,ngenes);
        s[i] = stdf(tv,ngenes);
        free(tv);
    }

    FILE *fo;
    fo = fopen(nclrf,"w");
    for(i=0; i<ngenes; i++){
        for (j=0; j<i; j++){
            float c = fabs(v[i*(i-1)/2 + j]);
            float zi = (c - m[i])/s[i];
            float zj = (c - m[j])/s[j];
            if (zi < 0) zi = 0; if (zj < 0) zj = 0;
            float z = sqrt(zi*zi + zj*zj);
            fprintf(fo,"%f",z);
            if (j<i) fprintf(fo," ");
        }
        if (i>0) fprintf(fo,"\n");
    }
    fclose(fo);
    free(m);free(s);
    free(v);
    printf("Done computing CLR matrix!\n");
}

/* Compute mutual information for lower triangular matrix from col fromCol to col toCol
 studied for 10 bins, 7 order spline

BMC Bioinformatics 2004
Estimating mutual information using B-spline functions – an
improved similarity measure for analysing gene expression data
Carsten O Daub*1,4, Ralf Steuer2, Joachim Selbig1 and Sebastian Kloska1,3
*/
void computeMI(int ngenes, int nexp, int numbins, int splineord, int fromRow, int toRow, const char* nmif, const char* nmicrf){
    int i, j;
    FILE *fo, *ft;
    char ftn[100];
    strcat(ftn,nmif);
    strcat(ftn,".mi_log.txt");
    double *X = readMicroarray(nmicrf,ngenes,nexp);
    /*delineate the knots*/
    int *knots = calloc(numbins + splineord, sizeof(int));
    double *entropies = calloc(ngenes, sizeof(double));
    /* allocate room for marginal weights one weight for each bin in each each condition point*/
    double **weights = (double**)calloc(ngenes, sizeof(double*));    /*each weight is a matrix num bins * num condiitons*/

    /*fills the knots array positions on hte curve 0 1 ... 10 11 11 11 ..7 times (10 bins, 7 order spline)*/
    SplineKnots(knots, numbins, splineord);

    /*computes the weights and fills the 1d entropies*/
    for (i = 0; i < ngenes; i++)//ngenes-toRow would have been better
    {
        /* allocate room for marginal weights one weight for each bin in each each condition point*/
        weights[i] = (double*) calloc(nexp * numbins, sizeof(double));
        /*this is the core function, it computes the weights using the spline smoothing, not from the expression itself but from an xToz transformaiton of it (step 1a in the article)*/
        findWeights(X + i*nexp, knots, weights[i], nexp, splineord, numbins);
        /*computes the 1d entropies suming up probabilities over the bins, each probability is the scaled sum of the weights for for all the conditions and the current bin*/
        entropies[i] = entropy1d(X + i*nexp, knots, weights[i], nexp, splineord, numbins);/*first parameter, data, is not used, since weights are already computed*/
        ft=fopen(ftn,"a"); fprintf(ft,"weight %d computed\n",i); fclose(ft);
    }
    printf("Done computing the weight arrays!\n");
    fo = fopen(nmif,"a");
    float e1, e2, e12, mi;
    for (i=fromRow; i<toRow; i++){
        e1 = entropies[i];
        for (j=0; j<i; j++){
            e2 = entropies[j];
            e12 = (float)entropy2d(X + i*nexp, X + j*nexp, knots, weights[i], weights[j], nexp, splineord, numbins);
            mi = e1 + e2 - e12;
            fprintf(fo,"%f",mi);
            if (j<i) fprintf(fo," ");
        }
        if (i>0) fprintf(fo,"\n");
        fclose(fo);
        fo = fopen(nmif,"a");
    }
    fclose(fo);
    /* free marginal weights */
    for (i = 0; i < ngenes; i++) free(weights[i]);
    free(weights);
    free(entropies);
    free(knots);
    free(X);
}

/*
variation of computeMI where the one-dimensional entropies are computed only for the commonly defined samples
*/
void computeMIundef(int ngenes, int nexp, int numbins, int splineord, int fromRow, int toRow, const char* nmif, const char* nmicrf){
    int i, j;
    FILE *fo;
    //FILE *ft;
    //char ftn[100];
    //strcat(ftn,nmif);
    //strcat(ftn,".mi_log.txt");
    double *X = readMicroarray(nmicrf,ngenes,nexp);
    /*delineate the knots*/
    int *knots = calloc(numbins + splineord, sizeof(int));
    /* allocate room for marginal weights one weight for each bin in each each condition point*/
    double **weights = (double**)calloc(ngenes, sizeof(double*));    /*each weight is a matrix num bins * num condiitons*/

    /*fills the knots array positions on hte curve 0 1 ... 10 11 11 11 ..7 times (10 bins, 7 order spline)*/
    SplineKnots(knots, numbins, splineord);

    /*computes the weights and fills the 1d entropies*/
    for (i = 0; i < ngenes; i++)//ngenes-toRow would have been better
    {
        /* allocate room for marginal weights one weight for each bin in each each condition point*/
        weights[i] = (double*) calloc(nexp * numbins, sizeof(double));
        /*this is the core function, it computes the weights using the spline smoothing, not from the expression itself but from an xToz transformaiton of it (step 1a in the article)*/
        findWeights(X + i*nexp, knots, weights[i], nexp, splineord, numbins);
        //ft=fopen(ftn,"a"); fprintf(ft,"weight %d computed\n",i); fclose(ft);
    }
    printf("Done computing the weight arrays!\n");
    fo = fopen(nmif,"a");
    float e1, e2, e12, mi;
    for (i=fromRow; i<toRow; i++){
        for (j=0; j<i; j++){
            mi = (float)miundefined(X + i*nexp, X + j*nexp, knots, weights[i], weights[j], nexp, splineord, numbins);
            fprintf(fo,"%f",mi);
            if (j<i) fprintf(fo," ");
        }
        if (i>0) fprintf(fo,"\n");
        fclose(fo);
        fo = fopen(nmif,"a");
    }
    fclose(fo);
    /* free marginal weights */
    for (i = 0; i < ngenes; i++) free(weights[i]);
    free(weights);
    free(knots);
    free(X);
}

void computeMIparalel(int ngenes, int nexp, int numbins, int splineord, const char* nmif, const char* nmicrf, int fromRow, int toRow, int fromCol, int toCol){
    int i, j;
    FILE *fo, *ft;
    //char ftn[100];
    char *ftn;
    strcat(ftn,nmif);
    strcat(ftn,".mi_log.txt");
    double *X = readMicroarray(nmicrf,ngenes,nexp);
    /*delineate the knots*/
    int *knots = calloc(numbins + splineord, sizeof(int));
    double *entropies = calloc(ngenes, sizeof(double));
    /* allocate room for marginal weights one weight for each bin in each each condition point*/
    double **weights = (double**)calloc(ngenes, sizeof(double*));    /*each weight is a matrix num bins * num condiitons*/
    /*fills the knots array positions on hte curve 0 1 ... 10 11 11 11 ..7 times (10 bins, 7 order spline)*/
    SplineKnots(knots, numbins, splineord);
    /*computes the weights and fills the 1d entropies*/
    for (i = fromRow; i < toRow; i++)//ngenes-toRow would have been better
    {
        /* allocate room for marginal weights one weight for each bin in each each condition point*/
        weights[i] = (double*) calloc(nexp * numbins, sizeof(double));
        /*this is the core function, it computes the weights using the spline smoothing, not from the expression itself but from an xToz transformaiton of it (step 1a in the article)*/
        findWeights(X + i*nexp, knots, weights[i], nexp, splineord, numbins);
        /*computes the 1d entropies suming up probabilities over the bins, each probability is the scaled sum of the weights for for all the conditions and the current bin*/
        entropies[i] = entropy1d(X + i*nexp, knots, weights[i], nexp, splineord, numbins);/*first parameter, data, is not used, since weights are already computed*/
        ft=fopen(ftn,"a"); fprintf(ft,"weight row %d computed\n",i); fclose(ft);
    }
    for (i = fromCol; i < toCol; i++)//ngenes-toRow would have been better
    {
        /* allocate room for marginal weights one weight for each bin in each each condition point*/
        weights[i] = (double*) calloc(nexp * numbins, sizeof(double));
        /*this is the core function, it computes the weights using the spline smoothing, not from the expression itself but from an xToz transformaiton of it (step 1a in the article)*/
        findWeights(X + i*nexp, knots, weights[i], nexp, splineord, numbins);
        /*computes the 1d entropies suming up probabilities over the bins, each probability is the scaled sum of the weights for for all the conditions and the current bin*/
        entropies[i] = entropy1d(X + i*nexp, knots, weights[i], nexp, splineord, numbins);/*first parameter, data, is not used, since weights are already computed*/
        ft=fopen(ftn,"a"); fprintf(ft,"weight col %d computed\n",i); fclose(ft);
    }

    printf("Done computing the weight arrays!\n");
    fo = fopen(nmif,"a");
    float e1, e2, e12, mi;
    for (i=fromRow; i<toRow; i++){
        e1 = entropies[i];
        for (j=fromCol; j<toCol; j++){
            if (j<i){
                e2 = entropies[j];
                e12 = (float)entropy2d(X + i*nexp, X + j*nexp, knots, weights[i], weights[j], nexp, splineord, numbins);
                mi = e1 + e2 - e12;
                fprintf(fo,"%f ",mi);
            }
        }
        if (i>0) fprintf(fo,"\n");
        fclose(fo);
        fo = fopen(nmif,"a");
    }
    fclose(fo);
    /* free marginal weights */
    for (i = 0; i < ngenes; i++) free(weights[i]);
    free(weights);
    free(entropies);
    free(knots);
    free(X);
}

void getMIrandomness(const char* nmif, int npairs, int nexp, double limit, int numbins, int splineord){
  int i,ipair;
  FILE *fo;
  fo = fopen(nmif,"a");
  for (ipair = 0; ipair<npairs; ipair++){
    double *x = calloc(nexp, sizeof(double));
    double *y = calloc(nexp, sizeof(double));
    //srand((unsigned)time(NULL));
    srand((unsigned)ipair*(unsigned)time(NULL));
    for (i = 0; i<nexp; i++){
      x[i] = ((double)rand()/(double)RAND_MAX)*limit;
      y[i] = ((double)rand()/(double)RAND_MAX)*limit;
      //printf("%f ",x[i]);
    }
    int *knots = calloc(numbins + splineord, sizeof(int));
    SplineKnots(knots, numbins, splineord);
    /* allocate room for marginal weights one weight for each bin in each each condition point*/
    double *weights1 = (double*) calloc(nexp * numbins, sizeof(double));
    double *weights2 = (double*) calloc(nexp * numbins, sizeof(double));
    /*this is the core function, it computes the weights using the spline smoothing, not from the expression itself but from an xToz transformaiton of it (step 1a in the article)*/
    findWeights(x, knots, weights1, nexp, splineord, numbins);
    findWeights(y, knots, weights2, nexp, splineord, numbins);
    float e1 = (float)entropy1d(x,knots,weights1,nexp,splineord,numbins);
    float e2 = (float)entropy1d(y,knots,weights2,nexp,splineord,numbins);
    float e12 = (float)entropy2d(x, y, knots, weights1, weights2, nexp, splineord, numbins);
    float mi = e1 + e2 - e12;
    fprintf(fo,"%f\n",mi);
    free(x);
    free(y);
    free(knots);
    free(weights1);
    free(weights2);
  }
  fclose(fo);
  return;
}
