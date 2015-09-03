#include "./genepair_f.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>


int main(int argc, char *argv[]){

  printf("Reading the options; %d\n",argc);
  int ntype = atoi(argv[1]);
  if (ntype == 1){//compute MIs
    char *nmicrf = argv[2];
    char *nmif = argv[3];
    size_t ngenes = atoi(argv[4]);
    int nexp = atoi(argv[5]);
    int gene1 = atoi(argv[6]);
    int gene2 = atoi(argv[7]);
    int numbins = atoi(argv[8]);
    int splineord = atoi(argv[9]);
    computeMI(ngenes, nexp, numbins, splineord, gene1, gene2, nmif, nmicrf);
  }

  if (ntype == 2){//compute clr
    char *nvalf = argv[2];
    size_t ngenes = atoi(argv[3]);
    char* nclrf = argv[4];
    int inputftype = atoi(argv[5]);
    computeCLR(nvalf,ngenes,nclrf,inputftype);
  }

	if (ntype == 3){//compute MI line
		char *nmilinef = argv[2];
		char *nmif = argv[3];
    size_t ngenes = atoi(argv[4]);
    int nexp = atoi(argv[5]);
    int numbins = atoi(argv[6]);
    int splineord = atoi(argv[7]);
		int linei = atoi(argv[8]);
		getMIline(nmilinef, nmif, ngenes, nexp, numbins, splineord, linei);
	}

  if (ntype == 4){//compute MIs in paralel
    char *nmicrf = argv[2];
    char *nmif = argv[3];
    size_t ngenes = atoi(argv[4]);
    int nexp = atoi(argv[5]);
    int numbins = atoi(argv[6]);
    int splineord = atoi(argv[7]);
    int geneR1 = atoi(argv[8]);
    int geneR2 = atoi(argv[9]);
    int geneC1 = atoi(argv[10]);
    int geneC2 = atoi(argv[11]);
    computeMIparalel(ngenes, nexp, numbins, splineord, nmif, nmicrf, geneR1, geneR2, geneC1, geneC2);
  }

//void getMIrandomness(const char* nmif, int npairs, int nexp, double limit, int numbins, int splineord){
  if (ntype == 5){//compute MIs for a number of random pairs between two numbers
    char *nmif = argv[2];
    int npairs = atoi(argv[3]);
    int nexp = atoi(argv[4]);
    double limit = (double) atof(argv[5]);
    int numbins = atoi(argv[6]);
    int splineord = atoi(argv[7]);
    getMIrandomness(nmif, npairs, nexp, limit, numbins, splineord);
  }

  if (ntype == 6){//compute MIs
    char *nmicrf = argv[2];
    char *nmif = argv[3];
    size_t ngenes = atoi(argv[4]);
    int nexp = atoi(argv[5]);
    int gene1 = atoi(argv[6]);
    int gene2 = atoi(argv[7]);
    int numbins = atoi(argv[8]);
    int splineord = atoi(argv[9]);
    computeMIundef(ngenes, nexp, numbins, splineord, gene1, gene2, nmif, nmicrf);
  }

//printf("test X reading: %f %e\n",X[0],X[5*nexp+2]);//X[10*nexp+0]);
//float m = getMIpair(X, ngenes, nexp, ngene1, ngene2, numbins, splineord);
//printf("mi pair: %f\n",m);
//getMIline(X, ngenes, nexp, numbins, splineord, ngene1);
//float z = getCLRpair(X, ngenes, nexp, ngene1, ngene2, numbins, splineord);
//printf("zscore: %f\n",z);
//getMIfile(X, ngenes, nexp, numbins, splineord);
//printf("CLR launched\n");
//computeCLR(ngenes);
//free(X);


}
