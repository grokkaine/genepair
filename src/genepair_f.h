#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

float* getMIline(const char*,const char*, size_t,int,int,int,int);
float getMIpair(const double*, size_t,int,int,int,int,int);
float getCLRpair(const double*, size_t,int,int,int,int,int);
void computeMI(size_t, int, int, int, int, int, const char*, const char*);
void computeCLR(const char*, size_t, const char*, int);
void computeMIparalel(size_t, int, int, int, const char*, const char*, int, int, int, int);
void getMIrandomness(const char*, int, int, double, int, int);
void computeMIundef(size_t, int, int, int, int, int, const char*, const char*);
//double* readMicroarray(const char*, size_t,int);
