#ifndef UTILE_H
#define UTILE_H

#define SUCCESS 0

#define FAILURE -1

#define MAX_ATTR 1000000
#define MIN_ATTR -1000000

#include <vector>
#include <iostream>
using namespace std;

#define PI 3.14159

void SetBit(unsigned int *a1, int objectid, int noCells, int noTotalObjects);
void UnsetBit(unsigned int *a1, int objectid, int noCells, int noTotalObjects);
int  GetBit(unsigned int *a1, int objectid, int noCells, int noTotalObjects);
void IntersectionArrayBits(unsigned int *res, unsigned int *a1, unsigned int *a2, int noCells);
void UnionArrayBits(unsigned int *res, unsigned int *a1, unsigned int *a2, int noCells);
void NegationArrayBits(unsigned int *res, unsigned int *a1, int noCells);
int GetNumberOneBits(unsigned int *a1, int noCells, int noTotalObjects);
vector<int> GetOneBits(unsigned int *a1, int noCells, int noTotalObjects);
void CopyArrayBits(unsigned int *res, unsigned int *a1, int noCells);
int TestEqualityArrayBits(unsigned int *a1, int noCells1, unsigned int *a2, int noCells2);


double Poisson(int x, double m);
int GetCriticalValuePoissonRight(double dMean, double dAlpha);
int GetCriticalValuePoissonLeft(double dMean, double dAlpha);

bool IsInf(double value);
bool IsNan(double val);
double Choose(int d, int p);
double Choose2(int n, int k);
double Binomial(int k, int n, double p);
double Binomial2(int k, int n, double p);
int GetCriticalValueBinomialRight(int n, double p, double dAlpha);
int GetCriticalValueBinomialRight2(int n, double p, double dAlpha);
int GetCriticalValueBinomialLeft(int n, double p, double dAlpha);
int GetCriticalValueBinomialLeft2(int n, double p, double dAlpha);

#ifndef __GNUG__
	double erf(double x);
#endif //__GNUG__

double gamma(double a, double x);

double Hypergeometric(int k, int M, int d);
int GetCriticalValueHypergeometricRight(int M, int d , double dAlpha);

double GetCriticalValueKolmogorov(double t_star, double dAlpha);

double GetRightCriticalValueStandardNormal(double dAlpha);

int ExistElem(int elem, vector<int> &v);
int ExistElem(int elem, int* v, int no_elem);
int AreVectorsEqual(int *v1, int l1, int *v2, int l2);
int AreVectorsEqual(vector<int> &v1, vector<int> &v2);
int CheckInclusionVectors(vector<int> &v1, vector<int> &v2);

#endif


