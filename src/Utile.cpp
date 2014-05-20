#include "Utile.h"

#include <math.h>
#include <limits.h>
#include <float.h>
#include <errno.h>
#include <assert.h>

/************
Functions for working with the representation on bits
*********** */

// Set the bit on position 'objectid' to 1 
void SetBit(unsigned int *a1, int objectid, int noCells, int noTotalObjects){
    assert( (objectid >= 0) && (objectid < noTotalObjects) );
	int cell = objectid / (8*sizeof(unsigned int));
	int pos_inside_cell = objectid - cell * (8*sizeof(unsigned int));

	unsigned int mask = 1 << pos_inside_cell;
	a1[cell] |= mask;
}

// Set the bit on position 'objectid' to 0
void UnsetBit(unsigned int *a1, int objectid, int noCells, int noTotalObjects){
    assert( (objectid >= 0) && (objectid < noTotalObjects) );
	int cell = objectid / (8*sizeof(unsigned int));
	int pos_inside_cell = objectid - cell * (8*sizeof(unsigned int));

	unsigned int mask = ~(1 << pos_inside_cell);
	a1[cell] &= mask;
}

// Returns the bit on position 'objectid' 
int GetBit(unsigned int *a1, int objectid, int noCells, int noTotalObjects){
    unsigned int value = 0;
    assert( (objectid >= 0) && (objectid < noTotalObjects) );
	
	int cell = objectid / (8*sizeof(unsigned int));
	int pos_inside_cell = objectid - cell * (8*sizeof(unsigned int));

	unsigned int mask = 1 << pos_inside_cell;
	value = a1[cell] & mask;
	value >>= pos_inside_cell;

	return value;
}

void IntersectionArrayBits(unsigned int *res, unsigned int *a1, unsigned int *a2, int noCells){
    int i;
    for(i=0; i<noCells; i++)
        res[i] = a1[i] & a2[i];
}

void UnionArrayBits(unsigned int *res, unsigned int *a1, unsigned int *a2, int noCells){
    int i;
    for(i=0; i<noCells; i++)
        res[i] = a1[i] | a2[i];
}

void NegationArrayBits(unsigned int *res, unsigned int *a1, int noCells){
    int i;
    for(i=0; i<noCells; i++){
        unsigned int negation = ~a1[i];
        res[i] = negation;
    }
}

int GetNumberOneBits(unsigned int *a1, int noCells, int noTotalObjects){
    int noCommonElems = 0;
    for(int i=0; i<noCells; i++){
        unsigned int mask = 1;
        for(unsigned int j=0; j<8*sizeof(unsigned int); j++){
            //if (temp & mask){
            //    cout << i*nocells+j << " ";
            //}
            int bit = i*8*sizeof(unsigned int) + j;
            if (bit < noTotalObjects)
                noCommonElems += ((a1[i] & mask) != 0);
            mask <<= 1;
        }
    }
    return noCommonElems;
}

vector<int> GetOneBits(unsigned int *a1, int noCells, int noTotalObjects){
    vector<int> res;
    for(int i=0; i<noCells; i++){
        unsigned int mask = 1;
        for(unsigned int j=0; j<8*sizeof(unsigned int); j++){
            //if (temp & mask){
            //    cout << i*nocells+j << " ";
            //}
            if ((a1[i] & mask) != 0){
                int bit = i*8*sizeof(unsigned int) + j;
                if (bit < noTotalObjects)
                    res.push_back(bit);
            }
            mask <<= 1;
        }
    }
    return res;
}

void CopyArrayBits(unsigned int *res, unsigned int *a1, int noCells){
    for(int i=0; i<noCells; i++){
        res[i] = a1[i];
    }
}

int TestEqualityArrayBits(unsigned int *a1, int noCells1, unsigned int *a2, int noCells2){

	if (noCells1 != noCells2)
		return 0;
	
	int bSame = 1;
	for(int i = 0; i < noCells1 && bSame; i++){
		if (a1[i] != a2[i]){
			bSame = 0;
		}
	}

	return bSame;
}

////////////////////////////////   POISSON RELATED FUNCTIONS //////////////////

/************ 
Compute Poisson p.d.f.
************/
double Poisson(int x, double m){
    double res; 
    
    //scalling trick
    //if (x>10000){
    //   x = (int)((double)x/10.0);
    //   m = m/10.0;
    //}
    
    double quantity = -m + (double)x * log(m);
    for(int i=1; i<=x; i++){
        quantity -= log((double)i);
    }
    res = exp(quantity);
    
    return res;    
}

/************ 
Compute the critical value at right of a Poisson distribution with mean = dMean at significance level dAlpha
************/
int GetCriticalValuePoissonRight(double dMean, double dAlpha){

	
	int iCriticalValue = 1; //always return a critical value of at least 1

	if (dAlpha <= 0){
		cerr << "Statistical significance level <= 0" << endl;
		return iCriticalValue;
	}
		

	int bSatisfyCondition = 0;
	double dSum = Poisson(0, dMean);
	while(!bSatisfyCondition){
	
		dSum += Poisson(iCriticalValue, dMean);
		if (dSum >= (1 - dAlpha)){
			bSatisfyCondition = 1;
		}
		else{
			iCriticalValue++;
		}
	}

	return iCriticalValue;

}

/************ 
Compute the critical value at left of a Poisson distribution with mean = dMean at significance level dAlpha
************/
int GetCriticalValuePoissonLeft(double dMean, double dAlpha){

	int iCriticalValue = 1; //always return a critical value of at least 1
	if (dAlpha <= 0){
		cerr << "Statistical significance level <= 0" << endl;
		return iCriticalValue;
	}
	
	int bSatisfyCondition = 1;
	double dSum = Poisson(0, dMean);
	while(bSatisfyCondition){
		
		dSum += Poisson(iCriticalValue, dMean);
		if (dSum >= dAlpha){
			bSatisfyCondition = 0;
		}
		else{
			iCriticalValue++;
		}
	}

	return iCriticalValue;

}

///////////////////////////////////////////////// BINOMIAL RELATED FUNCTIONS ///////////////////////

/************ 
Checks if a value is infinity
************/
bool IsInf(double value){

	return (value >= DBL_MAX);
}

/************ 
Checks if a value is nan
************/
bool IsNan(double val){

	return val != val;
}

/************
Computes choose(d,p); combinations of d objects taken as groups of p
************/
double Choose(int d, int p){
	
	double iResult = 1;
	
	for(int i = 1; i <= p; i++){
		iResult *= (double)(d - (i-1)) / (double)i;
	}

	//if (IsInf(iResult))
	//	cerr << "Choose(" << d << "," << p << ")" << " is INFINITY" << endl;

	return iResult;
}

/************
Computes choose(n,k); a more robust version
************/
double Choose2(int n, int k){

	double dResult = 0;

	int i = 0;

	errno = 0;

	double dSum1 = 0;
	double dSum2 = 0;

	for(i = n; i >= n - (k - 1); i--){
		errno = 0;
		dSum1 += log((double)i);
		if ((errno == EDOM) || (errno == ERANGE)){
			cerr << "log(i) error message: " << strerror(errno) << endl;
			cerr << "i = " << i << endl;
			cerr << "k = " << k << " n = " << n << endl;
		}
		errno = 0;
		dSum2 += log((double)(n-i+1));
		if ((errno == EDOM) || (errno == ERANGE)){
			cerr << "log(n-i+1) error message: " << strerror(errno) << endl;
			cerr << "n-i+1 = " << n - i + 1 << endl;
			cerr << "k = " << k << " n = " << n << endl;
		}
	}
	
	dResult = dSum1 - dSum2;
	dResult = exp(dResult);
	if (errno == ERANGE){
		cerr << "exp error message: " << strerror(errno) << endl;
		cerr << "dSum1 = " << dSum1 << " dSum2 = " << dSum2 << endl;
		cerr << "n = " << n << " k = " << k << " dResult = " << dResult << endl;
	}

	dResult = ceil(dResult);

	return dResult;
}


/************ 
Compute Binomial p.d.f.
************/
bool bElapsed = false;
bool bDisplay = false;
double Binomial(int k, int n, double p){

	double dProb = 1;

	//dProb = Choose(n,k);
	//if (IsInf(dProb)){
	//	return dProb;
	//}
	//dProb *= pow(p, (double)k);
	//dProb *= pow(1-p, (double)(n-k));

	for(int i = 1; i <= k; i++){
		dProb *= (double)(n - (i-1)) * p / (double)i;
	}
	if (bElapsed && bDisplay) {
		cout << "!!!(0) p=" << p << ", n="<<n<<", k="<<k<< endl;
		cout << "!!!(1) " << dProb << endl;
	}

	double temp = pow(1-p, n-k);
	if (bElapsed && bDisplay)
		cout << "!!!(2) " << temp << endl;
	dProb *= temp;
	if (bElapsed && bDisplay)
		cout << "!!!(3) " << dProb << endl;

	return dProb;
}

/************ 
Compute Binomial p.d.f. in a different way
************/
double Binomial2(int k, int n, double p){

	//this should be the case; but probably because of numerical problems
	//p can be slightly below 0 or slightly above 1
	//assert((p >= 0) && (p <= 1));
	if (p < 0){
		p = 0;
	}
	if (p > 1){
		p = 1;
	}

	//deal with the extreme cases
	if (p == 0){
		if (k == 0)
			return 1;
		else
			return 0;
	}

	if (p == 1){
		if (k == n)
			return 1;
		else
			return 0;
	}
	
	double dProbab = 0;
	int i = 0;
	
	errno = 0;

	double dSum1 = 0;
	double dSum2 = 0;
	for(i = n; i >= n - (k - 1); i--){
		errno = 0;
		dSum1 += log((double)i);
		if ((errno == EDOM) || (errno == ERANGE)){
			cerr << "log(i) error message: " << strerror(errno) << endl;
			cerr << "i = " << i << endl;
			cerr << "k = " << k << " n = " << n << " p = " << p << endl;
		}
		errno = 0;
		dSum2 += log((double)(n-i+1));
		if ((errno == EDOM) || (errno == ERANGE)){
			cerr << "log(n-i+1) error message: " << strerror(errno) << endl;
			cerr << "n-i+1 = " << n - i + 1 << endl;
			cerr << "k = " << k << " n = " << n << " p = " << p << endl;
		}
	}
	
	errno = 0;

	//EDOM = when the argument of log is < 0;
	//ERANGE = when the argument of log is == 0
	double dTerm3 = k * log(p);
	if ((errno == EDOM) || (errno == ERANGE)){
		cerr << "log(p) error message: " << strerror(errno) << endl;
		cerr << "k = " << k << " n = " << n << " p = " << p << endl;
		dTerm3 = -DBL_MAX;
	}

	errno = 0;
	
	double dTerm4 = (n-k) * log(1-p);
	if ((errno == EDOM) || (errno == ERANGE)){
		cerr << "log(1-p) error message: " << strerror(errno) << endl;
		cerr << "k = " << k << " n = " << n << " p = " << p << endl;
		dTerm4 = -DBL_MAX;
	}

	errno = 0;
	
	//cerr.precision(25);
	dProbab = dSum1 - dSum2 + dTerm3 + dTerm4;
	//cerr << dSum1 << " - " << dSum2 << " + " << dTerm3 << " + " << dTerm4 << " = " << dProbab << endl;
	dProbab = exp(dProbab);
	if (errno == ERANGE){
		cerr << "exp error message: " << strerror(errno) << endl;
		cerr << "dSum1 = " << dSum1 << " dSum2 = " << dSum2 << " dTerm3 = " << dTerm3 << " dTerm4 = " << dTerm4 << endl;
		cerr << "k = " << k << " n = " << n << " p = " << p << " dProbab = " << dProbab << endl;
		dProbab = 0;
	}

	if ((dProbab < 0) || (dProbab > 1)){
		cerr << dSum1 << " - " << dSum2 << " + " << dTerm3 << " + " << dTerm4 << " = " << dProbab << endl;
		cerr << "k = " << k << " n = " << n << " p = " << p << endl;
	}
	assert(dProbab >= 0 && dProbab <= 1);
	
	return dProbab;
}


/************ 
Compute the critical value at right of a Binomial distribution with parameters n and p at significance level dAlpha
************/
int GetCriticalValueBinomialRight(int n, double p, double dAlpha){

	int iCriticalValue = 1; //always return a right critical value of 1 
	
	if (dAlpha <= 0){
		cerr << "#Statistical significance level <= 0" << endl;
		return iCriticalValue;
	}

	//just for numerical reasons
	if (p < 0){
		p = 0;
	}
	if (p > 1){
		p = 1;
	}

	//cerr.precision(25);
	int bSatisfyCondition = 0;
	double dSum = Binomial2(0,n,p);
	double dPrevSum = 0;
	while((!bSatisfyCondition) && (iCriticalValue <= n)){
		
		//if ((iCriticalValue == 1500) || (iCriticalValue == 2500) || (iCriticalValue == 3500))
		//	getchar();

		dPrevSum = dSum;
		double dBinom = Binomial2(iCriticalValue, n, p);
		dSum += dBinom;
		if (bElapsed && bDisplay)
			cout << dSum << "(" << dBinom << ")" << endl;
		//if ((dSum > 0.9) && (fabs(dSum - dPrevSum) < 1.0E-15)){
		//	bSatisfyCondition = 1;
			//cerr << "End computation: " << dSum << " " << dPrevSum << endl;
		//	continue;
		//}

		if (dSum >= (1.0 - dAlpha)) {
			bSatisfyCondition = 1;
		}
		else{
			iCriticalValue++;
		}		
	}
	if (bElapsed && bDisplay)
		cout << "-----" << endl;

	//dAlpha may be so small that 1 - dAlpha is so close to 1
	//I may end up with iCriticalValue == n+1
	if (iCriticalValue == (n+1)){
		cerr << "Critical value reached n" << endl;
		cerr << " n = " << n << " p = " << p << " dAlpha = " << dAlpha << endl;
	}

	return iCriticalValue;
}

/************ 
Compute the critical value at right of a Binomial distribution with parameters n and p at significance level dAlpha
Use the normal approximation when possible
************/
int GetCriticalValueBinomialRight2(int n, double p, double dAlpha){

	int iCriticalValue = 1; //always return a critical value of at least 1
	
	if (dAlpha <= 0){
		cerr << "Statistical significance level <= 0" << endl;
		return iCriticalValue;
	}

	int bNormalApprox = 0;
	if ((n * p >= 5) && (n * (1-p) >= 5)){
		bNormalApprox = 1;
	}

	if (bNormalApprox == 0){ //cannot use the normal approximation
		
		int bSatisfyCondition = 0;
		double dSum = Binomial(0, n, p);
		while((!bSatisfyCondition) && (iCriticalValue <= n)){
			
			double dBinom = Binomial(iCriticalValue, n, p);
			dSum += dBinom;
			if (dSum >= (1.0 - dAlpha)){
				bSatisfyCondition = 1;
			}
			else{
				iCriticalValue++;
			}
			
		}

		if (iCriticalValue == (n+1)){
			cerr << "Right critical value: Binomial approx: Critical value reached n" << endl;
			iCriticalValue = (int)ceil(n * p);
		}
		
		return iCriticalValue;
	}
	else{

		//Probab(x <= iCriticalValue) for Binomial(n,p) = Probab(x <= iCriticalValue + 0.5) for Normal(np,np(1-p))
		double mu = n * p;
		double sigma = sqrt(n*p*(1-p));

		iCriticalValue = (int)ceil(n * p); // I can start from here, since I am searching for the critical value at right
		
		int bSatisfyCondition = 0;
		while((!bSatisfyCondition) && (iCriticalValue <= n)){
		
			double zscore = (iCriticalValue + 0.5 - mu) / sigma;
				
			//dCDF = Probab(Z <= zscore), for Z ~ N(0,1)
			double dCDF = 0; 
			if (zscore > 0){
				double x = zscore / sqrt(2.0);
				dCDF = 0.5 + 0.5 * erf(x);
			}
			else{
				zscore = -zscore;
				double x = zscore / sqrt(2.0);
				dCDF = 1 - (0.5 + 0.5 * erf(x));
			}

			if (dCDF >= (1.0 - dAlpha)){
				bSatisfyCondition = 1;
			}
			else{
				iCriticalValue++;
			}
		}

		if (iCriticalValue == (n+1)){
			cerr << "Right critical value: Normal approx: Critical value reached n" << endl;
			iCriticalValue = (int)ceil(n * p);
		}
		
		return iCriticalValue;
	}
}


/************ 
Compute the critical value at left of a Binomial distribution with parameters n and p at significance level dAlpha
************/
int GetCriticalValueBinomialLeft(int n, double p, double dAlpha){

	int iCriticalValue = 0;
	
	if (dAlpha <= 0){
		cerr << "#Statistical significance level <= 0" << endl;
		iCriticalValue++;
		return iCriticalValue;
	}

	int bSatisfyCondition = 1;
	double dSum = 0;
	while((bSatisfyCondition) && (iCriticalValue <= n)){
		
		double dBinom = Binomial(iCriticalValue, n, p);
		dSum += dBinom;
		if (dSum >= dAlpha){
			bSatisfyCondition = 0;
		}
		else{
			iCriticalValue++;
		}
	}


	return iCriticalValue;
}

/************ 
Compute the critical value at left of a Binomial distribution with parameters n and p at significance level dAlpha
Use the normal approximation when possible
************/
int GetCriticalValueBinomialLeft2(int n, double p, double dAlpha){

	int iCriticalValue = 1; //always return a critical value of at least 1
	
	if (dAlpha <= 0){
		cerr << "Statistical significance level <= 0" << endl;
		return iCriticalValue;
	}

	int bNormalApprox = 0;
	if ((n * p >= 5) && (n * (1-p) >= 5)){
		bNormalApprox = 1;
	}

	if (bNormalApprox == 0){ //cannot use the normal approximation
	
		int bSatisfyCondition = 1;
		double dSum = Binomial(0, n, p);
		while((bSatisfyCondition) && (iCriticalValue <= n)){

			double dBinom = Binomial(iCriticalValue, n, p);
			dSum += dBinom;
			if (dSum >= dAlpha){
				bSatisfyCondition = 0;
			}
			else{
				iCriticalValue++;
			}
		}

		return iCriticalValue;
	}
	else{
	
		double mu = n * p;
		double sigma = sqrt(n*p*(1-p));
		int bSatisfyCondition = 1;
		while((bSatisfyCondition) && (iCriticalValue <= n)){
		
			double zscore = (iCriticalValue + 0.5 - mu) / sigma;
				
			//dCDF = Probab(Z <= zscore), for Z ~ N(0,1)
			double dCDF = 0; 
			if (zscore > 0){
				double x = zscore / sqrt(2.0);
				dCDF = 0.5 + 0.5 * erf(x);
			}
			else{
				zscore = -zscore;
				double x = zscore / sqrt(2.0);
				dCDF = 1 - (0.5 + 0.5 * erf(x));
			}

			if (dCDF >= dAlpha){
				bSatisfyCondition = 0;
			}
			else{
				iCriticalValue++;
			}	
		}

		return iCriticalValue;
	}
}


#ifndef __GNUG__
	/************ 
	error function
	************/
	double erf(double x){

		assert(x >= 0);
		double dResult = (1 / sqrt(PI)) * gamma(0.5, x*x);	
		return dResult;
	}
#endif //__GNUG__


/************ 
lower incomplete gammma function
************/
#define ITMAX 25
double gamma(double a, double x){
	
	double dResult = 0;
	if (x < 0 || a <= 0){
		cerr << "Invalid arguments in gamma function" << endl;
		return dResult;
	}
	
	double dSum = 1.0 / a;
	double dTerm = 1.0 / a;

	for(int i = 1; i <= ITMAX; i++){
		a = a + 1;
		dTerm *= (x / a);
		dSum += dTerm;
	}

	dResult = pow(x,a) * exp(-x) * dSum;

	return dResult;
}

///////////////////////////////////////////////// Hypergeometric RELATED FUNCTIONS ///////////////////////
/************ 
compute the hyper-geometric distribution
************/
double Hypergeometric(int k, int M, int d){

	double dProbab = 0;

	long int A = (d-1)*(d-2)/2;
	long int B = d*(d-1)/2;
	
	if (k > (d-1)){
		return dProbab;
	}
	if ((M-k) > A){
		return dProbab;
	}
	
	int i;
	long int j;

	errno = 0;

	double dSum1 = 0;
	double dSum2 = 0;
	for(i = d-1; i >= d-1-(k-1); i--){
		errno = 0;
		dSum1 += log((double)i);
		if ((errno == EDOM) || (errno == ERANGE)){
			cerr << "log(i) error message: " << strerror(errno) << endl;
			cerr << "i = " << i << endl;
			cerr << "k = " << k << " M = " << M << " d = " << d << endl;
		}
		errno = 0;
		dSum2 += log((double)(d-i));
		if ((errno == EDOM) || (errno == ERANGE)){
			cerr << "log(d-i) error message: " << strerror(errno) << endl;
			cerr << "i = " << i << endl;
			cerr << "k = " << k << " M = " << M << " d = " << d << endl;
		}
	}
	
	double dSum3 = 0;
	double dSum4 = 0;
	for(j = A; j >= A-(M-k-1); j--){
		errno = 0;
		dSum3 += log(double(j));
		if ((errno == EDOM) || (errno == ERANGE)){
			cerr << "log(j) error message: " << strerror(errno) << endl;
			cerr << "j = " << j << endl;
			cerr << "k = " << k << " M = " << M << " d = " << d << endl;
		}
		errno = 0;
		dSum4 += log((double)(A-j+1));
		if ((errno == EDOM) || (errno == ERANGE)){
			cerr << "log(A-j+1) error message: " << strerror(errno) << endl;
			cerr << "j = " << j << endl;
			cerr << "k = " << k << " M = " << M << " d = " << d << endl;
		}
	}
	
	double dSum5 = 0;
	double dSum6 = 0;
	for(j = B; j >= B-(M-1); j--){
		errno = 0;
		dSum5 += log(double(j));
		if ((errno == EDOM) || (errno == ERANGE)){
			cerr << "log(j) error message: " << strerror(errno) << endl;
			cerr << "j = " << j << endl;
			cerr << "k = " << k << " M = " << M << " d = " << d << endl;
		}
		errno = 0;
		dSum6 += log((double)(B-j+1));
		if ((errno == EDOM) || (errno == ERANGE)){
			cerr << "log(B-j+1) error message: " << strerror(errno) << endl;
			cerr << "j = " << j << endl;
			cerr << "k = " << k << " M = " << M << " d = " << d << endl;
		}
	}
	
	dProbab = dSum1 - dSum2 + dSum3 - dSum4 - dSum5 + dSum6;
	errno = 0;
	dProbab = exp(dProbab);
	if (errno == ERANGE){
		cerr << "exp error message: " << strerror(errno) << endl;
		cerr << "dSum1 = " << dSum1 << " dSum2 = " << dSum2 << "dSum3 = " << dSum3 << "dSum4 = " << dSum4 << "dSum5 = " << dSum5 << "dSum6 = " << dSum6 << endl;
		cerr << "k = " << k << " M = " << M << " d = " << d << " dProbab = " << dProbab << endl;
		dProbab = 0;
	}

	assert((dProbab >= 0) && (dProbab <= 1));

	return dProbab;

}

/************ 
compute the right critical value of the hyper-geometric distribution
************/
int GetCriticalValueHypergeometricRight(int M, int d , double dAlpha){

	int iCriticalValue = 1; //always return a right critical value of 1 
	
	if (dAlpha <= 0){
		cerr << "#Statistical significance level <= 0" << endl;
		return iCriticalValue;
	}

	int bSatisfyCondition = 0;
	//double dDenominator = Choose2(d*(d-1)/2, M);
	double dSum = Hypergeometric(0,M,d);
	while((!bSatisfyCondition) && (iCriticalValue <= M)){
		
		double dHyper = Hypergeometric(iCriticalValue,M,d);
		dSum += dHyper;
		
		if (dSum >= (1.0 - dAlpha)) {
			bSatisfyCondition = 1;
		}
		else{
			iCriticalValue++;
		}		
	}
	
	//dAlpha may be so small that 1 - dAlpha is so close to 1
	//I may end up with iCriticalValue == M+1
	if (iCriticalValue == (M+1)){
		cerr << "Critical value reached M" << endl;
		cerr << " M = " << M << " d = " << d << " dAlpha = " << dAlpha << endl;
	}


	return iCriticalValue;

}

///////////////////////////////////////////////// KOLMOGOROV RELATED FUNCTIONS ///////////////////////
/************
Computes the critical value of the Kolmogorov distribution given the signficance level dAlpha
t_star is the critical value for alpha = 0.01
************/
double GetCriticalValueKolmogorov(double t_star, double dAlpha){

	double dCriticalValue = t_star;
	
	double dConvergence = 1.0E-10;
	double dStep = 0.001;

	int bSatisfyCondition = 0;
	while(!bSatisfyCondition){
		
		double dSum = 0;
		double dPrevSum = -1;
		int i = 1;
		while (fabs(dSum - dPrevSum) > dConvergence){
		
			dPrevSum = dSum;
			dSum = dPrevSum + pow(-1, (double)i-1) * exp(-2 * i * i * dCriticalValue * dCriticalValue);
			i++;
		}
			
		double dVal = 1 - 2 * dSum;

		if (dVal >= 1 - dAlpha){
			bSatisfyCondition = 1;
		}
		else{
			dCriticalValue += dStep;
		}
	}

	return dCriticalValue;
}

/************
Computes the right critical value of standard Gaussian at significance level dAlpha
************/
double GetRightCriticalValueStandardNormal(double dAlpha){

	double dCriticalValue = 0;

	if (dAlpha == 0.001){
		dCriticalValue = 3.09;
		return dCriticalValue;
	}
	if (dAlpha == 0.005){
		dCriticalValue = 2.576;
		return dCriticalValue;
	}
	if (dAlpha == 0.01){
		dCriticalValue = 2.326;
		return dCriticalValue;
	}
	if (dAlpha == 0.025){
		dCriticalValue = 1.96;
		return dCriticalValue;
	}
	if (dAlpha == 0.05){
		dCriticalValue = 1.645;
		return dCriticalValue;
	}
	if (dAlpha == 0.1){
		dCriticalValue = 1.282;
		return dCriticalValue;
	}

	double dStart = 0;
	double dAreaSoFar = 0;
	if (dAlpha < 0.001){
		dStart = 3.09;
		dAreaSoFar = 1.0 - 0.001;
	}
	else{
		if (dAlpha < 0.005){
			dStart = 2.576;
			dAreaSoFar = 1.0 - 0.005;
		}
		else{
			if (dAlpha < 0.01){
				dStart = 2.326;
				dAreaSoFar = 1.0 - 0.01;
			}
			else{
				if (dAlpha < 0.025){
					dStart = 1.96;
					dAreaSoFar = 1.0 - 0.025;
				}
				else{
					if (dAlpha < 0.05){
						dStart = 1.645;
						dAreaSoFar = 1.0 - 0.05;
					}
					else{
						if (dAlpha < 0.1){
							dStart = 1.282;
							dAreaSoFar = 1.0 - 0.1;
						}
						else{
							dStart = 0;
							dAreaSoFar = 0.5;
						}
					}
				}
			}
		}
	}

	double dStep = 0.01;
	double dCurrentVal = dStart;
	int bSatisfyCondition = 0;
	while(!bSatisfyCondition){
		
		dCurrentVal += dStep;
		double dValToAdd = (1.0 / (sqrt(2*PI))) * exp(- (dCurrentVal * dCurrentVal) / 2.0);
		dAreaSoFar += dValToAdd;
		if (dAreaSoFar >= 1 - dAlpha){
			dCriticalValue = dCurrentVal;
			bSatisfyCondition = 1;
		}
	}
	
	return dCriticalValue;
}

////////////////////////////////////////////////	VARIOUS USEFUL FUNCTIONS	////////////////////////////////////////////

/************
Checks whether "elem" is an element of vector "v"
************/
int ExistElem(int elem, vector<int> &v){

    int bExist = 0;

	int iLengthVector = (int)v.size(); 
    for(int i = 0;i < iLengthVector && !bExist;i++)
        if (elem == v[i])
            bExist = 1;

    return bExist;
}

/************
Checks whether "elem" is an element of vector "v"
************/
int ExistElem(int elem, int* v, int no_elem){
	
	int bExist = 0;

	for(int i = 0;i < no_elem && !bExist;i++)
        if (elem == v[i])
            bExist = 1;

    return bExist;
}

/************ 
Check if vector v1 of length l1 equals vector v2 of length l2
************/
int AreVectorsEqual(int *v1, int l1, int *v2, int l2){

    if (l1 != l2)
        return 0;
    else{
        for(int i = 0; i < l1; i++){
            int bExist = ExistElem(v1[i], v2, l2);
			if (bExist == 0)
				return 0;
        }
    
    }

    return 1;
}

/************ 
Check if vector v1 equals vector v2
************/
int AreVectorsEqual(vector<int> &v1, vector<int> &v2){
	
	int l1 = (int)v1.size();
	int l2 = (int)v2.size();

    if (l1 != l2)
        return 0;
    else{
        for(int i = 0; i < l1; i++){
            int bExist = ExistElem(v1[i], v2);
			if (bExist == 0)
				return 0;
        }
    
    }

    return 1;
}

/************ 
Check if vector v1 is included in vector v2
************/
int CheckInclusionVectors(vector<int> &v1, vector<int> &v2){

	int l1 = (int)v1.size();
	int l2 = (int)v2.size();
	
	if (l1 > l2)
		return 0;
	else{
		for(int i = 0; i < l1; i++){
			int bExist = ExistElem(v1[i], v2);
			if (bExist == 0)
				return 0;
		}
	}
	
	return 1;
}

