#ifndef MyMath_h___
#define MyMath_h___
#include <stdlib.h>
#include <math.h>

#include <random>

extern std::mt19937 gen;

//static std::uniform_real_distribution<> dist(0, 1);
// -----------------------------------------------------------------------------
//static unsigned long long rand_hint = 0;

static inline float RUnif()
/* ========================================================================
* Returns a uniformly distributed real number.
* NOTE: generates a random number in the interval [0;1)
* Uses the standard routine from the ISO c++11 <random> library
* ========================================================================
*/
{
//return drand48();
return std::generate_canonical<float, 20>(gen);
}


static inline int RUnifInt(int max)
/* ========================================================================
* Returns a uniformly distributed integer number.
* NOTE: generates a random number in the interval [0;max)
* Uses the standard routine from the ISO c++11 <random> library
* ========================================================================
*/
{
std::uniform_int_distribution<int> dis(0, max-1);
return dis(gen);
}


static inline float RUnif( float lower, float upper )
/* ========================================================================
* Returns a uniformly distributed real number.
* NOTE: generates a random number in the interval [lower;upper)
* Uses the standard routine from the ISO c++11 <random> library
* ========================================================================
*/
{
	std::uniform_real_distribution<float> unifdist(lower,upper);
	return unifdist(gen);
}


static inline float RExp( float mean )
/* ========================================================================
* Returns a exponential distributed real number.
* NOTE: use mean > 0.0
*       the standard exp distribution uses factor lambda = 1/mean     
* Uses the standard routine from the ISO c++11 <random> library
* ========================================================================
*/
{
        std::exponential_distribution<float> expdist(1/mean);
	//return -log(RUnif())*mean;
	return expdist(gen);
}


static inline double RNormal(double m, double s)
/* ========================================================================
* Returns a normal (Gaussian) distributed real number.
* NOTE: use s > 0.0
*      Normal(m, s)     mean =  m    variance =   s*s
* Uses the standard routine from the ISO c++11 <random> library
* ========================================================================
*/
{
	std::normal_distribution<float> d(m,s);
	return d(gen);
}

static inline double RGamma(double a, double b)
/* ========================================================================
* Returns a Gamma distributed real number.
* NOTE: use a,b > 0.0
*      GammaDistribution(alpha, beta)     
* Uses the standard routine from the ISO c++11 <random> library
* ========================================================================
*/
{
	std::gamma_distribution<double> d(a,b);
	return d(gen);
}

// -----------------------------------------------------------------------------


static inline double Mylog(double x) {return(log(x));}
static inline double Myexp(double x) {return(exp(x));}
static inline double Mypow(double x, double p) {return( pow(x,p));}
static inline double Mysqrt(double x) {return( sqrt(x));}


// -----------------------------------------------------------------------------

static inline double MyMin(double x, double y)
{
	return x<y ? x : y;
}

// -----------------------------------------------------------------------------

static inline double MyMax(double x, double y)
{
	return x>y ? x : y;
}

// -----------------------------------------------------------------------------

static inline float DegToRad(float lat)
{
	return lat*M_PI/180.0;
}

// -----------------------------------------------------------------------------

static inline double MyAbs( double x )
{
	return x<0 ? -x : x;
}




// -----------------------------------------------------------------------------
// OBSOLETE
// double rNormal( double x, double mu, double sigma )
// {
// 	return 1./(sqrt(2.*M_PI)*sigma)*exp(-pow(x-mu,2)/(2.*pow(sigma,2)));
// }

// -----------------------------------------------------------------------------

static inline long MyRound( float x )
{
	return lround(x);
/*    if ( x-floor(x) >= 0.5 )
        return (int) ceil(x);
    else
        return (int) floor(x);*/
}

// -----------------------------------------------------------------------------

static inline double SlidingAverage( double average, double new_val, double old_val, int num )
{
	return ( (double)num*average+new_val-old_val )/(double)num;
}

// -----------------------------------------------------------------------------

static inline float MySignum( float x ) { return x>0 ? 1. : 0.; }

// -----------------------------------------------------------------------------




#endif



