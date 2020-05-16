// randomizer.h
// (C) 2013-2020 Nicholas G Davies

#ifndef RANDOMIZER_H
#define RANDOMIZER_H

#include <vector>
#include <limits>

#ifdef __WITH_MKL__
#include <mkl.h>
#include <cassert>
const MKL_INT SINGLE_VALUE = 1;
#define CHECK_ERROR() assert(err == VSL_STATUS_OK)
#else
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#endif

#include <cmath>
#include <stdexcept>
#include <string>

#include <Rcpp.h>

unsigned long num_multinomial = 0;


class Randomizer
{
public:
    Randomizer(unsigned long int seed = 0);
    ~Randomizer();

    void Reset();

    double Uniform(double min = 0.0, double max = 1.0);
    double RoundedUniform(double min = 0.0, double max = 1.0, double shoulder = 0.01);
    double Normal(double mean = 0.0, double sd = 1.0);
    double Normal(double mean, double sd, double clamp);
    double Cauchy(double x0 = 0.0, double gamma = 1.0);
    double LogNormal(double zeta = 0.0, double sd = 1.0);
    double Exponential(double rate = 1.0);
    double Gamma(double shape, double scale);
    double Beta(double alpha, double beta);
    unsigned int Discrete(unsigned int size);
    int Discrete(int min, int max);
    // int Discrete(std::vector<unsigned int>& cumulative_weights);
    // int Discrete(std::vector<double>& cumulative_weights);
    void Multinomial(unsigned int N, std::vector<double>& p, std::vector<unsigned int>& n);
    bool Bernoulli(double p);
    unsigned int Binomial(unsigned int n, double p);
    unsigned int BetaBinomial(unsigned int n, double p, double a_plus_b);
    int Poisson(double mean);
    int Geometric(double p);
    int Round(double x);

    template <typename T>
    void Shuffle(std::vector<T>& vec);

#ifdef WITH_MKL
    inline VSLStreamStatePtr MKL_RNG()
#else
    inline gsl_rng* GSL_RNG()
#endif
    {
      return r;
    }

private:
    unsigned long int seed;
#ifdef WITH_MKL
    VSLStreamStatePtr r;
    int err;
    double result[1];
#else
    gsl_rng* r;
#endif
};

Randomizer::Randomizer(unsigned long int s)
#ifdef WITH_MKL
 : seed(s)
#else
 : seed(s), r(gsl_rng_alloc(gsl_rng_mt19937))
#endif
{
#ifdef WITH_MKL
    err = vslNewStream(&r, VSL_BRNG_MT19937, s);
    CHECK_ERROR();
#endif
    Reset();
}

Randomizer::~Randomizer()
{
#ifdef WITH_MKL
    err = vslDeleteStream(&r);
    CHECK_ERROR();
#else
    gsl_rng_free(r);
#endif
}

void Randomizer::Reset()
{
#ifdef WITH_MKL
    err = vslDeleteStream(&r);
    CHECK_ERROR();
    err = vslNewStream(&r, VSL_BRNG_MT19937, s);
    CHECK_ERROR();
#else
    gsl_rng_set(r, seed);
#endif
}

double Randomizer::Uniform(double min, double max)
{
#ifdef WITH_MKL
    err = vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, r, SINGLE_VALUE, result, min, max);
    CHECK_ERROR();
    return result[0];
#else
    return min + gsl_rng_uniform(r) * (max - min);
#endif
}

double Randomizer::RoundedUniform(double min, double max, double shoulder)
{
    if (min >= max)
        return min;
    double z = Uniform();
    double sd = shoulder * (max - min) / ((1 - shoulder) * 2.50662827463);
    if (z < shoulder / 2)
        return min - abs(Normal(0, sd));
    else if (z < shoulder)
        return max + abs(Normal(0, sd));
    else
        return Uniform(min, max);
}

double Randomizer::Normal(double mean, double sd)
{
#ifdef WITH_MKL
    err = vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, r, SINGLE_VALUE, result, mean, sd);
    CHECK_ERROR();
    return result[0];
#else
    return mean + gsl_ran_gaussian_ziggurat(r, sd);
#endif
}

double Randomizer::Normal(double mean, double sd, double clamp)
{
    double n;
    do n = Normal(mean, sd); while (std::fabs(n - mean) > clamp);
    return n;
}

double Randomizer::LogNormal(double zeta, double sd)
{
#ifdef WITH_MKL
    err = vdRngLognormal(VSL_RNG_METHOD_LOGNORMAL_BOXMULLER2, r, SINGLE_VALUE, result, zeta, sd, 0, 1);
    CHECK_ERROR();
    return result[0];
#else
    return gsl_ran_lognormal(r, zeta, sd);
#endif
}

double Randomizer::Cauchy(double x0, double gamma)
{
#ifdef WITH_MKL
    err = vdRngCauchy(VSL_RNG_METHOD_CAUCHY_ICDF, r, SINGLE_VALUE, result, x0, gamma);
    CHECK_ERROR;
    return result[0];
#else
    return x0 + gsl_ran_cauchy(r, gamma);
#endif
}

double Randomizer::Exponential(double rate)
{
#ifdef WITH_MKL
    err = vdRngExponential(VSL_RNG_METHOD_EXPONENTIAL_ICDF, r, SINGLE_VALUE, result, 0, rate);
    CHECK_ERROR();
    return result[0];
#else
    return gsl_ran_exponential(r, 1. / rate);
#endif
}

double Randomizer::Gamma(double shape, double scale)
{
#ifdef WITH_MKL
    err = vdRngGamma(VSL_RNG_METHOD_GAMMA_GNORM, r, SINGLE_VALUE, result, shape, 0, scale);
    CHECK_ERROR();
    return result[0];
#else
    return gsl_ran_gamma(r, shape, scale);
#endif
}

double Randomizer::Beta(double alpha, double beta)
{
#ifdef WITH_MKL
    err = vdRngBeta(VSL_RNG_METHOD_BETA_CJA, r, SINGLE_VALUE, result, alpha, beta, 0, 1);
    CHECK_ERROR();
    return result[0];
#else
    return gsl_ran_beta(r, alpha, beta);
#endif
}

unsigned int Randomizer::Discrete(unsigned int size)
{
#ifdef WITH_MKL
    if (size > INT_MAX) {
      throw std::runtime_error("Generator only works with signed integers; this number is too large.");
    }
    err = viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, r, SINGLE_VALUE, result, 0, size);
    CHECK_ERROR();
    return result[0];
#else
    // TODO improve...
    if (size > gsl_rng_max(r))
        throw std::runtime_error("Generator cannot produce integers larger than " + std::to_string(gsl_rng_max(r)));
    return gsl_rng_get(r) % size;
#endif
}

int Randomizer::Discrete(int min, int max)
{
#ifdef WITH_MKL
    err = viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, r, SINGLE_VALUE, result, min, max);
    CHECK_ERROR();
    return result[0];
#else
    return min + gsl_rng_uniform_int(r, max - min + 1);
#endif
}
// ///
// int Randomizer::Discrete(std::vector<unsigned int>& cumulative_weights)
// {

//     if (cumulative_weights.back() == 0)
//         return Discrete(cumulative_weights.size());
//     return std::distance(cumulative_weights.begin(),
//                 std::lower_bound(cumulative_weights.begin(),
//                     cumulative_weights.end(),
//                     Discrete(cumulative_weights.back())));
// }
// ///
// int Randomizer::Discrete(std::vector<double>& cumulative_weights)
// {
//     if (cumulative_weights.back() == 0)
//         return Discrete(cumulative_weights.size());
//     return std::distance(cumulative_weights.begin(),
//                 std::lower_bound(cumulative_weights.begin(),
//                     cumulative_weights.end(),
//                     Uniform(0.0, cumulative_weights.back())));
// }

void Randomizer::Multinomial(unsigned int N, std::vector<double>& p, std::vector<unsigned int>& n)
{
  /* if (N > 5) {
    Rcpp::Rcout << "K: " << p.size() << std::endl;
    Rcpp::Rcout << "N: " << N << std::endl;
    Rcpp::Rcout << "p: " << std::endl;
    for (int i=0; i<p.size(); i++) {
        Rcpp::Rcout << p[i] << std::endl;
    }
    Rcpp::Rcout << "n: " << std::endl;
    for (int i=0; i<p.size(); i++) {
        Rcpp::Rcout << n[i] << std::endl;
    }
    Rcpp::Rcout << "====" << std::endl;
    } */
  num_multinomial += 1;
#ifdef WITH_MKL
    err = viRngMultinomial(VSL_RNG_METHOD_MULTINOMIAL_MULTPOISSON, r, SINGLE_VALUE, &n[0], N, p.size(), &p[0]);
    CHECK_ERROR();
#else
    gsl_ran_multinomial(r, p.size(), N, &p[0], &n[0]);
#endif
}

bool Randomizer::Bernoulli(double p)
{
#ifdef WITH_MKL
    err = viRngBernoulli(VSL_RNG_METHOD_BERNOULLI_ICDF, r, SINGLE_VALUE, result, p);
    CHECK_ERROR();
    return result[0];
#else
    if (p <= 0) return false;
    if (p >= 1) return true;
    return gsl_rng_uniform(r) < p;
#endif
}

unsigned int Randomizer::Binomial(unsigned int n, double p)
{
#ifdef WITH_MKL
    err = viRngBernoulli(VSL_RNG_METHOD_BINOMIAL_BTPE, r, SINGLE_VALUE, result, n, p);
    CHECK_ERROR();
    return result[0];
#else
    if (p <= 0) return 0;
    return gsl_ran_binomial(r, p, n);
#endif
}

unsigned int Randomizer::BetaBinomial(unsigned int n, double p, double a_plus_b)
{
    if (a_plus_b > 0)
       p = Beta(a_plus_b * p, a_plus_b * (1 - p));
    return Binomial(n, p);
}

int Randomizer::Poisson(double mean)
{
#ifdef WITH_MKL
    err = viRngPoisson(VSL_RNG_METHOD_POISSON_PTPE, r, SINGLE_VALUE, result, mean);
    CHECK_ERROR();
    return result[0];
#else
    if (mean <= 0) return 0;
    return gsl_ran_poisson(r, mean);
#endif
}

int Randomizer::Geometric(double p)
{
#ifdef WITH_MKL
    err = viRngGeometric(VSL_RNG_METHOD_GEOMETRIC_ICDF, r, SINGLE_VALUE, result, p);
    CHECK_ERROR();
    return result[0];
#else
    if (p <= 0) return 0;
    return gsl_ran_geometric(r, p);
#endif
}

int Randomizer::Round(double x)
{
    int sign = x < 0 ? -1 : 1;
    double intpart, fracpart;
    fracpart = std::modf(std::fabs(x), &intpart);
    return sign * (intpart + Bernoulli(fracpart));
}


#ifdef WITH_MKL

/* MKL doesn't provide an equivalent to gsl_ran_shuffle, so
   we borrow the implementation from randist/shuffle.c at
   https://github.com/CURTLab/gsl-mkl, but replace the
   GSL uniform call with an MKL one. */

static inline void 
copy (void * dest, size_t i, void * src, size_t j, size_t size)
{
  register char * a = size * i + (char *) dest ;
  register char * b = size * j + (char *) src ;
  register size_t s = size ;
  
  do                                            
    {                                           
      *a++ = *b++;                              
    } 
  while (--s > 0);                              
}

/* Randomly permute (shuffle) N indices
   Supply an array x[N] with nmemb members, each of size size and on
   return it will be shuffled into a random order.  The algorithm is
   from Knuth, SemiNumerical Algorithms, v2, p139, who cites Moses and
   Oakford, and Durstenfeld */

void
gsl_ran_shuffle (const gsl_rng * r, void * base, size_t n, size_t size)
{
  size_t i ;

  for (i = n - 1; i > 0; i--)
    {
      size_t j = Discrete(0, i);
      /* previously gsl_rng_uniform_int(r, i+1); */
      /* originally (i + 1) * gsl_rng_uniform (r) */

      swap (base, size, i, j) ;
    }
}
#endif


template <typename T>
void Randomizer::Shuffle(std::vector<T>& vec)
{
    gsl_ran_shuffle(r, &vec[0], vec.size(), sizeof(vec[0]));
}


#endif
