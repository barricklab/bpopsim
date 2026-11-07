#include "random.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <random>

namespace bpopsim {

  struct RandomGenerator::Impl
  {
    Impl()
    : rng(gsl_rng_alloc(gsl_rng_taus2))
    {}
    ~Impl() { gsl_rng_free(rng); }
    gsl_rng* rng;
  };

  RandomGenerator::RandomGenerator(unsigned long seed)
  : impl_(new Impl)
  , seed_(seed)
  {
    gsl_rng_set(impl_->rng, seed);
  }

  RandomGenerator::~RandomGenerator() = default;

  void RandomGenerator::Seed(unsigned long seed)
  {
    seed_ = seed;
    gsl_rng_set(impl_->rng, seed);
  }

  double RandomGenerator::Uniform()
  {
    return gsl_rng_uniform(impl_->rng);
  }

  unsigned long RandomGenerator::UniformInt(unsigned long n)
  {
    return gsl_rng_uniform_int(impl_->rng, n);
  }

  unsigned int RandomGenerator::Binomial(double p, unsigned int n)
  {
    return gsl_ran_binomial(impl_->rng, p, n);
  }

  double RandomGenerator::Exponential(double mean)
  {
    return gsl_ran_exponential(impl_->rng, mean);
  }

  unsigned long RandomGenerator::RandomSeed()
  {
    std::random_device device;
    unsigned long seed = 0;
    while (seed == 0) {
      seed = device();
    }
    return seed;
  }

}  // namespace bpopsim
