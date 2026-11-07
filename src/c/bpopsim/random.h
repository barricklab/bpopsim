#pragma once

#include <memory>

namespace bpopsim {

  //! The single source of randomness in a simulation.
  //!
  //! Wraps a GSL Tausworthe (taus2) generator. GSL is kept, rather than the
  //! C++ <random> library, because its generators and distribution samplers
  //! produce identical streams on every platform for a given seed and GSL
  //! version, which the regression tests depend on. The distribution methods
  //! map one-to-one onto the GSL calls the simulation has always used, so a
  //! seed reproduces output from before this wrapper existed.
  //!
  //! Only random.cc includes GSL headers.
  class RandomGenerator
  {
  public:
    //! Seeds the generator. A seed of 0 selects GSL's default seed for taus2;
    //! callers wanting a random seed should draw one with RandomSeed().
    explicit RandomGenerator(unsigned long seed = 0);
    ~RandomGenerator();

    RandomGenerator(const RandomGenerator&) = delete;
    RandomGenerator& operator=(const RandomGenerator&) = delete;

    void Seed(unsigned long seed);
    unsigned long seed() const { return seed_; }

    //! Uniform double in [0, 1).
    double Uniform();

    //! Uniform integer in [0, n). n must be > 0.
    unsigned long UniformInt(unsigned long n);

    //! Number of successes in n Bernoulli trials with success probability p.
    unsigned int Binomial(double p, unsigned int n);

    //! Exponentially distributed double with the given mean.
    double Exponential(double mean);

    //! A non-deterministic seed from the operating system, never 0.
    static unsigned long RandomSeed();

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
    unsigned long seed_;
  };

}  // namespace bpopsim
