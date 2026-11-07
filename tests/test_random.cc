#include "doctest.h"
#include "random.h"

#include <cmath>
#include <set>
#include <vector>

using namespace std;
using namespace bpopsim;

// Reference values recorded 2026-11-07 from GSL 2.8 taus2 with seed 12345,
// identical from conda-forge GSL (clang 21, arm64) and Homebrew GSL (Apple
// clang 17, arm64). This pins the cross-platform reproducibility that the
// golden regression test and every published seeded run depend on. If GSL
// changes a sampler in a future version, this test is what will tell you.
TEST_CASE("seed 12345 reproduces the recorded GSL 2.8 reference sequence")
{
  RandomGenerator r(12345);

  const double uniform[5] = {0.14079645113088191, 0.85450767702423036, 0.54992264253087342, 0.48398289736360312,
                             0.38425721903331578};
  for (double expected : uniform) CHECK(r.Uniform() == expected);

  const unsigned long uniform_int[5] = {639, 150, 918, 329, 136};
  for (unsigned long expected : uniform_int) CHECK(r.UniformInt(1000) == expected);

  const unsigned int binomial[5] = {1053, 1037, 1058, 981, 1001};
  for (unsigned int expected : binomial) CHECK(r.Binomial(0.01, 100000) == expected);

  const double exponential[5] = {2110867.7712601679, 4196622.1659500916, 1151990.7039830985, 983597.8671306466,
                                 421827.07449777075};
  // Exponential draws go through the C library log(), whose last bit differs
  // between glibc and Apple libm (seen in CI: the 5th value differed by 1 ulp on
  // Ubuntu). Uniform, UniformInt and Binomial are integer algorithms and are
  // exact everywhere. Allow a few ulps here.
  for (double expected : exponential) CHECK(r.Exponential(1e6) == doctest::Approx(expected).epsilon(1e-13));
}

TEST_CASE("the same seed reproduces the same stream; different seeds differ")
{
  RandomGenerator a(12345), b(12345), c(12346);
  for (int i = 0; i < 100; ++i) {
    double ua = a.Uniform();
    CHECK(ua == b.Uniform());
    CHECK(ua >= 0.0);
    CHECK(ua < 1.0);
  }
  bool any_differ = false;
  for (int i = 0; i < 100; ++i) {
    if (a.Uniform() != c.Uniform()) any_differ = true;
  }
  CHECK(any_differ);
}

TEST_CASE("Seed() restarts the stream")
{
  RandomGenerator r(7);
  double first = r.Uniform();
  r.Uniform();
  r.Uniform();
  r.Seed(7);
  CHECK(r.Uniform() == first);
  CHECK(r.seed() == 7);
}

TEST_CASE("UniformInt stays in range and hits every value")
{
  RandomGenerator r(1);
  set<unsigned long> seen;
  for (int i = 0; i < 10000; ++i) {
    unsigned long v = r.UniformInt(10);
    REQUIRE(v < 10);
    seen.insert(v);
  }
  CHECK(seen.size() == 10);
}

TEST_CASE("Binomial has the right mean and variance and edge cases")
{
  RandomGenerator r(42);
  const int draws = 100000;
  const unsigned int n = 1000;
  const double p = 0.1;
  double sum = 0, sumsq = 0;
  for (int i = 0; i < draws; ++i) {
    unsigned int k = r.Binomial(p, n);
    REQUIRE(k <= n);
    sum += k;
    sumsq += static_cast<double>(k) * k;
  }
  double mean = sum / draws;
  double var = sumsq / draws - mean * mean;
  CHECK(mean == doctest::Approx(n * p).epsilon(0.01));
  CHECK(var == doctest::Approx(n * p * (1 - p)).epsilon(0.05));

  CHECK(r.Binomial(0.0, 100) == 0);
  CHECK(r.Binomial(1.0, 100) == 100);
  CHECK(r.Binomial(0.5, 0) == 0);
}

TEST_CASE("Exponential has the requested mean")
{
  RandomGenerator r(99);
  const int draws = 200000;
  const double mean = 1.0 / 6e-7;
  double sum = 0;
  for (int i = 0; i < draws; ++i) {
    double x = r.Exponential(mean);
    REQUIRE(x >= 0.0);
    sum += x;
  }
  CHECK(sum / draws == doctest::Approx(mean).epsilon(0.01));
}

TEST_CASE("RandomSeed is never zero and varies")
{
  unsigned long a = RandomGenerator::RandomSeed();
  unsigned long b = RandomGenerator::RandomSeed();
  CHECK(a != 0);
  CHECK(b != 0);
  // Two consecutive OS draws colliding is astronomically unlikely.
  CHECK(a != b);
}
