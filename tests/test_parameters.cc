#include "doctest.h"
#include "common.h"
#include "cPopulation.h"
#include "options.h"

#include <cmath>
#include <memory>
#include <string>
#include <vector>

using namespace std;
using namespace bpopsim;

namespace {

// Build an AnyOption exactly as main() does and parse the given arguments.
// Until the parameter structs are decoupled from AnyOption (roadmap Phase 3)
// this is the only way to construct SimulationParameters.
struct ParsedOptions {
  explicit ParsedOptions(const vector<string>& args)
    : options("test")
  {
    argv_storage.push_back("bpopsim");
    for (const string& a : args) argv_storage.push_back(a);
    for (string& a : argv_storage) argv.push_back(const_cast<char*>(a.c_str()));
    DefineBpopsimOptions(options);
    options.processCommandArgs(static_cast<int>(argv.size()), argv.data());
  }

  vector<string> argv_storage;
  vector<char*> argv;
  AnyOption options;
};

struct Rng {
  Rng(unsigned long seed) : rng(gsl_rng_alloc(gsl_rng_taus2)) { gsl_rng_set(rng, seed); }
  ~Rng() { gsl_rng_free(rng); }
  gsl_rng* rng;
};

} // namespace

TEST_CASE("defaults are parsed when no arguments are given")
{
  ParsedOptions p({});
  Rng r(1);
  SimulationParameters sp(p.options, r.rng);

  CHECK(sp.maximum_number_of_transfers == 50);
  CHECK(sp.initial_population_size == doctest::Approx(5e6));
  CHECK(sp.initial_population_size_after_transfer == doctest::Approx(5e6));
  CHECK(sp.generations_per_transfer == doctest::Approx(6.64));
  CHECK(sp.initial_fitness == doctest::Approx(1.0));
  CHECK_FALSE(sp.multiplicative_fitness_effects);
  CHECK(sp.mutation_fitness_effect_model == "u");
  REQUIRE(sp.mutation_rates_per_division.size() == 1);
  CHECK(sp.mutation_rates_per_division[0] == doctest::Approx(1e-7));
  REQUIRE(sp.mutation_fitness_effects.size() == 1);
  CHECK(sp.mutation_fitness_effects[0] == doctest::Approx(0.1));
}

TEST_CASE("dilution factor and transfer probability derive from generations per transfer")
{
  ParsedOptions p({"-T", "6.64", "-N", "1E5"});
  Rng r(1);
  SimulationParameters sp(p.options, r.rng);

  double expected_factor = pow(2.0, 6.64);
  CHECK(sp.transfer_dilution_factor == doctest::Approx(expected_factor));
  CHECK(sp.binomial_sampling_transfer_probability == doctest::Approx(1.0 / expected_factor));
  CHECK(sp.final_population_size_at_transfer == doctest::Approx(expected_factor * 1e5));
}

TEST_CASE("integer generations per transfer give an exact power of two")
{
  ParsedOptions p({"-T", "4", "-N", "1000"});
  Rng r(1);
  SimulationParameters sp(p.options, r.rng);

  CHECK(sp.transfer_dilution_factor == doctest::Approx(16.0));
  CHECK(sp.final_population_size_at_transfer == doctest::Approx(16000.0));
}

TEST_CASE("repeated -u and -s define mutation categories")
{
  ParsedOptions p({"-u", "1E-7", "-u", "2E-7", "-u", "3E-7", "-s", "0.1", "-s", "0.2", "-s", "0.3"});
  Rng r(1);
  SimulationParameters sp(p.options, r.rng);

  REQUIRE(sp.mutation_rates_per_division.size() == 3);
  REQUIRE(sp.mutation_fitness_effects.size() == 3);
  CHECK(sp.total_mutation_rate_per_division == doctest::Approx(6e-7));
  CHECK(sp.total_mutation_rate_exponential_mean == doctest::Approx(1.0 / 6e-7));

  REQUIRE(sp.fractional_chances_of_mutation_categories.size() == 3);
  CHECK(sp.fractional_chances_of_mutation_categories[0] == doctest::Approx(1.0 / 6.0));
  CHECK(sp.fractional_chances_of_mutation_categories[1] == doctest::Approx(2.0 / 6.0));
  CHECK(sp.fractional_chances_of_mutation_categories[2] == doctest::Approx(3.0 / 6.0));

  double sum = 0;
  for (double f : sp.fractional_chances_of_mutation_categories) sum += f;
  CHECK(sum == doctest::Approx(1.0));
}

TEST_CASE("DetermineMutationCategory samples categories in proportion to their rates")
{
  ParsedOptions p({"-u", "1E-7", "-u", "3E-7", "-s", "0.1", "-s", "0.1"});
  Rng r(12345);
  SimulationParameters sp(p.options, r.rng);

  const int draws = 100000;
  int counts[2] = {0, 0};
  for (int i = 0; i < draws; ++i) {
    uint32_t c = sp.DetermineMutationCategory(r.rng);
    REQUIRE(c < 2);
    counts[c]++;
  }
  // Expect 25% / 75%; allow a generous tolerance for a 1e5-draw sample.
  CHECK(static_cast<double>(counts[0]) / draws == doctest::Approx(0.25).epsilon(0.03));
  CHECK(static_cast<double>(counts[1]) / draws == doctest::Approx(0.75).epsilon(0.03));
}

TEST_CASE("multiplicative fitness flag and first-mutation effects are parsed")
{
  ParsedOptions p({"-m", "-1", "0.05", "-1", "0.02", "-f", "e"});
  Rng r(1);
  SimulationParameters sp(p.options, r.rng);

  CHECK(sp.multiplicative_fitness_effects);
  CHECK(sp.mutation_fitness_effect_model == "e");
  REQUIRE(sp.first_mutation_fitness_effects.size() == 2);
  CHECK(sp.first_mutation_fitness_effects[0] == doctest::Approx(0.05));
  CHECK(sp.first_mutation_fitness_effects[1] == doctest::Approx(0.02));
}

TEST_CASE("cGenotype::AddOneMutation applies additive and multiplicative models")
{
  Rng r(7);

  SUBCASE("additive")
  {
    ParsedOptions p({"-s", "0.1", "-f", "u"});
    SimulationParameters sp(p.options, r.rng);
    cGenotype g(sp, 0);
    CHECK(g.fitness == doctest::Approx(1.0));
    REQUIRE(g.AddOneMutation(sp, r.rng));
    CHECK(g.fitness == doctest::Approx(1.1));
    CHECK(g.total_mutation_count == 1);
    REQUIRE(g.mutation_counts.size() == 1);
    CHECK(g.mutation_counts[0] == 1);
  }

  SUBCASE("multiplicative")
  {
    ParsedOptions p({"-m", "-s", "0.1", "-f", "u"});
    SimulationParameters sp(p.options, r.rng);
    cGenotype g(sp, 0);
    REQUIRE(g.AddOneMutation(sp, r.rng));
    REQUIRE(g.AddOneMutation(sp, r.rng));
    CHECK(g.fitness == doctest::Approx(1.1 * 1.1));
    CHECK(g.total_mutation_count == 2);
  }

  SUBCASE("one-time model refuses a second mutation in the same category")
  {
    ParsedOptions p({"-s", "0.1", "-f", "o"});
    SimulationParameters sp(p.options, r.rng);
    cGenotype g(sp, 0);
    REQUIRE(g.AddOneMutation(sp, r.rng));
    CHECK(g.fitness == doctest::Approx(1.1));
    CHECK_FALSE(g.AddOneMutation(sp, r.rng));
    CHECK(g.fitness == doctest::Approx(1.1));
    CHECK(g.total_mutation_count == 1);
  }

  SUBCASE("first-mutation effects override the category effect")
  {
    ParsedOptions p({"-s", "0.1", "-1", "0.5"});
    SimulationParameters sp(p.options, r.rng);
    cGenotype g(sp, 0);
    REQUIRE(g.AddOneMutation(sp, r.rng));
    CHECK(g.fitness == doctest::Approx(1.5));
    // Second mutation falls back to the category effect.
    REQUIRE(g.AddOneMutation(sp, r.rng));
    CHECK(g.fitness == doctest::Approx(1.6));
  }
}
