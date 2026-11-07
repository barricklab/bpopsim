#include "doctest.h"
#include "common.h"
#include "cPopulation.h"
#include "options.h"
#include "random.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;
using namespace bpopsim;

namespace {

  // Build an OptionParser exactly as main() does and parse the given arguments.
  struct ParsedOptions
  {
    explicit ParsedOptions(const vector<string>& args)
    : options("test")
    {
      argv_storage.push_back("bpopsim");
      for (const string& a : args) argv_storage.push_back(a);
      for (string& a : argv_storage) argv.push_back(const_cast<char*>(a.c_str()));
      DefineBpopsimOptions(options);
      argv.push_back(nullptr);
      options.Parse(static_cast<int>(argv_storage.size()), argv.data());
    }

    vector<string> argv_storage;
    vector<char*> argv;
    OptionParser options;
  };

}  // namespace

// ---- Direct construction (no command line) ----

TEST_CASE("SimulationParameters defaults match the command-line defaults")
{
  SimulationParameters sp;
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
  // Derived fields are computed by the default constructor too.
  CHECK(sp.transfer_dilution_factor == doctest::Approx(pow(2.0, 6.64)));
  CHECK_NOTHROW(sp.Validate());
}

TEST_CASE("Finalize derives dilution factor and transfer probability")
{
  SimulationParameters sp;
  sp.generations_per_transfer = 6.64;
  sp.initial_population_size_after_transfer = 1e5;
  sp.Finalize();

  double expected_factor = pow(2.0, 6.64);
  CHECK(sp.transfer_dilution_factor == doctest::Approx(expected_factor));
  CHECK(sp.binomial_sampling_transfer_probability == doctest::Approx(1.0 / expected_factor));
  CHECK(sp.final_population_size_at_transfer == doctest::Approx(expected_factor * 1e5));

  SUBCASE("integer generations give an exact power of two")
  {
    sp.generations_per_transfer = 4;
    sp.initial_population_size_after_transfer = 1000;
    sp.Finalize();
    CHECK(sp.transfer_dilution_factor == doctest::Approx(16.0));
    CHECK(sp.final_population_size_at_transfer == doctest::Approx(16000.0));
  }
}

TEST_CASE("Finalize derives category totals and fractions, and can be re-run")
{
  SimulationParameters sp;
  sp.mutation_rates_per_division = {1e-7, 2e-7, 3e-7};
  sp.mutation_fitness_effects = {0.1, 0.2, 0.3};
  sp.Finalize();

  CHECK(sp.total_mutation_rate_per_division == doctest::Approx(6e-7));
  CHECK(sp.total_mutation_rate_exponential_mean == doctest::Approx(1.0 / 6e-7));
  REQUIRE(sp.fractional_chances_of_mutation_categories.size() == 3);
  CHECK(sp.fractional_chances_of_mutation_categories[0] == doctest::Approx(1.0 / 6.0));
  CHECK(sp.fractional_chances_of_mutation_categories[1] == doctest::Approx(2.0 / 6.0));
  CHECK(sp.fractional_chances_of_mutation_categories[2] == doctest::Approx(3.0 / 6.0));

  double sum = 0;
  for (double f : sp.fractional_chances_of_mutation_categories) sum += f;
  CHECK(sum == doctest::Approx(1.0));

  // Re-running Finalize after a change must not accumulate stale entries.
  sp.mutation_rates_per_division = {5e-8};
  sp.mutation_fitness_effects = {0.1};
  sp.Finalize();
  REQUIRE(sp.fractional_chances_of_mutation_categories.size() == 1);
  CHECK(sp.fractional_chances_of_mutation_categories[0] == doctest::Approx(1.0));
}

TEST_CASE("Validate rejects mismatched -u and -s counts")
{
  SimulationParameters sp;
  sp.mutation_rates_per_division = {1e-7, 2e-7};
  sp.mutation_fitness_effects = {0.1};
  sp.Finalize();
  CHECK_THROWS_AS(sp.Validate(), std::invalid_argument);

  sp.mutation_fitness_effects = {0.1, 0.2};
  CHECK_NOTHROW(sp.Validate());
}

TEST_CASE("OutputParameters and MarkerDivergence defaults")
{
  OutputParameters op;
  CHECK(op.output_directory_name == ".");
  CHECK(op.coarse_graining == 1);
  CHECK(op.burn_in == 0);
  CHECK(op.diverged_mutation_depth == 10);
  CHECK(op.minimum_output_frequency == doctest::Approx(0.01));
  CHECK_FALSE(op.verbose);
  CHECK_FALSE(op.debug);

  MarkerDivergence md;
  CHECK_FALSE(md.marker_states);
}

// ---- Parsing from a command line ----

TEST_CASE("ParseSimulationParameters with no arguments equals the defaults")
{
  ParsedOptions p({});
  SimulationParameters parsed = ParseSimulationParameters(p.options);
  SimulationParameters defaults;

  CHECK(parsed.maximum_number_of_transfers == defaults.maximum_number_of_transfers);
  CHECK(parsed.initial_population_size == doctest::Approx(defaults.initial_population_size));
  CHECK(parsed.generations_per_transfer == doctest::Approx(defaults.generations_per_transfer));
  CHECK(parsed.mutation_fitness_effect_model == defaults.mutation_fitness_effect_model);
  CHECK(parsed.mutation_rates_per_division == defaults.mutation_rates_per_division);
  CHECK(parsed.mutation_fitness_effects == defaults.mutation_fitness_effects);
  CHECK(parsed.transfer_dilution_factor == doctest::Approx(defaults.transfer_dilution_factor));
}

TEST_CASE("ParseSimulationParameters reads repeated -u/-s and the model flags")
{
  ParsedOptions p({"-T", "4",   "-N", "1000", "-p", "500", "-u",   "1E-7", "-u",   "3E-7", "-s", "0.1",
                   "-s", "0.2", "-m", "-f",   "e",  "-1",  "0.05", "-1",   "0.02", "-n",   "12"});
  SimulationParameters sp = ParseSimulationParameters(p.options);

  CHECK(sp.maximum_number_of_transfers == 12);
  CHECK(sp.initial_population_size == doctest::Approx(500));
  CHECK(sp.initial_population_size_after_transfer == doctest::Approx(1000));
  CHECK(sp.transfer_dilution_factor == doctest::Approx(16.0));
  REQUIRE(sp.mutation_rates_per_division.size() == 2);
  CHECK(sp.mutation_rates_per_division[1] == doctest::Approx(3e-7));
  REQUIRE(sp.mutation_fitness_effects.size() == 2);
  CHECK(sp.mutation_fitness_effects[1] == doctest::Approx(0.2));
  CHECK(sp.multiplicative_fitness_effects);
  CHECK(sp.mutation_fitness_effect_model == "e");
  REQUIRE(sp.first_mutation_fitness_effects.size() == 2);
  CHECK(sp.first_mutation_fitness_effects[0] == doctest::Approx(0.05));
  CHECK(sp.first_mutation_fitness_effects[1] == doctest::Approx(0.02));
  CHECK_NOTHROW(sp.Validate());
}

TEST_CASE("ParseSimulationParameters with mismatched -u/-s fails Validate")
{
  ParsedOptions p({"-u", "1E-7", "-u", "2E-7", "-s", "0.1"});
  SimulationParameters sp = ParseSimulationParameters(p.options);
  CHECK_THROWS_AS(sp.Validate(), std::invalid_argument);
}

TEST_CASE("ParseOutputParameters and ParseMarkerDivergence read their options")
{
  ParsedOptions p({"-o", "out_dir", "-c", "5", "-b", "3", "--diverged-mutation-depth", "4",
                   "--minimum-output-frequency", "0.05", "--output-diverged-frequencies", "--output-longest-sweep-time",
                   "-v", "--debug", "-k", "2"});
  OutputParameters op = ParseOutputParameters(p.options);
  CHECK(op.output_directory_name == "out_dir");
  CHECK(op.coarse_graining == 5);
  CHECK(op.burn_in == 3);
  CHECK(op.diverged_mutation_depth == 4);
  CHECK(op.minimum_output_frequency == doctest::Approx(0.05));
  CHECK(op.output_diverged_frequencies);
  CHECK(op.output_longest_sweep_time);
  CHECK(op.verbose);
  CHECK(op.debug);

  MarkerDivergence md = ParseMarkerDivergence(p.options);
  CHECK(md.marker_states == 2);
}

// ---- Behavior built on the parameters ----

TEST_CASE("DetermineMutationCategory samples categories in proportion to their rates")
{
  SimulationParameters sp;
  sp.mutation_rates_per_division = {1e-7, 3e-7};
  sp.mutation_fitness_effects = {0.1, 0.1};
  sp.Finalize();
  RandomGenerator rng(12345);

  const int draws = 100000;
  int counts[2] = {0, 0};
  for (int i = 0; i < draws; ++i) {
    uint32_t c = sp.DetermineMutationCategory(rng);
    REQUIRE(c < 2);
    counts[c]++;
  }
  // Expect 25% / 75%; allow a generous tolerance for a 1e5-draw sample.
  CHECK(static_cast<double>(counts[0]) / draws == doctest::Approx(0.25).epsilon(0.03));
  CHECK(static_cast<double>(counts[1]) / draws == doctest::Approx(0.75).epsilon(0.03));
}

TEST_CASE("cGenotype::AddOneMutation applies the fitness-effect models")
{
  RandomGenerator rng(7);
  SimulationParameters sp;
  sp.mutation_fitness_effects = {0.1};

  SUBCASE("additive")
  {
    sp.mutation_fitness_effect_model = "u";
    sp.Finalize();
    cGenotype g(sp, 0);
    CHECK(g.fitness == doctest::Approx(1.0));
    REQUIRE(g.AddOneMutation(sp, rng));
    CHECK(g.fitness == doctest::Approx(1.1));
    CHECK(g.total_mutation_count == 1);
    REQUIRE(g.mutation_counts.size() == 1);
    CHECK(g.mutation_counts[0] == 1);
  }

  SUBCASE("multiplicative")
  {
    sp.multiplicative_fitness_effects = true;
    sp.Finalize();
    cGenotype g(sp, 0);
    REQUIRE(g.AddOneMutation(sp, rng));
    REQUIRE(g.AddOneMutation(sp, rng));
    CHECK(g.fitness == doctest::Approx(1.1 * 1.1));
    CHECK(g.total_mutation_count == 2);
  }

  SUBCASE("one-time model refuses a second mutation in the same category")
  {
    sp.mutation_fitness_effect_model = "o";
    sp.Finalize();
    cGenotype g(sp, 0);
    REQUIRE(g.AddOneMutation(sp, rng));
    CHECK(g.fitness == doctest::Approx(1.1));
    CHECK_FALSE(g.AddOneMutation(sp, rng));
    CHECK(g.fitness == doctest::Approx(1.1));
    CHECK(g.total_mutation_count == 1);
  }

  SUBCASE("first-mutation effects override the category effect")
  {
    sp.first_mutation_fitness_effects = {0.5};
    sp.Finalize();
    cGenotype g(sp, 0);
    REQUIRE(g.AddOneMutation(sp, rng));
    CHECK(g.fitness == doctest::Approx(1.5));
    // Second mutation falls back to the category effect.
    REQUIRE(g.AddOneMutation(sp, rng));
    CHECK(g.fitness == doctest::Approx(1.6));
  }

  SUBCASE("child genotype inherits fitness and counts from its ancestor")
  {
    sp.Finalize();
    cGenotype parent(sp, 0);
    REQUIRE(parent.AddOneMutation(sp, rng));
    cGenotype child(sp, 1, parent);
    CHECK(child.fitness == doctest::Approx(parent.fitness));
    CHECK(child.total_mutation_count == 1);
    CHECK(child.unique_node_id == 1);
  }
}
