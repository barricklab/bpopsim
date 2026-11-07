#include "doctest.h"
#include "common.h"
#include "cPopulation.h"
#include "random.h"

#include <cmath>
#include <vector>

using namespace std;
using namespace bpopsim;

namespace {

// A small, fast population: 1000 cells, 3 generations per transfer (8x dilution),
// a high mutation rate so mutations actually happen in a few transfers.
SimulationParameters SmallPopulation(uint32_t transfers)
{
  SimulationParameters sp;
  sp.initial_population_size = 1000;
  sp.initial_population_size_after_transfer = 1000;
  sp.generations_per_transfer = 3;
  sp.mutation_rates_per_division = {1e-3};
  sp.mutation_fitness_effects = {0.1};
  sp.maximum_number_of_transfers = transfers;
  sp.Finalize();
  return sp;
}

OutputParameters DebugOutput()
{
  OutputParameters op;
  op.debug = true;   // turns on the internal population-size consistency asserts
  return op;
}

} // namespace

TEST_CASE("a seeded simulation runs to the requested number of transfers")
{
  const uint32_t transfers = 8;
  SimulationParameters sp = SmallPopulation(transfers);
  RandomGenerator r(2024);
  cPopulation pop(sp, DebugOutput(), MarkerDivergence(), r, 1);
  pop.RunSimulation();

  const cReplicateStatistics& s = pop.replicate_statistics;
  // One record at transfer 0 plus one per completed transfer.
  CHECK(s.average_population_fitness.size() == transfers + 1);
  CHECK(s.average_total_mutation_count.size() == transfers + 1);
  CHECK(s.genotype_frequencies.size() == transfers + 1);
  CHECK(s.clade_frequencies.size() == transfers + 1);
}

TEST_CASE("fitness starts at the initial fitness and never decreases under beneficial mutations")
{
  SimulationParameters sp = SmallPopulation(10);
  RandomGenerator r(99);
  cPopulation pop(sp, DebugOutput(), MarkerDivergence(), r, 1);
  pop.RunSimulation();

  const vector<double>& f = pop.replicate_statistics.average_population_fitness;
  REQUIRE(!f.empty());
  CHECK(f.front() == doctest::Approx(sp.initial_fitness));
  for (double v : f) CHECK(v >= sp.initial_fitness - 1e-12);
}

TEST_CASE("genotype frequencies sum to one at every recorded transfer")
{
  SimulationParameters sp = SmallPopulation(10);
  RandomGenerator r(31337);
  cPopulation pop(sp, DebugOutput(), MarkerDivergence(), r, 1);
  pop.RunSimulation();

  for (const GenotypeFrequencyMap& gf : pop.replicate_statistics.genotype_frequencies) {
    double sum = 0;
    for (const auto& kv : gf) sum += kv.second;
    CHECK(sum == doctest::Approx(1.0));
  }
}

TEST_CASE("the ancestral clade always has frequency one and mutations accumulate")
{
  SimulationParameters sp = SmallPopulation(10);
  RandomGenerator r(4242);
  cPopulation pop(sp, DebugOutput(), MarkerDivergence(), r, 1);
  pop.RunSimulation();

  const cReplicateStatistics& s = pop.replicate_statistics;
  for (GenotypeFrequencyMap cf : s.clade_frequencies) {
    // Genotype 0 is the ancestor; every cell descends from it.
    CHECK(cf.GetFrequency(0) == doctest::Approx(1.0));
  }
  // With a 1e-3 per-division rate and ~7000 divisions per transfer, mutations
  // certainly arose, so the final average mutation count is positive.
  CHECK(s.average_total_mutation_count.back() > 0);
  CHECK(s.total_mutations > 0);
}

TEST_CASE("the same seed gives identical statistics; different seeds differ")
{
  SimulationParameters sp = SmallPopulation(10);

  RandomGenerator a(555);
  cPopulation pa(sp, OutputParameters(), MarkerDivergence(), a, 1);
  pa.RunSimulation();

  RandomGenerator b(555);
  cPopulation pb(sp, OutputParameters(), MarkerDivergence(), b, 1);
  pb.RunSimulation();

  RandomGenerator c(556);
  cPopulation pc(sp, OutputParameters(), MarkerDivergence(), c, 1);
  pc.RunSimulation();

  CHECK(pa.replicate_statistics.average_population_fitness == pb.replicate_statistics.average_population_fitness);
  CHECK(pa.replicate_statistics.total_mutations == pb.replicate_statistics.total_mutations);
  CHECK(pa.replicate_statistics.average_population_fitness != pc.replicate_statistics.average_population_fitness);
}

TEST_CASE("average fitness end condition stops the run early")
{
  SimulationParameters sp = SmallPopulation(200);
  sp.average_fitness_end_condition = 1.05;
  sp.Finalize();
  RandomGenerator r(8);
  cPopulation pop(sp, OutputParameters(), MarkerDivergence(), r, 1);
  pop.RunSimulation();

  const vector<double>& f = pop.replicate_statistics.average_population_fitness;
  REQUIRE(f.size() >= 2);
  CHECK(f.size() < 201);
  CHECK(f.back() >= 1.05);
  // Every earlier record was below the threshold.
  for (size_t i = 0; i + 1 < f.size(); ++i) CHECK(f[i] < 1.05);
}

TEST_CASE("burn-in transfers are not recorded")
{
  SimulationParameters sp = SmallPopulation(6);
  OutputParameters op;
  op.burn_in = 3;
  RandomGenerator r(1);
  cPopulation pop(sp, op, MarkerDivergence(), r, 1);
  pop.RunSimulation();

  // Transfers are numbered from -burn_in and the run continues until the
  // transfer number reaches maximum_number_of_transfers, so burn-in transfers
  // are performed in addition to the requested ones. Records start at
  // transfer 0: with -n 6 -b 3 the run does 9 transfers and records 0..6.
  CHECK(pop.replicate_statistics.average_population_fitness.size() == 7);
}

TEST_CASE("coarse graining records every Nth transfer")
{
  SimulationParameters sp = SmallPopulation(10);
  OutputParameters op;
  op.coarse_graining = 5;
  RandomGenerator r(1);
  cPopulation pop(sp, op, MarkerDivergence(), r, 1);
  pop.RunSimulation();

  // Records at transfers 0, 5, 10.
  CHECK(pop.replicate_statistics.average_population_fitness.size() == 3);
}
