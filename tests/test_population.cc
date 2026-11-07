#include "doctest.h"
#include "common.h"
#include "cPopulation.h"
#include "random.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <string>
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
    op.debug = true;  // turns on the internal population-size consistency asserts
    return op;
  }

}  // namespace

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

namespace {

  // Parse a frequency table written by OutputGenotypeFrequencies/OutputCladeFrequencies:
  // header "transfer\t<prefix>_<id>..." then one row per transfer. Returns id -> max frequency.
  map<uint32_t, double> ReadMaxFrequencyPerColumn(const string& path, const string& prefix)
  {
    ifstream in(path);
    REQUIRE(in.is_open());
    string line;
    vector<uint32_t> ids;
    while (getline(in, line)) {
      if (line.rfind("transfer", 0) == 0) {
        for (const string& field : split(line, "\t")) {
          if (field.rfind(prefix, 0) == 0) ids.push_back(from_string<uint32_t>(field.substr(prefix.size())));
        }
        break;
      }
    }
    map<uint32_t, double> max_frequency;
    for (uint32_t id : ids) max_frequency[id] = 0;
    while (getline(in, line)) {
      vector<string> fields = split(line, "\t");
      REQUIRE(fields.size() == ids.size() + 1);
      for (size_t i = 0; i < ids.size(); ++i) {
        max_frequency[ids[i]] = std::max(max_frequency[ids[i]], from_string<double>(fields[i + 1]));
      }
    }
    return max_frequency;
  }

}  // namespace

TEST_CASE("genotype frequency output selects columns by genotype, not clade, frequency")
{
  // Long enough, with a high mutation rate, that ancestors fall below the
  // threshold while their clades sweep: exactly the case where selecting by
  // clade frequency would add spurious columns.
  SimulationParameters sp = SmallPopulation(40);
  sp.mutation_rates_per_division = {5e-3};
  sp.Finalize();

  OutputParameters op;
  op.minimum_output_frequency = 0.05;
  std::filesystem::path dir = std::filesystem::temp_directory_path() / "bpopsim_test_genotype_columns";
  std::filesystem::remove_all(dir);
  std::filesystem::create_directories(dir);
  op.output_directory_name = dir.string();

  RandomGenerator r(777);
  cPopulation pop(sp, op, MarkerDivergence(), r, 1);
  pop.RunSimulation();
  pop.OutputGenotypeFrequencies();
  pop.OutputCladeFrequencies();

  // Expected column set: every genotype whose own frequency ever reached the threshold.
  set<uint32_t> expected;
  for (const GenotypeFrequencyMap& gf : pop.replicate_statistics.genotype_frequencies) {
    for (const auto& kv : gf) {
      if (kv.second >= op.minimum_output_frequency) expected.insert(kv.first);
    }
  }
  // And the clade-selected set, which must be a strict superset for this test to be meaningful.
  set<uint32_t> clade_selected;
  for (const GenotypeFrequencyMap& cf : pop.replicate_statistics.clade_frequencies) {
    for (const auto& kv : cf) {
      if (kv.second >= op.minimum_output_frequency) clade_selected.insert(kv.first);
    }
  }
  REQUIRE(!expected.empty());
  REQUIRE(clade_selected.size() > expected.size());

  map<uint32_t, double> written =
      ReadMaxFrequencyPerColumn((dir / "genotype_frequencies_0.0500_1.dat").string(), "genotype_");
  set<uint32_t> written_ids;
  for (const auto& kv : written) {
    written_ids.insert(kv.first);
    CHECK(kv.second >= op.minimum_output_frequency);  // no column that never met the threshold
  }
  CHECK(written_ids == expected);

  // The clade file keeps selecting by clade frequency.
  map<uint32_t, double> clade_written =
      ReadMaxFrequencyPerColumn((dir / "clade_frequencies_0.050000_1.dat").string(), "clade_");
  set<uint32_t> clade_written_ids;
  for (const auto& kv : clade_written) clade_written_ids.insert(kv.first);
  CHECK(clade_written_ids == clade_selected);

  std::filesystem::remove_all(dir);
}

// ---- Marker mode ----

namespace {

  MarkerDivergence Markers(uint32_t states, double end_fraction = 0)
  {
    MarkerDivergence md;
    md.marker_states = states;
    md.divergence_end_fraction = end_fraction;
    return md;
  }

}  // namespace

TEST_CASE("marker mode is off by default and records nothing")
{
  SimulationParameters sp = SmallPopulation(5);
  RandomGenerator r(1);
  cPopulation pop(sp, DebugOutput(), MarkerDivergence(), r, 1);
  pop.RunSimulation();
  CHECK(pop.replicate_statistics.marker_frequencies.empty());
}

TEST_CASE("marker mode seeds equal fractions and keeps them summing to one")
{
  for (uint32_t states : {2u, 3u, 5u}) {
    CAPTURE(states);
    SimulationParameters sp = SmallPopulation(8);
    sp.initial_population_size = 1500;  // divisible by 2, 3 and 5
    sp.initial_population_size_after_transfer = 1500;
    sp.Finalize();
    RandomGenerator r(11);
    cPopulation pop(sp, DebugOutput(), Markers(states), r, 1);
    pop.RunSimulation();

    const vector<vector<double>>& mf = pop.replicate_statistics.marker_frequencies;
    REQUIRE(mf.size() == pop.replicate_statistics.average_population_fitness.size());
    REQUIRE(mf.front().size() == states);
    for (double f : mf.front()) CHECK(f == doctest::Approx(1.0 / states));
    for (const vector<double>& at_transfer : mf) {
      REQUIRE(at_transfer.size() == states);
      double sum = 0;
      for (double f : at_transfer) {
        CHECK(f >= 0.0);
        CHECK(f <= 1.0);
        sum += f;
      }
      CHECK(sum == doctest::Approx(1.0));
    }
  }
}

TEST_CASE("each marker fraction equals the clade frequency of its root genotype")
{
  SimulationParameters sp = SmallPopulation(10);
  RandomGenerator r(99);
  cPopulation pop(sp, DebugOutput(), Markers(2), r, 1);
  pop.RunSimulation();

  const cReplicateStatistics& s = pop.replicate_statistics;
  REQUIRE(s.marker_frequencies.size() == s.clade_frequencies.size());
  for (size_t t = 0; t < s.marker_frequencies.size(); ++t) {
    GenotypeFrequencyMap clades = s.clade_frequencies[t];
    // Roots are genotype ids 0 and 1, in marker order.
    CHECK(s.marker_frequencies[t][0] == doctest::Approx(clades.GetFrequency(0)));
    CHECK(s.marker_frequencies[t][1] == doctest::Approx(clades.GetFrequency(1)));
  }
}

TEST_CASE("without mutations marker fractions start at one half and only drift")
{
  SimulationParameters sp = SmallPopulation(10);
  sp.mutation_rates_per_division = {0.0};
  sp.Finalize();
  RandomGenerator r(5);
  cPopulation pop(sp, OutputParameters(), Markers(2), r, 1);
  pop.RunSimulation();

  const vector<vector<double>>& mf = pop.replicate_statistics.marker_frequencies;
  REQUIRE(mf.size() == 11);
  CHECK(mf[0][0] == doctest::Approx(0.5));
  CHECK(mf[0][1] == doctest::Approx(0.5));
  for (const vector<double>& at_transfer : mf) {
    CHECK(at_transfer[0] > 0.0);
    CHECK(at_transfer[0] < 1.0);
  }
}

TEST_CASE("the divergence end fraction stops a replicate once one marker dominates")
{
  // Strong selection in a small population so one marker's lineage sweeps.
  SimulationParameters sp = SmallPopulation(400);
  sp.mutation_rates_per_division = {5e-3};
  sp.mutation_fitness_effects = {0.5};
  sp.Finalize();
  const double end_fraction = 0.9;

  RandomGenerator r(2024);
  cPopulation pop(sp, OutputParameters(), Markers(2, end_fraction), r, 1);
  pop.RunSimulation();

  const vector<vector<double>>& mf = pop.replicate_statistics.marker_frequencies;
  REQUIRE(mf.size() >= 2);
  CHECK(mf.size() < 401);  // ended early
  auto largest = [](const vector<double>& v) { return std::max(v[0], v[1]); };
  CHECK(largest(mf.back()) >= end_fraction);
  for (size_t t = 0; t + 1 < mf.size(); ++t) CHECK(largest(mf[t]) < end_fraction);

  // The same seed with the end fraction off runs the full length.
  RandomGenerator r2(2024);
  cPopulation full(sp, OutputParameters(), Markers(2, 0), r2, 1);
  full.RunSimulation();
  CHECK(full.replicate_statistics.marker_frequencies.size() == 401);
}

TEST_CASE("marker_frequencies.tab has one row per replicate and marker")
{
  SimulationParameters sp = SmallPopulation(6);
  OutputParameters op;
  std::filesystem::path dir = std::filesystem::temp_directory_path() / "bpopsim_test_marker_frequencies";
  std::filesystem::remove_all(dir);
  std::filesystem::create_directories(dir);
  op.output_directory_name = dir.string();

  cStatistics stats(op);
  RandomGenerator r(3);
  for (uint32_t replicate = 1; replicate <= 2; ++replicate) {
    cPopulation pop(sp, op, Markers(3), r, replicate);
    pop.RunSimulation();
    stats.push_back(pop.replicate_statistics);
  }
  stats.OutputMarkerFrequencies();

  ifstream in(dir / "marker_frequencies.tab");
  REQUIRE(in.is_open());
  string line;
  REQUIRE(getline(in, line));
  vector<string> header = split(line, "\t");
  REQUIRE(header.size() == 2 + 7);  // replicate, marker, transfers 0..6
  CHECK(header[0] == "replicate");
  CHECK(header[1] == "marker");
  CHECK(header[2] == "0");
  CHECK(header.back() == "6");

  int rows = 0;
  while (getline(in, line)) {
    vector<string> fields = split(line, "\t");
    REQUIRE(fields.size() == header.size());
    uint32_t replicate = from_string<uint32_t>(fields[0]);
    uint32_t marker = from_string<uint32_t>(fields[1]);
    REQUIRE(replicate >= 1);
    REQUIRE(replicate <= 2);
    REQUIRE(marker >= 1);
    REQUIRE(marker <= 3);
    for (size_t t = 0; t < 7; ++t) {
      CHECK(from_string<double>(fields[2 + t]) ==
            doctest::Approx(stats[replicate - 1].marker_frequencies[t][marker - 1]));
    }
    ++rows;
  }
  CHECK(rows == 2 * 3);
  std::filesystem::remove_all(dir);
}

TEST_CASE("a run with two groups of categories and mixed models records one count per category")
{
  SimulationParameters sp = SmallPopulation(10);
  sp.mutation_rates_per_division = {1e-3, 1e-4, 5e-4};
  sp.mutation_fitness_effects = {0.1, 0.3, 0.05};
  sp.mutation_fitness_effect_models = {'u', 'u', 'e'};
  sp.mutation_category_group_sizes = {2, 1};
  sp.Finalize();
  CHECK_NOTHROW(sp.Validate());

  RandomGenerator r(321);
  cPopulation pop(sp, DebugOutput(), MarkerDivergence(), r, 1);
  pop.RunSimulation();

  const cReplicateStatistics& s = pop.replicate_statistics;
  REQUIRE(!s.average_mutation_counts.empty());
  for (const vector<double>& per_category : s.average_mutation_counts) REQUIRE(per_category.size() == 3);
  for (const GenotypeFrequencyMap& gf : s.genotype_frequencies) {
    double sum = 0;
    for (const auto& kv : gf) sum += kv.second;
    CHECK(sum == doctest::Approx(1.0));
  }
  CHECK(s.total_mutations > 0);
}
