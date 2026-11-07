#include "common.h"
#include "cPopulation.h"
#include "option_parser.h"
#include "options.h"
#include "random.h"

using namespace bpopsim;
using namespace std;

int bpopsim_default_action(int argc, char* argv[])
{
  OptionParser options("Usage: bpopsim [options] --output-KIND [--output-KIND ...]");
  DefineBpopsimOptions(options);

  try {
    options.Parse(argc, argv);
  } catch (const std::invalid_argument& e) {
    options.AddUsage("");
    options.AddUsage(e.what());
    options.PrintUsage();
    return -1;
  }

  // Print help if requested.
  if (options.count("help")) {
    options.PrintUsage();
    return 0;
  }

  // Require one output option
  if (!options.count("output-average-fitness") && !options.count("output-average-mutations") &&
      !options.count("output-clade-frequencies") && !options.count("output-genotype-frequencies") &&
      !options.count("output-diverged-frequencies") && !options.count("output-longest-sweep-time") &&
      !options.count("output-muller")) {

    options.AddUsage("");
    options.AddUsage("You must provide at least one output option.");
    options.PrintUsage();
    return -1;
  }

  // Convert the command line into parameter structs; everything downstream uses these.
  SimulationParameters simulation_parameters = ParseSimulationParameters(options);
  OutputParameters output_parameters = ParseOutputParameters(options);
  MarkerDivergence marker_divergence = ParseMarkerDivergence(options);

  try {
    simulation_parameters.Validate();
  } catch (const std::invalid_argument& e) {
    options.AddUsage("");
    options.AddUsage(e.what());
    options.PrintUsage();
    return -1;
  }

  // Create output directory
  create_path(output_parameters.output_directory_name);

  // Use one main random number generator for entire program (don't reset for each replicate)
  auto rng_seed = from_string<unsigned long>(options["rng-seed"]);
  if (rng_seed == 0) rng_seed = RandomGenerator::RandomSeed();
  RandomGenerator rng(rng_seed);

  auto replicates = from_string<uint32_t>(options["replicates"]);

  simulation_parameters.Print();

  cStatistics final_statistics(output_parameters);

  cout << '\n' << "***** Beginning Simulations *****" << '\n';
  cout << "(Random number generator seed = " << rng_seed << ")" << '\n';

  for (uint32_t on_replicate = 1; on_replicate <= replicates; on_replicate++) {
    // Initialize Population object
    cPopulation population(simulation_parameters, output_parameters, marker_divergence, rng, on_replicate);

    // Perform Simulation
    population.RunSimulation();

    // Save the statistics
    final_statistics.push_back(population.replicate_statistics);

    cout << '\n' << "***** Replicate Output *****" << '\n';

    // Output Per-Simulation Files
    if (options.count("output-clade-frequencies")) {
      population.OutputCladeFrequencies();
    }

    if (options.count("output-genotype-frequencies")) {
      population.OutputGenotypeFrequencies();
    }

    if (options.count("output-muller")) {
      population.OutputMullerMatrix(from_string<uint32_t>(options["muller-resolution"]), options.count("muller-lod"));
    }
  }

  // Output Per-Execution Files

  cout << '\n' << "***** Summary Output *****" << '\n';

  if (options.count("output-average-fitness")) {
    final_statistics.OutputAveragePopulationFitness();
  }

  if (options.count("output-average-mutations")) {
    final_statistics.OutputAveragePopulationMutationCounts();
  }

  if (options.count("output-diverged-frequencies")) {
    final_statistics.OutputDivergedFrequenciesByDepth();
  }

  if (options.count("output-longest-sweep-time")) {
    final_statistics.OutputLongestSweepTime();
  }

  return 0;
}

int main(int argc, char* argv[])
{
  return bpopsim_default_action(argc, argv);
}
