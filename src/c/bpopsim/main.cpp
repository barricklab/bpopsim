#include "common.h"
#include "cPopulation.h"
#include "options.h"
#include "random.h"

using namespace bpopsim;
using namespace std;

int bpopsim_default_action(int argc, char* argv[])
{
  AnyOption options("Usage Example: bpopsim -T 6.64 -n 50 -i");
  DefineBpopsimOptions(options);
  options.processCommandArgs(argc, argv);

  // Print help if requested.
  if (options.count("help")) {
    options.printUsage();
    return -1;
  }

  // Require one output option
  if (!options.count("output-average-fitness") && !options.count("output-average-mutations") &&
      !options.count("output-clade-frequencies") && !options.count("output-genotype-frequencies") &&
      !options.count("output-diverged-frequencies") && !options.count("output-largest-sweep-size") &&
      !options.count("output-longest-sweep-time") && !options.count("output-muller")) {

    options.addUsage("");
    options.addUsage("You must provide at least one output option.");
    options.printUsage();
    return -1;
  }

  // Convert the command line into parameter structs; everything downstream uses these.
  SimulationParameters simulation_parameters = ParseSimulationParameters(options);
  OutputParameters output_parameters = ParseOutputParameters(options);
  MarkerDivergence marker_divergence = ParseMarkerDivergence(options);

  try {
    simulation_parameters.Validate();
  } catch (const std::invalid_argument& e) {
    options.addUsage("");
    options.addUsage(e.what());
    options.printUsage();
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

    /*
    if( use_mute_num ) {
      cout << "Printing expectation value for set number of mutations..." << endl;
      population.PrintExpectationValue(options["output-folder"]);
    }
    */

    /*
    if( print_max_diff ) {
      cout << "Printing max difference of relevant mutations..." << endl;
      population.CalculateSimilarity(options["output-folder"]);
    }
     */

    /*
    if( print_time_to_sweep ) {
      cout << "Printing time to sweep..." << endl;
      population.TimeToSweep(options["output-folder"]);
    }
     */

    /*
    if( print_sweeping_descent_fitness ) {
      cout << "Printing fitness of the winning line of descent..." << endl;
      population.PrintWinningFitness(options["output-folder"], on_run);
    }
     */
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

  if (options.count("output-largest-sweep-size")) {
    final_statistics.OutputLargestSweepSize();
  }

  if (options.count("output-longest-sweep-time")) {
    final_statistics.OutputLongestSweepTime();
  }

  /*
  if (g_ro_only) {
    //Initialize Population object
    cPopulation population;
    cout << endl << "Printing to r/w ratio file.... " << endl;
    population.PrintOut_RedWhiteOnly(options["output-folder"], &red_white_ratios,
  from_string<uint16_t>(options["transfer-interval-to-print"]));
  }
  */

  /*
        }

        while( population.GetTransfers() < population.GetTotalTransfers() ) {

          if (g_verbose) {
            for(vector<cSubpopulation>::iterator this_time = population.GetPopulation().begin(); this_time !=
    population.GetPopulation().end(); this_time++) { cout << "Genotype3: " << this_time->GetNode_id() << " Frequency3: "
    << this_time->GetNumber() << endl;
            }
          }

          // Calculate the number of divisions until the next mutation
          if( use_mute_num && number_of_mutations < from_string<uint32_t>(options["mut_num"]) &&
    from_string<uint32_t>(options["mut_num"]) == 1) { population.SetDivisionsUntilMutation( on_run );
          }

          else if( use_mute_num && number_of_mutations >= from_string<uint32_t>(options["mut_num"]) &&
    from_string<uint32_t>(options["mut_num"]) == 1 ) { population.SetDivisionsUntilMutation( pow(2,
    from_string<double>(options["generations-per-transfer"])) );
          }

          else if ( use_mute_num && from_string<uint32_t>(options["mut_num"]) != 1 ) {
            population.SetDivisionsUntilMutation( mutation_division[number_of_mutations+1] -
    mutation_division[number_of_mutations] );

            //population.SetDivisionsUntilMutation( pow(2,
    from_string<double>(options["generations-per-transfer"]))/pow(2,
    from_string<double>(options["generations-per-transfer"])/2) );
            //cout << population.GetPopulationSize() << " " << population.GetDivisionsUntilMutation() << " " <<
    mutation_division.size() << " " << mutation_division[number_of_mutations] << " " <<
    pow(2,from_string<double>(options["generations-per-transfer"])) << endl;
          }

          else if ( !use_mute_num ) {
            population.SetDivisionsUntilMutation(population.GetDivisionsUntilMutation() +
    round(gsl_ran_exponential(randgen, population.GetLambda())));
          }

          if (g_verbose) {
            for(vector<cSubpopulation>::iterator this_time = population.GetPopulation().begin(); this_time !=
    population.GetPopulation().end(); this_time++) { cout << "Genotype4: " << this_time->GetNode_id() << " Frequency4: "
    << this_time->GetNumber() << endl;
            }
          }

          while( population.GetDivisionsUntilMutation() > 0 && (population.GetTransfers() <
    population.GetTotalTransfers()))
          {

            if (g_verbose) {
              for(vector<cSubpopulation>::iterator this_time = population.GetPopulation().begin(); this_time !=
    population.GetPopulation().end(); this_time++) { cout << "Genotype5: " << this_time->GetNode_id() << " Frequency5: "
    << this_time->GetNumber() << endl;
              }
            }

            //cout << population.GetTransfers() << " " << population.GetTotalTransfers() << " " <<
    population.GetDivisionsUntilMutation() << endl;

            population.CalculateDivisions();

            if (g_verbose) {
              for(vector<cSubpopulation>::iterator this_time = population.GetPopulation().begin(); this_time !=
    population.GetPopulation().end(); this_time++) { cout << "Genotype7: " << this_time->GetNode_id() << " Frequency7: "
    << this_time->GetNumber() << endl;
              }
            }

            if( population.GetDivisionsUntilMutation() <= 0 ) {
              if( !use_mute_num ) population.Mutate();

              else if( use_mute_num && number_of_mutations <  from_string<uint32_t>(options["mut_num"]) ) {
                population.Mutate();
                number_of_mutations++;
                cout << "Number of mutations: " << number_of_mutations << endl;
              }

            }

            if( population.GetPopulationSize() >= population.GetPopSizeBeforeDilution()) {

              population.FrequenciesPerTransferPerNode();
              population.CalculateAverageFitness();

              if( print_single_fit )
                population.Deterministic_Resample();
              else {
                population.Resample();
                population.CullPopulations();
              }

              if( print_single_fit && !use_mute_num) {
                population.PrintSingleFitness(options["output-folder"]);
                //cout << "Population size: " << population.GetPopulationSize() << endl;
              }

              if( g_verbose ) cout << "Total pop size3: " << population.GetPopulationSize() <<endl;

              count++;
              cout << "Passing.... " << count << endl;

              if ( population.GetTransfers() % from_string<uint16_t>(options["transfer-interval-to-print"]) == 0 ) {
                current_ro_ratio.push_back(population.GetRatio());

                if (g_verbose)
                {
                  cout << "Transfer " << population.GetTransfers() << " : " <<
                  "=>" << population.GetPopulationSize() << "  R/W Ratio: " << population.GetRatio() << endl;
                  cout << "Total mutations: " << population.GetTotalMutations() << " Maximum Fitness: " <<
    population.GetMaxW() << endl; cout << "Size = " << current_ro_ratio.size() << endl;
                }
              }
            }
            //if (g_verbose) population.PrintTree();
          }
        }

        //cout << endl << endl;
        //population.RunSummary();

        if (g_verbose) {
          population.PrintTree();
        }

        if (g_ro_only) {
          red_white_ratios.push_back(current_ro_ratio);
        }

        //cout << endl << endl;
        //population.PrintFreqsQuick();
        //population.PrintTree();

      }



    }

  */
  return 0;
}

int do_convert_tree(int argc, char* argv[])
{
  AnyOption options("Usage: bpopsim CONVERT-TREE -i ");
  options("help,h", "produce help message", TAKES_NO_ARGUMENT)("verbose,v", "Verbose", TAKES_NO_ARGUMENT)(
      "input-file,i", "Input file for converting external tree")("output-folder,o", "Output folder", ".")(
      "muller-resolution,w", "Muller plot vertical resolution.", 2500);
  options.processCommandArgs(argc, argv);

  if (options.count("help") || !options.count("input_file") || !options.count("output-folder") ||
      !options.count("input-file")) {
    options.printUsage();
    return -1;
  }

  /*
  cPopulation access_to_functions();
  access_to_functions.SetMullerRez(from_string<uint32_t>(options["muller-resolution"]));
  access_to_functions.ConvertExternalData(options["input_file"]);

  cerr << endl << "Generating Muller Matrix.... " << endl;
  vector< vector<int> > muller_matrix;
  access_to_functions.DrawMullerMatrix(options["output-folder"], muller_matrix);
  */
  return 0;
}

int main(int argc, char* argv[])
{

  if (argc <= 1) {
    return bpopsim_default_action(argc, argv);  // Gives default usage in this case.
  }

  // Extract the sub-command argument; argv_new is argv without it (and stays
  // NULL-terminated like a real argv).
  string command = to_upper(argv[1]);
  vector<char*> argv_new;
  argv_new.push_back(argv[0]);
  for (int i = 2; i < argc; i++) argv_new.push_back(argv[i]);
  argv_new.push_back(nullptr);
  int argc_new = argc - 1;

  // Pass the command to the proper handler
  if (command == "CONVERT-TREE") {
    return do_convert_tree(argc_new, argv_new.data());
  }
  return bpopsim_default_action(argc, argv);
}
