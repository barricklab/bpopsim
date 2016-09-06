#include "common.h"
#include "cPopulation.h"

using namespace bpopsim;
using namespace std;

bool g_verbose;

gsl_rng* initialize_random_number_generator(uint16_t& seed)
{
  const gsl_rng_type *T;
  gsl_rng * randgen;
  gsl_rng_env_setup();
  T = gsl_rng_taus2;
  randgen = gsl_rng_alloc(T);
  
  if (!seed)
    seed = time(NULL) * getpid();
  
  gsl_rng_set(randgen, seed);
  return randgen;
}

int bpopsim_default_action(int argc, char* argv[])
{
  
  // @JEB: doubles must be put in quotes if in scientific notation
  AnyOption options("Usage Example: bpopsim -T 6.64 -n 50 -i");
  options("help,h", "Produce this help message", TAKES_NO_ARGUMENT);
  options("verbose,v", "Output extra details of runs", TAKES_NO_ARGUMENT);
  options("debug", "Perform extra checks on simulation", TAKES_NO_ARGUMENT);

  options.addUsage("");
  options.addUsage("==== Simulation Options ====");
  options.addUsage("");
  options("rng-seed,d", "Seed for random number generator. [0 = time x pid]", 0);
  options("replicates,r", "Replicates", 1);
  options("initial-population-size,p", "Initial population size.", "5E6");
  options("generations-per-transfer,T", "Generations per transfer", 6.64);
  options("population-size-after-transfer,N", "Population size after transfer", "5E6");
  options("number-of-transfers,n", "Max number of transfer to replicate", 50);
  options("initial-fitness,z", "Initial fitness.", 1.0);
  options("mutation-rates,u", "Mutation rate per cell division. Supply option multiple times to define categories of mutations.", "1E-7");
  options("fitness-effects,s", "Mean fitness increment per mutation. Supply option multiple times to define categories of mutations.", 0.1);
  options("fitness-effect-model,f", "Distribution of mutation fitness effects. Choices are 'u' for uniform, 'e' for exponential distributions, 'o' for only one beneficial mutation possible in that category.", 'u');
  options("first-mutation-fitness-effects,1", "Use these fitness effects for the first mutations. Supply option multiple times to define the fitness effects of the 1st, 2nd, 3rd, etc. mutations. These values will override the -s option until the number of mutations exceeds the number of these values.");
  options("average-fitness-end-condition", "End a replicate when the average fitness of the population is this or greater (0=OFF)", 0);
  
  options.addUsage("");
  options.addUsage("==== Marker Divergence Options ====");
  options.addUsage("");
  options("marker-states,k", "Begin with the initial population divided equally into this many subpopulations each with a different state of a neutral genetic marker. Do not track genotypes. Only track the neutral marker state of each subpopulation. Default (0=OFF).", 0);
  options("max-marker-divergence-ratio,m", "Stop if divergence factor between the subpopulations with the most abundant and next-most abundant marker state exceeds this factor.", 100);
  options("transfer-interval-to-print,t", "Marker divergence printing intervals.", 1);
  options.addUsage("");
  options.addUsage("==== Statistical Output Files and Options ====");  
  options.addUsage("");  
  options("output-folder,o", "Base output folder. All output files will be created here.", ".");
  options("burn-in,b", "Perform this number of transfers before beginning to output statistics (before transfer numbered 0).", 0);
  options("coarse-graining,c", "Only print stats every this many transfers to all output files.", 1);
  options("output-average-fitness", "Output average fitness to file 'average_fitness.tab'",TAKES_NO_ARGUMENT);
  options("output-average-mutations", "Output mutation numbers to file 'average_mutations_count.tab'",TAKES_NO_ARGUMENT);
  options("output-clade-frequencies", "Output clade (mutation) frequencies (above minimum) for replicate X to file 'clade_frequencies_X.tab'", TAKES_NO_ARGUMENT);
  options("output-genotype-frequencies", "Output genotype frequencies (above minimum) for replicate X to file 'genotype_frequencies_X.tab'.", TAKES_NO_ARGUMENT);
  options("minimum-output-frequency", "Only output clade/genotype frequencies that are above this frequency at some time point.", 0.01);
  options("output-diverged-frequencies", "Output frequencies of clades that differ by at least X mutations from the line of descent to the final dominant to 'diverged_frequencies_X.tab'",TAKES_NO_ARGUMENT);
  options("diverged-mutation-depth", "Maximum depth of diverged mutations to output.", 10);
  options("output-largest-sweep-size", "Output the largest number of mutations sweeping simultaneously along the line of descent to the final dominant to 'largest_sweep_size.tab'",TAKES_NO_ARGUMENT);
  options("minimum-sweep-frequency-difference", "Adjacent mutations along the line of descent must have frequncies that are different by at least this threshold at one recorded time point to be considered part of different sweeps.",0.05);
  options("output-longest-sweep-time", "Output the longest time that a mutation along the line of descent to the final dominant took to sweep from ≥5% to ≥95% frequency in the population to 'longest_sweep_time.tab'",TAKES_NO_ARGUMENT);
  options.addUsage("");
  options.addUsage("==== Muller Plot Output Files and Options ====");  
  options.addUsage("");  
  options("output-muller", "Output Muller matrix for replicate X to file 'muller_matrix_X.dat'. Cell values unique genotypes. This file can be used to make plots of population dynamics.", TAKES_NO_ARGUMENT);  
  options("muller-resolution", "Muller matrix vertical resolution.", 200);
  options("muller-lod", "In the Muller matrix, give genotypes on the line of descent to the final dominant negative numbers, so that they can be plotted in separate colors.", TAKES_NO_ARGUMENT);

  // Not (re)implemented
//  options("output-dominant-genotypes", "Output dominant genotype frequencies to file 'dominant_genotypes.csv'", TAKES_NO_ARGUMENT);
//  ("output-time-to-sweep", "Print time to sweep for each mutation", TAKES_NO_ARGUMENT)
//  ("output-max-diff", "Print max difference of sweeping mutations", TAKES_NO_ARGUMENT)
//  ("output-single-fitness", "Print the fitness of the single cell.", TAKES_NO_ARGUMENT)
//  ("output-sweeping-descent-fitness", "Output the average fitness of sweeping descent per time.", TAKES_NO_ARGUMENT)
  
 // options.addUsage("");
 // options.addUsage("==== Other Options ====");
 // options.addUsage("");
 // options("exact-mutations-per-transfer,y", "Cause exactly this number of mutations each transfer --beneficial-mutation-rate will not be used. (0=OFF)", 0);
 // options("new-simulator", "Use the new simulator type.", TAKES_NO_ARGUMENT);
  
  options.processCommandArgs(argc, argv);
  
  // Print help if requested.
  if( options.count("help") ) {
		options.printUsage();
		return -1;
	}    
    
  // Require one output option
  if ( !options.count("output-average-fitness")
    && !options.count("output-average-mutations")
    && !options.count("output-clade-frequencies")
    && !options.count("output-genotype-frequencies")
    && !options.count("output-diverged-frequencies")
    && !options.count("output-largest-sweep-size")
    && !options.count("output-longest-sweep-time")
    && !options.count("output-muller")
  ) {
    
    options.addUsage("");
    options.addUsage("You must provide at least one output option.");
    options.printUsage();
		return -1;
  }
  
  // Create output directory
  create_path(options["output-folder"]);
  
  g_verbose = options.count("verbose");

  // Use one main random number generator for entire program (don't reset for each replicate)
  uint16_t rng_seed = from_string<int16_t>(options["rng-seed"]);
  gsl_rng* rng = initialize_random_number_generator(rng_seed);
        
  uint32_t replicates = from_string<uint32_t>(options["replicates"]);
    
  // This is just to display settings
  {
    cPopulation population(options, rng, 0);
    population.DisplayParameters();
  }
  
  cStatistics final_statistics(options);
  
  cout << endl << "***** Beginning Simulations *****" << endl;
  cout << "(Random number generator seed = " << rng_seed << ")" << endl;
  
  for (uint32_t on_replicate=1; on_replicate <= replicates; on_replicate++)
  {
    //Initialize Population object
    cPopulation population(options, rng, on_replicate);
    
    // Perform Simulation
    population.RunSimulation();
    
    // Save the statistics
    final_statistics.push_back(population.replicate_statistics);
    
    
    cout << endl << "***** Replicate Output *****" << endl;

    // Output Per-Simulation Files
    if( options.count("output-clade-frequencies") ) {
      population.OutputCladeFrequencies();
    }
    
    if( options.count("output-genotype-frequencies") ) {
      population.OutputGenotypeFrequencies();
    }
    
    if( options.count("output-muller") ) {
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
  
  cout << endl << "***** Summary Output *****" << endl;
  
  if( options.count("output-average-fitness")  ) {
    final_statistics.OutputAveragePopulationFitness();
  }
  
  if( options.count("output-average-mutations")  ) {
    final_statistics.OutputAveragePopulationMutationCounts();
  }
  
  if (options.count("output-diverged-frequencies") ) {
    final_statistics.OutputDivergedFrequenciesByDepth();
  }

  if (options.count("output-largest-sweep-size") ) {
    final_statistics.OutputLargestSweepSize();
  }
  
  if (options.count("output-longest-sweep-time") ) {
    final_statistics.OutputLongestSweepTime();
  }
  
  /*
  if (g_ro_only) {
    //Initialize Population object
    cPopulation population;
    cout << endl << "Printing to r/w ratio file.... " << endl;
    population.PrintOut_RedWhiteOnly(options["output-folder"], &red_white_ratios, from_string<uint16_t>(options["transfer-interval-to-print"]));
  }
  */
      

/*
      }
      
      while( population.GetTransfers() < population.GetTotalTransfers() ) {
        
        if (g_verbose) {
          for(vector<cSubpopulation>::iterator this_time = population.GetPopulation().begin(); this_time != population.GetPopulation().end(); this_time++) {
            cout << "Genotype3: " << this_time->GetNode_id() << " Frequency3: " << this_time->GetNumber() << endl;
          }
        }
        
        // Calculate the number of divisions until the next mutation 
        if( use_mute_num && number_of_mutations < from_string<uint32_t>(options["mut_num"]) && from_string<uint32_t>(options["mut_num"]) == 1) {
          population.SetDivisionsUntilMutation( on_run );
        }
        
        else if( use_mute_num && number_of_mutations >= from_string<uint32_t>(options["mut_num"]) && from_string<uint32_t>(options["mut_num"]) == 1 ) {
          population.SetDivisionsUntilMutation( pow(2, from_string<double>(options["generations-per-transfer"])) );
        }
        
        else if ( use_mute_num && from_string<uint32_t>(options["mut_num"]) != 1 ) {
          population.SetDivisionsUntilMutation( mutation_division[number_of_mutations+1] - mutation_division[number_of_mutations] );
          
          //population.SetDivisionsUntilMutation( pow(2, from_string<double>(options["generations-per-transfer"]))/pow(2, from_string<double>(options["generations-per-transfer"])/2) );
          //cout << population.GetPopulationSize() << " " << population.GetDivisionsUntilMutation() << " " << mutation_division.size() << " " << mutation_division[number_of_mutations] << " " << pow(2,from_string<double>(options["generations-per-transfer"])) << endl;
        }
        
        else if ( !use_mute_num ) {
          population.SetDivisionsUntilMutation(population.GetDivisionsUntilMutation() + round(gsl_ran_exponential(randgen, population.GetLambda())));
        }
        
        if (g_verbose) {
          for(vector<cSubpopulation>::iterator this_time = population.GetPopulation().begin(); this_time != population.GetPopulation().end(); this_time++) {
            cout << "Genotype4: " << this_time->GetNode_id() << " Frequency4: " << this_time->GetNumber() << endl;
          }
        }
        
        while( population.GetDivisionsUntilMutation() > 0 && (population.GetTransfers() < population.GetTotalTransfers())) 
        {
          
          if (g_verbose) {
            for(vector<cSubpopulation>::iterator this_time = population.GetPopulation().begin(); this_time != population.GetPopulation().end(); this_time++) {
              cout << "Genotype5: " << this_time->GetNode_id() << " Frequency5: " << this_time->GetNumber() << endl;
            }
          }
          
          //cout << population.GetTransfers() << " " << population.GetTotalTransfers() << " " << population.GetDivisionsUntilMutation() << endl;
          
          population.CalculateDivisions();
          
          if (g_verbose) {
            for(vector<cSubpopulation>::iterator this_time = population.GetPopulation().begin(); this_time != population.GetPopulation().end(); this_time++) {
              cout << "Genotype7: " << this_time->GetNode_id() << " Frequency7: " << this_time->GetNumber() << endl;
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
                cout << "Total mutations: " << population.GetTotalMutations() << " Maximum Fitness: " << population.GetMaxW() << endl;
                cout << "Size = " << current_ro_ratio.size() << endl;
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
  options
  ("help,h", "produce help message", TAKES_NO_ARGUMENT)
  ("verbose,v", "Verbose", TAKES_NO_ARGUMENT)  
  ("input-file,i", "Input file for converting external tree")
  ("output-folder,o", "Output folder", ".")
  ("muller-resolution,w", "Muller plot vertical resolution.", 2500)
  ;
  options.processCommandArgs(argc, argv);

  if (options.count("help") 
      || !options.count("input_file")
      || !options.count("output-folder")
      || !options.count("input-file")
      )
  {
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
  
  // Extract the sub-command argument
	string command;
	char* argv_new[argc];
	int argc_new = argc - 1;
  
  if (argc > 1) {
		command = argv[1];
		argv_new[0] = argv[0];
		for (int32_t i = 1; i < argc; i++)
			argv_new[i] = argv[i + 1];
	} else {
    return bpopsim_default_action(argc, argv); // Gives default usage in this case.
    return -1; 
	}
  
  // Pass the command to the proper handler
	command = to_upper(command);
  
  if (command == "CONVERT-TREE") {
    do_convert_tree(argc_new, argv_new);
  } else {
    bpopsim_default_action(argc, argv);
  }
}