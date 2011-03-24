#include "cPopulation.h"

// setup and parse configuration options:
void get_cmdline_options(variables_map &options, int argc, char* argv[]) {

  options_description cmdline_options("Allowed options");
  cmdline_options.add_options()
  ("help,h", "produce this help message")
  ("generations-per-transfer,T", value<double>(), "Generations per transfer")
  ("population-size-after-transfer,N", value<uint64_t>(), "Population size after transfer")
  ("number-of-transfers,n", value<int>(), "Max number of transfer to replicate")
  ("output-file,o", value<std::string>(), "Output file")
  ("mutation-rate-per-division,u", value<double>(), "Mutation rate per division")
  ("average-selection-coefficient,s", value<double>(), "Average selection coefficient")
  ("time-interval,i", value<int>(), "Time interval")
  ("replicates,r", value<int>(), "Replicates")
  ("marker-divergence,m", value<int>(), "Max divergence factor")
  ("type-of-mutations,f", value<char>(), "Type of mutations")
  ("verbose,v", value<int>(), "Verbose")
  ("lineage-tree,l", value<int>(), "Lineage Tree")
  ("seed,d", value<long>(), "Seed for random number generator")
  ;

/* Need to add these as options...
  'fitnesses|f=s' => \@fitnesses,
  'detailed' => \$detailed,
  'minimum-data-points|z=s' => \$minimum_data_points,
  'maximum-data-points|x=s' => \$maximum_data_points,
  'skip-generations|k=s' => \$skip_generations,
  'input_initial_w|0=s' => \$input_initial_w,
  'drop_frequency|2=s' => \$drop_frequency,
  'multiplicative' => \$multiplicative_selection_coefficients
*/

  store(parse_command_line(argc, argv, cmdline_options), options);
  notify(options);

  // check here for required options
  if(options.count("help")
      || !options.count("output-file")) {
      std::cerr << "Usage: GSL_RNG_SEED=123 bpopsim -o output" << std::endl << std::endl;
      std::cerr << cmdline_options << std::endl;
      exit(0);
  }
}

int main(int argc, char* argv[])
{
	 tree<cGenotype>::iterator_base loc;
	 u_int64_t node_id;
	 long seed;
	
   //set up command line options
   variables_map cmdline_options;
   get_cmdline_options(cmdline_options, argc, argv);
   std::string output_file = cmdline_options["output-file"].as<std::string>();
	
   //Initialize Population object
	 cPopulation population;
   population.SetParameters(cmdline_options);
   population.DisplayParameters();
	
	 //create generator and seed
	 //@agm defaults to system time seed
	 const gsl_rng_type *T;
	 gsl_rng * randgen;
	 T = gsl_rng_mt19937;
	 randgen = gsl_rng_alloc (T);
	 if ( population.GetSeed() == 0 ) { seed = time (NULL) * getpid(); }
	 else { seed = population.GetSeed(); }
	 gsl_rng_set(randgen, seed);
	
	 //Initialize Tree object 
	 cLineageTree newtree;
	
	 std::vector< std::vector<cGenotypeFrequency> > frequencies;
	 
   for (u_int64_t on_run=0; on_run < population.GetReplicates(); on_run++)
   {
      population.ClearRuns(newtree);
		 
      u_int64_t count(0);
      std::cout << "Replicate " << on_run+1;   
		 
      population.NewSeedSubpopulation(newtree, node_id);
		  //std::cout << node_id << std::endl;
		 
      population.ResetRunStats();
		 
      while( (population.GetTransfers() < population.GetTotalTransfers()) && 
              population.GetKeepTransferring() ) 
			{
				
         // Calculate the number of divisions until the next mutation 
				 population.SetDivisionsUntilMutation(population.GetDivisionsUntilMutation() + round(gsl_ran_exponential(randgen, population.GetLambda())));
				
				 // 
				if (population.GetVerbose()) { 
					std::cout << "  New divisions before next mutation: " << population.GetDivisionsUntilMutation() << std::endl; }
         
			 	 while( population.GetDivisionsUntilMutation() > 0 && 
						population.GetTransfers() < population.GetTotalTransfers() && 
						population.GetKeepTransferring() ) 
				 {
					 population.CalculateDivisions();
					 
					 if( population.GetDivisionsUntilMutation() <= 0) { population.NewMutate(randgen, newtree, node_id); }
					 
					 if( population.GetPopulationSize() >= population.GetPopSizeBeforeDilution()) {
						 population.FrequenciesPerTransferPerNode(newtree, frequencies);
						 population.Resample(randgen); 
						 count++;
						 Cout << Endl << "Passing.... " << count << Endl;
					 }
                 }
            }
      Cout << Endl << Endl;
      population.RunSummary();
      population.PushBackRuns();
      //kptree::print_tree_bracketed(newtree);
      Cout << Endl << Endl << "Printing to screen.... " << Endl;
      population.PrintFrequenciesToScreen(frequencies);
      population.PrintOut(output_file, frequencies);
   }
	 
   //population.PrintOut(output_file, frequencies);
}
