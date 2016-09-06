#include "common.h"
#include "cPopulation.h"

using namespace bpopsim;
using namespace std;


cPopulation::cPopulation(AnyOption& options, gsl_rng* in_rng, uint32_t in_replicate)
: simulation_parameters(options, in_rng)
, output_parameters(options)
, marker_divergence(options)
, replicate(in_replicate)
, unique_genotype_count(0)
, existing_genotype_count(0)
, maximum_subpopulation_fitness(simulation_parameters.initial_fitness)
, current_population_size(0)
, num_completed_transfers(-output_parameters.burn_in)   // This prevents recording
, average_subpopulation_fitness(0)
, run_end_condition_met(false)
, num_divisions_until_next_mutation(0)
, num_completed_divisions(0)
, rng(in_rng)
, log_2(log(2))
{
  debug = options.count("debug");
}

void cPopulation::DisplayParameters()
{
  simulation_parameters.Print();
}

void cPopulation::RunSimulation()
{
  //Initialize the population
  if (marker_divergence.marker_states) {
    SeedPopulationWithMarkedGenotypes(marker_divergence.marker_states); 
  } else {
    SeedPopulationWithOneGenotype();
  }
  
  //Get an initial time point
  RecordStatisticsAtTransfer();
  
  cerr << "    Replicate: " << setw(3) << replicate << "  Transfer: " << setw(4) << num_completed_transfers;
  cerr << "  Fitness: " << setw(7) << fixed << setprecision(4) << average_subpopulation_fitness << "  Genotypes: " << setw(6) << existing_genotype_count; 
  cerr << endl;
  
  // Print the initial tree
  //if (g_verbose) population.PrintTree();
  
  // Get us started with a mutation
  num_divisions_until_next_mutation += CalculateDivisionsUntilNextBeneficialMutation();

  run_end_condition_met = false;
  while(!run_end_condition_met && (num_completed_transfers < simulation_parameters.maximum_number_of_transfers)  ) {
        
    ProcessCellDivisionTimeStepExactWithFractionalCells();
    
    //cerr << current_population_size << endl;
    
    // Debug check on mutate function
    double cells_before_mutate;
    if (debug) cells_before_mutate = CalculatePopulationSize();
    
    // Slightly greater than one for floating point errors.
    if ( num_divisions_until_next_mutation <= 0.00001 ) {
      MutateExactWithFractionalCells();
    }
    num_divisions_until_next_mutation += CalculateDivisionsUntilNextBeneficialMutation();

    
    if (debug) {
      double cells_after_mutate = CalculatePopulationSize();
      assert( abs(cells_before_mutate - cells_after_mutate) < 0.00001);
    }
      
    if( current_population_size >= simulation_parameters.final_population_size_at_transfer ) {

      // @JEB: Would need to re-implement mutation accumulation version
      /*
        if( print_single_fit )
          TransferResampleExactlyOne();
        else {
      */ 
      TransferResampleDilution();
      CullSubpopulationsThatDidNotEstablish();
      //  }
      
      // Record statistics after the transfer
      RecordStatisticsAtTransfer();
        
      /* @JEB: Would need to re-implement...
       
        if( print_single_fit && !use_mute_num) {
          population.PrintSingleFitness(options["output-folder"]);
          //cout << "Population size: " << population.GetPopulationSize() << endl;
        }
       */
      
      cerr << "    Replicate: " << setw(3) << replicate << "  Transfer: " << setw(4) << num_completed_transfers;
      cerr << "  Fitness: " << setw(7) << fixed << setprecision(4) << average_subpopulation_fitness << "  Genotypes: " << setw(6) << existing_genotype_count; 
      cerr << endl;
      
      if ((simulation_parameters.average_fitness_end_condition != 0) && (average_subpopulation_fitness >= simulation_parameters.average_fitness_end_condition))
        run_end_condition_met = true;
    }
  }
  
  RecordStatisticsAtEnd();
}

void cPopulation::DisplaySimulationSummary()
{
  cerr << "Total mutations: " << replicate_statistics.total_mutations << endl;
  cerr << "Total subpopulations lost: " << replicate_statistics.total_subpopulations_lost << endl;
  cerr << "Transfers: " << num_completed_transfers << endl;
  cerr << "Maximum Fitness: " << maximum_subpopulation_fitness   << endl;
}

//@agm I comandeered this function to print stuff out in the manner I see fit

/***** Important for post-hoc use with R ******
 The output file will have a header and it
 will have the largest number of columns in file.
 Thus, in R, you can simply set header to 
 true, and fill NA spaces as the read.table
 function allows. */

void cPopulation::OutputCladeFrequencies()
{  
  //vector<uint32_t> all_sweep_ids = GenotypesFromAncestorToFinalDominant(frequency_threshold);
  
  vector<uint32_t> all_sweep_ids = CladesAboveThreshold(output_parameters.minimum_output_frequency);

  set<uint32_t> lod = GenotypesFromAncestorToFinalDominant();
  
  //Print each sweeping genotype and its frequency per time
	ofstream output_file;
  string output_file_name = output_parameters.output_directory_name + "/clade_frequencies_" 
    + to_string(output_parameters.minimum_output_frequency) + "_" + to_string(replicate) + ".dat";
	output_file.open(output_file_name.c_str(),ios::out);///removed ",4" from the minimum_output_frequency to_string to distinguish frequencies lower than 1E-4
    
  cerr << "Output: " << output_file_name << endl;
  
  output_file << "#on_dominant_line_of_descent";
  for (uint32_t a_node=0; a_node<all_sweep_ids.size(); a_node++) {
    output_file << "\t" <<  (lod.count(all_sweep_ids[a_node]) ? 1 : 0);
  }
  output_file << endl;
  
  output_file << "#fitness_effect";
  for (uint32_t a_node=0; a_node<all_sweep_ids.size(); a_node++) {
    output_file << "\t" <<  genotype_id_map_into_tree[all_sweep_ids[a_node]]->this_mutation_fitness_effect;
  }
  output_file << endl;
  
  output_file << "transfer";
  for (uint32_t a_node=0; a_node<all_sweep_ids.size(); a_node++) {
    output_file << "\t" << "clade_" << all_sweep_ids[a_node];
  }
  output_file << endl;
  
  for (uint32_t transfer=0; transfer<replicate_statistics.clade_frequencies.size(); transfer++) {
    output_file << (transfer * output_parameters.coarse_graining);
    GenotypeFrequencyMap& transfer_clade_frequencies = replicate_statistics.clade_frequencies[transfer];
    for (uint32_t a_node=0; a_node<all_sweep_ids.size(); a_node++) {
      output_file << "\t" << transfer_clade_frequencies.GetFrequency(all_sweep_ids[a_node]);
    }
    output_file << endl;
  }
  
  output_file.close();
}

void cPopulation::OutputGenotypeFrequencies()
{  
  vector<uint32_t> all_sweep_ids = GenotypesAboveThreshold(output_parameters.minimum_output_frequency);
  
  //Print each sweeping genotype and its frequency per time
	ofstream output_file;
  string output_file_name = output_parameters.output_directory_name + "/genotype_frequencies_" 
  + to_string(output_parameters.minimum_output_frequency, 4) + "_" + to_string(replicate) + ".dat";
	output_file.open(output_file_name.c_str(),ios::out);
  
  cerr << "Output: " << output_file_name << endl;
  
  output_file << "transfer";
  for (uint32_t a_node=0; a_node<all_sweep_ids.size(); a_node++) {
    output_file << "\t" << "genotype_" << all_sweep_ids[a_node];
  }
  output_file << endl;
  
  for (uint32_t transfer=0; transfer<replicate_statistics.genotype_frequencies.size(); transfer++) {
    output_file << (transfer * output_parameters.coarse_graining);
    GenotypeFrequencyMap& transfer_genotype_frequencies = replicate_statistics.genotype_frequencies[transfer];
    for (uint32_t a_node=0; a_node<all_sweep_ids.size(); a_node++) {
      output_file << "\t" << transfer_genotype_frequencies.GetFrequency(all_sweep_ids[a_node]);
    }
    output_file << endl;
  }
  
  output_file.close();
}

void cPopulation::OutputMullerMatrix(uint32_t frequency_resolution, bool line_of_descent_separate)
{
  ofstream output_file;
  string output_file_name = output_parameters.output_directory_name + "/muller_matrix_" 
    /* + to_string(frequency_resolution) + "_" */ + to_string(replicate) + ".dat";
	output_file.open(output_file_name.c_str(),ios_base::out);
  
  cerr << "Output: " << output_file_name << endl;
  
  map<uint32_t, int32_t> renumber;
  int32_t renumber_value(0);
  int32_t renumber_lod_value(-1);
  
  // algorithm below is the same no matter what, just depends on whether lod_genotypes is populated
  set<uint32_t> lod_genotypes;
  if (line_of_descent_separate) {
    lod_genotypes = GenotypesFromAncestorToFinalDominant();
  }
  
  uint32_t time = 0;
  
  //step through simulation time
  for (vector<GenotypeFrequencyMap>::iterator this_time_freq = replicate_statistics.clade_frequencies.begin(); 
       this_time_freq != replicate_statistics.clade_frequencies.end(); ++this_time_freq) {
    cerr << " Writing Transfer: " << setw(4) << (time * output_parameters.coarse_graining) << endl;
    time++;
    vector<cFrequencySlice> child_freqs;
    
    tree<cGenotype>::sibling_iterator location;
    location = genotype_tree.begin();
    
    AssignChildFrequency(location, 0, 1, &child_freqs, *this_time_freq);
    sort(child_freqs.begin(), child_freqs.end(), cSortByLow());
        
    // Note: we look at frequencies that are !!centered!! on the locations of the pixels.
    // if you look at frequencies on the bottom or top of pixels then you get a row that
    // has essentailly random genotypes in it on the edge
    double pixel_step = 1.0/frequency_resolution;
    
    // create the pixels
    vector<int32_t> output_pixels(frequency_resolution, 0);
    
    for (uint32_t j=0; j<child_freqs.size(); j++) {
      
      // We look at frequencies that are !!centered!! on the locations of the pixels.
      // this is why we are subtracting half here. 
      int32_t low = ceil(child_freqs[j].low/pixel_step - 0.5);
      int32_t high = ceil(child_freqs[j].high/pixel_step - 0.5);
            
      for(int32_t pixel = low; pixel < high; pixel++) {
        
        if( renumber.count(child_freqs[j].unique_node_id) == 0 ) {
          
          if ( lod_genotypes.count(child_freqs[j].unique_node_id) != 0 ) {
            renumber[child_freqs[j].unique_node_id] = renumber_lod_value--;
          } else {
            renumber[child_freqs[j].unique_node_id] = renumber_value++;
          }
        }
        output_pixels[pixel] = renumber[child_freqs[j].unique_node_id];
      }
    }
    
    for(uint32_t pixel = 0; pixel < output_pixels.size(); pixel++) {
      if (pixel != 0) output_file << " ";
      output_file<< left << setw(6) << output_pixels[pixel];
    }
    
    output_file << endl;
  }
}

void cStatistics::OutputAveragePopulationFitness() {
  
	ofstream output_file;
  string output_file_name = output_directory_name + "/average_fitness.tab";
	output_file.open(output_file_name.c_str(),ios::out);
  
  cerr << "Output: " << output_file_name << endl;
    
  int32_t max_transfers_printed = 0;
  for (uint32_t replicate = 0; replicate < this->size(); ++replicate) {
    cReplicateStatistics& replicate_statistics = (*this)[replicate];
    if (replicate_statistics.average_population_fitness.size() > max_transfers_printed)
      max_transfers_printed = replicate_statistics.average_population_fitness.size();
  }
  
  output_file << "replicate";
  for (uint32_t transfer = 0; transfer < max_transfers_printed; ++transfer) {
    output_file << "\t" << (transfer * coarse_graining);
  }
  output_file << endl;
  
  for (uint32_t replicate = 0; replicate < this->size(); ++replicate) {
    
    cReplicateStatistics& replicate_statistics = (*this)[replicate];
    output_file << (replicate+1);
    
    for (uint32_t transfer = 0; transfer < max_transfers_printed; ++transfer) {
              
      if (transfer < replicate_statistics.average_population_fitness.size() )
        output_file << "\t" << replicate_statistics.average_population_fitness[transfer];
      else
        output_file << "\t" << "NA";

    }
    output_file << endl;
  }
  output_file.close();
}


void cStatistics::OutputAveragePopulationMutationCounts() {
  
  //
  // Main output file of overall totals
  //
  
	ofstream output_file;
  string output_file_name = output_directory_name + "/average_mutations.tab";
	output_file.open(output_file_name.c_str(),ios::out);
  
  cerr << "Output: " << output_file_name << endl;
  
  int32_t num_entries = (*this)[0].average_total_mutation_count.size();
  
  int32_t max_transfers_printed = 0;
  for (uint32_t replicate = 0; replicate < this->size(); ++replicate) {
    cReplicateStatistics& replicate_statistics = (*this)[replicate];
    if (replicate_statistics.average_total_mutation_count.size() > max_transfers_printed)
      max_transfers_printed = replicate_statistics.average_total_mutation_count.size();
  }
  
  output_file << "replicate";
  for (uint32_t transfer = 0; transfer < num_entries; ++transfer) {
    output_file << "\t" << (transfer * coarse_graining);
  }
  output_file << endl;
  
  for (uint32_t replicate = 0; replicate < this->size(); ++replicate) {
    
    cReplicateStatistics& replicate_statistics = (*this)[replicate];
    output_file << (replicate+1);
    
    for (uint32_t transfer = 0; transfer < num_entries; ++transfer) {
      
      if (transfer < replicate_statistics.average_total_mutation_count.size() )
        output_file << "\t" << replicate_statistics.average_total_mutation_count[transfer];
      else
        output_file << "\t" << "NA" << endl;
      
    }
    output_file << endl;
  }
  output_file.close();
  
  //
  // Output files for each category of mutation
  //
  
  for(uint32_t mutation_category=0; mutation_category<(*this)[0].average_mutation_counts[0].size() ; mutation_category++) {
    
    ofstream output_file;
    string output_file_name = output_directory_name + "/average_mutations_category_" + to_string(mutation_category+1) + ".tab";
    output_file.open(output_file_name.c_str(),ios::out);
    
    cerr << "Output: " << output_file_name << endl;
    
    output_file << "transfer";
    for (uint32_t transfer = 0; transfer < num_entries; ++transfer) {
      output_file << "\t" << (transfer * coarse_graining);
    }
    output_file << endl;
    
    for (uint32_t replicate = 0; replicate < this->size(); ++replicate) {
      
      cReplicateStatistics& replicate_statistics = (*this)[replicate];
      output_file << (replicate+1);
      
      for (uint32_t transfer = 0; transfer < num_entries; ++transfer) {
                  
        if (transfer < replicate_statistics.average_mutation_counts.size() )
          output_file << "\t" << replicate_statistics.average_mutation_counts[transfer][mutation_category];
        else
          output_file << "\t" << "NA" << endl;          
      }
      output_file << endl;
    }
    output_file.close();
  }
  
}

void cStatistics::OutputDivergedFrequenciesByDepth()
{
  //
  // Output files for each depth of mutational divergence
  //
  
  int32_t num_entries = (*this)[0].diverged_frequencies_by_depth.size();
  
  int32_t max_transfers_printed = 0;
  for (uint32_t replicate = 0; replicate < this->size(); ++replicate) {
    cReplicateStatistics& replicate_statistics = (*this)[replicate];
    if (replicate_statistics.diverged_frequencies_by_depth.size() > max_transfers_printed)
      max_transfers_printed = replicate_statistics.diverged_frequencies_by_depth.size();
  }
  
  for(uint32_t divergence_depth=0; divergence_depth<(*this)[0].diverged_frequencies_by_depth[0].size() ; divergence_depth++) {
    
    ofstream output_file;
    string output_file_name = output_directory_name + "/mutation_divergence_depth_frequency_" + to_string(divergence_depth+1) + ".tab";
    output_file.open(output_file_name.c_str(),ios::out);
    
    cerr << "Output: " << output_file_name << endl;
    
    output_file << "replicate";
    for (uint32_t transfer = 0; transfer < num_entries; ++transfer) {
      output_file << "\t" << (transfer * coarse_graining);
    }
    output_file << endl;
    
    for (uint32_t replicate = 0; replicate < this->size(); ++replicate) {
      
      cReplicateStatistics& replicate_statistics = (*this)[replicate];
      output_file << (replicate+1);
      
      for (uint32_t transfer = 0; transfer < num_entries; ++transfer) {
        
        if (transfer < replicate_statistics.average_mutation_counts.size() )
          output_file << "\t" << replicate_statistics.diverged_frequencies_by_depth[transfer][divergence_depth];
        else
          output_file << "\t" << "NA" << endl;
        
      }
      output_file << endl;
    }
    output_file.close();
  }
}

// Not implemented
void cStatistics::OutputLargestSweepSize() 
{
  ofstream output_file;
  string output_file_name = output_directory_name + "/largest_sweep_size.tab";
  output_file.open(output_file_name.c_str(),ios::out);
 /* 
  int32_t num_entries = (*this)[0].diverged_frequencies_by_depth.size();

  for (uint32_t replicate = 0; replicate < this->size(); ++replicate) {
    
    cReplicateStatistics& replicate_statistics = (*this)[replicate];
    
    // Loop through the final line of descent
    
    
    output_file << endl;
  }
*/  
  output_file.close();
}

void cStatistics::OutputLongestSweepTime()
{
  ofstream output_file;
  string output_file_name = output_directory_name + "/longest_sweep_time.tab";
  output_file.open(output_file_name.c_str(),ios::out);
  
  for (uint32_t replicate = 0; replicate < this->size(); ++replicate) {
    cReplicateStatistics& replicate_statistics = (*this)[replicate];
    output_file << replicate_statistics.longest_sweep_time << endl;
  }
  
  output_file.close();
}

void cPopulation::UpdateSubpopulationsForGrowthExactWithFractionalCells(double update_time) 
{
  uint32_t i=-1;
  current_population_size = 0; // Update the population size.
  for (vector<cSubpopulation>::iterator it = current_subpopulations.begin(); it!=current_subpopulations.end(); ++it) {
    i++; // must advance iterator before continue statement
    
    // N = No * exp(log(2) * relative_growth_rate * t) 
    double new_number = it->GetNumber() * exp(log_2 * it->GetFitness() * update_time);     
    it->SetNumber(new_number);
    current_population_size += new_number;
  }
}

double cPopulation::TimeToNextWholeCell() 
{
   double time_to_next_whole_cell = -1;
   for (vector<cSubpopulation>::iterator it = current_subpopulations.begin(); it!=current_subpopulations.end(); ++it) { 
   
      if(it->GetNumber() == 0) continue;

      //what is the time to get to the next whole number of cells?     
      double current_cells(it->GetNumber());
      double next_whole_cells( floor(it->GetNumber())+1 );
      
      //cout << "Current cells: " << current_cells << " Next whole cells: " << next_whole_cells << endl;
      // N = No * exp(log(2) * growth_rate * t) 
      double this_time_to_next_whole_cell = log2(next_whole_cells / current_cells) / (it->GetFitness());   

      if ( time_to_next_whole_cell == -1 || (this_time_to_next_whole_cell < time_to_next_whole_cell) ) {
        time_to_next_whole_cell = this_time_to_next_whole_cell;
      }
  }
  
  return time_to_next_whole_cell;
}

 
//
// void cPopulation::ProcessCellDivisionTimeStep()
//
// We would like to move forward by as many divisions as it takes to
// get to either (1) the next mutation OR (2) the next transfer.

void cPopulation::ProcessCellDivisionTimeStepExactWithFractionalCells()
{
  
  if (g_verbose) {
    for(vector<cSubpopulation>::iterator this_time = current_subpopulations.begin(); this_time != current_subpopulations.end(); this_time++) {
      cout << "Genotype: " << this_time->GetNode_id() << " Fitness: " << this_time->GetFitness() << " Number: " << this_time->GetNumber() << endl;
    }
  }
  
  double desired_divisions = num_divisions_until_next_mutation;
  
  // Adjust divisions to be to the end of the transfer if it is larger
  if (desired_divisions + current_population_size > simulation_parameters.final_population_size_at_transfer) {
    desired_divisions = simulation_parameters.final_population_size_at_transfer - current_population_size;
  }
  
  // Note: we underestimate by a few divisions so that we can step forward by single division increments
  // as we get close to the one where the mutation happened (or right before a transfer).      
  
  if (desired_divisions < 1) {
    desired_divisions = 1;
  }
  
  if (g_verbose) {
    cout << "Total pop size: " << current_population_size <<endl;
    cout << "Desired divisions: " << desired_divisions <<endl;
  }
  
  // How much time would we like to pass to achieve the desired number of divisions?
  // (Assuming the entire population has the maximum fitness, makes us underestimate by a few)
  
  double update_time(0);
  
  
  // This is the exact amount of time to be sure that this many fractional cells were produced
  update_time = log2((desired_divisions+current_population_size) / current_population_size) / CalculateAverageSubpopulationFitness();  

  if (g_verbose) {
    cout << "Update time: " << update_time <<endl;    
  }
  
  //Now update all lineages by the time that actually passed  
  double previous_population_size = current_population_size;
  UpdateSubpopulationsForGrowthExactWithFractionalCells(update_time);
  double completed_divisions = current_population_size - previous_population_size;
  //assert(completed_divisions != 0);
  
  // Check approximation -- it's very good.
  //cerr << completed_divisions << " " << desired_divisions << endl;
  
  if (g_verbose) {
    cerr << "Completed divisions: " << completed_divisions << endl;
    for(vector<cSubpopulation>::iterator this_time = current_subpopulations.begin(); this_time != current_subpopulations.end(); this_time++) {
      cerr << "Genotype: " << this_time->GetNode_id() << " Frequency: " << this_time->GetNumber() << " Fitness: " << this_time->GetFitness() << endl;
    }
  }
  
  num_divisions_until_next_mutation -= completed_divisions;
  
  /*
   if (num_divisions_until_next_mutation <= 0)
   if (completed_divisions > 3)
   cerr << completed_divisions << endl;
   */ // We should really backtrack so that all pops divided simultaneously, this is slightly off.
  // We could move time forward by the number of mutations and then weight which subpop it occurs in by the number of divisions in each one
  
  if (g_verbose) {
    cerr << "Divisions until mutation: " << num_divisions_until_next_mutation << endl;
  }
}

void cPopulation::TransferResampleDilution() 
{
  //When it is time for a transfer, resample population
  double by_color[0];
  
  // @JEB: Move to main simulation loop
  num_completed_transfers++;
  
  if (g_verbose) cerr << "  Population size before transfer: " << current_population_size << endl;
  
  current_population_size = 0; // recalculate population size
  for (vector<cSubpopulation>::iterator it = current_subpopulations.begin(); it!=current_subpopulations.end(); ++it) {
    
    // Perform binomial sampling
    it->Transfer(simulation_parameters.binomial_sampling_transfer_probability, rng);
    
    if( it->GetMarker() == 'r' ) by_color[RED] += it->GetNumber();
    else if( it->GetMarker() == 'w' ) by_color[WHITE] += it->GetNumber();
    
    //@agm This section is to delete new mutations from the tree that do not get passed, 
    //     it also deletes subpopulations from the list that have zero population
    
    current_population_size += it->GetNumber();
    
    if( it->GetNumber() == 0 ) {
      if (new_genotype_ids_since_last_transfer.count(it->GetNode_id())) {
        it->GetGenotype().marked_for_deletion = true;
      }
      it = current_subpopulations.erase(it);
      --it;
      replicate_statistics.total_subpopulations_lost++;
      existing_genotype_count--;
    }
  }
  
  if (g_verbose) cerr << "  Population size after transfer: " << current_population_size << endl;
  
  // Check our calculated population size
  // Check how we have updated the population size
  if (debug) {
    double calculated_population_size = CalculatePopulationSize();
    assert(calculated_population_size == current_population_size);
  }
  
  /*
   if (g_verbose) cout << "Colors: " << by_color[RED] << " / " << by_color[WHITE] << endl;
   SetRatio( (double) by_color[RED]/by_color[WHITE] );
   
   // Checks for stopping early if in marker divergence mode
   if (g_ro_only) {
   
   // One color was lost -- bail
   if ( (by_color[RED] == 0) || (by_color[WHITE] == 0) ) {
   run_end_condition_met = true;
   }
   
   // We have the minimum number of transfers to be printed and have diverged sufficiently -- bail
   else {
   if  ( (GetRatio() > GetMaxDivergenceFactor()) || 
   (GetRatio() < 1/GetMaxDivergenceFactor()) )   
   {
   if (g_verbose) cout << "DIVERGENCE CONDITION MET" << endl;
   run_end_condition_met = true;
   }  
   }
   }
   */
  
  //cout << "Pop size: " << m_population_size << endl;
}

//Currently this only works for picking one cell after resample
void cPopulation::TransferResampleExactlyOne() {
  long int starting_size(0), random_number;
  //vector<cSubpopulation> temp_populations;
  
  assert(rng);
  
  num_completed_transfers++; // @JEB should be moved to main loop
  current_population_size = 0; // recalculate population size
  
  for (vector<cSubpopulation>::iterator it = current_subpopulations.begin(); it!=current_subpopulations.end(); ++it) {
    starting_size += it->GetNumber();
  }
  
  //cout << "Population size before dilution: " << starting_size << endl;
  
  //for (uint32_t num_cells = 0; num_cells < m_initial_population_size_per_transfer; num_cells++) {
  random_number = gsl_rng_uniform_int(rng, starting_size);
  
  vector<cSubpopulation>::iterator it;
  for (it = current_subpopulations.begin(); it!=current_subpopulations.end(); ++it) {
    if( random_number < it->GetNumber() ) {
      it->SetNumber(1);
      
      // Erase empty subpopulations
      ++it;
      while(it != current_subpopulations.end()) {
        current_subpopulations.erase(it);
      }
      break;
    }
    else {
      random_number -= it->GetNumber();
      current_subpopulations.erase(it);
      --it;
    }
  }
  
  //}
  
  /*for (vector<cSubpopulation>::iterator it = m_current_subpopulations.begin(); it!=m_current_subpopulations.end(); ++it) {
   cout << it->GetNumber() << " " << it->GetNode_id() << " " << it->GetFitness() << endl;
   }*/
  
  current_population_size = 1;
}


void cPopulation::SeedPopulationWithMarkedGenotypes(uint32_t genetic_marker_states) 
{	
  // For now assumes...
  assert(genetic_marker_states == 2);
  
	cGenotype r(simulation_parameters, unique_genotype_count++);
  cGenotype w(simulation_parameters, unique_genotype_count++);
	tree<cGenotype>::iterator top, red_side, white_side;
	double starting_fitness = simulation_parameters.initial_fitness;
	
	//initialize object of cSubpopulation type
	cSubpopulation red, white;
	
	/*@agm The unique code and fitness is set. */
  
  /*@jeb note that we increment genotype count here.
   It should really happen in AddSubpopulation?
   */
  
	r.fitness = starting_fitness;
  w.fitness = starting_fitness;
  
	//Start building the tree
	red_side = genotype_tree.insert(genotype_tree.begin(), r);
	white_side = genotype_tree.insert(genotype_tree.begin(), w);
	
	red.SetNumber(simulation_parameters.initial_population_size / genetic_marker_states);
	red.SetGenotype(red_side);
	red.SetMarker('r');
	
	white.SetNumber(simulation_parameters.initial_population_size / genetic_marker_states);
	white.SetGenotype(white_side);
	white.SetMarker('w');
	
	AddSubpopulation(red);
	AddSubpopulation(white);	
  
  current_population_size = simulation_parameters.initial_population_size;

}

//@agm This function seeds the population with only one colony to avoid the red/white problem
//     when possible. For this to work properly, I need to get rid of victory conditions.

//update: Victory conditions are now a command line option.  If not set to false the simulation
//        will run to the max number of iterations.

void cPopulation::SeedPopulationWithOneGenotype() {
  tree<cGenotype>::iterator start_position;
  double starting_fitness(simulation_parameters.initial_fitness);
  
  cGenotype neutral(simulation_parameters, unique_genotype_count++);
  neutral.fitness = starting_fitness;
  start_position = genotype_tree.insert(genotype_tree.begin(), neutral);
  
  cSubpopulation begin_here;
  begin_here.SetNumber(simulation_parameters.initial_population_size);
  begin_here.SetGenotype(start_position);
  begin_here.SetMarker('n');
  AddSubpopulation(begin_here);
  
  current_population_size = simulation_parameters.initial_population_size;
  
  if (g_verbose) {
    cout << "Completed divisions: " << num_completed_divisions << endl;
    for(vector<cSubpopulation>::iterator this_time = current_subpopulations.begin(); this_time != current_subpopulations.end(); this_time++) {
      cout << "Genotype1: " << this_time->GetNode_id() << " Frequency1: " << this_time->GetNumber() << endl;
    }
  }
  
}

// Note that if you call add subpopulation, you may need to update the population size!
void cPopulation::AddSubpopulation(cSubpopulation& subpop) 
{
  genotype_id_map_into_tree[subpop.GetGenotype().unique_node_id] = subpop.GetGenotypeIter();
	current_subpopulations.push_back(subpop);
  existing_genotype_count++;
  
  new_genotype_ids_since_last_transfer.insert(subpop.GetNode_id()); // Keep track of all mutations since last transfer
	//cout << subpop.GetNode_id() << " " << subpop.GetFitness() << endl;
}


// Don't cull all populations that disappeared because we need the final tree of all genotypes
void cPopulation::CullSubpopulationsThatDidNotEstablish() {
  
  // Be sure to use leaf_begin() and leaf_end()

  for(tree<cGenotype>::leaf_iterator it = genotype_tree.begin_leaf(); it != genotype_tree.end_leaf();) {
    if (it->marked_for_deletion) {
      tree<cGenotype>::leaf_iterator ite = it;
      it++;
      genotype_tree.erase(ite);
      genotype_id_map_into_tree.erase(ite->unique_node_id);
    } else {
      it++;
    }
  }
  new_genotype_ids_since_last_transfer.clear();
}



//Generates new mutant and adds it to the tree
void cPopulation::MutateExactWithFractionalCells() 
{	
	if (g_verbose) cout << "* Mutating!" << endl;
  if (g_verbose) cout << "Total population: " << CalculatePopulationSize() << endl;
	
  // Mutation depends on instantaneuous Fitness * population size of subpopulations.
  vector<double> fractional_chances_per_subpopulation(current_subpopulations.size());
  
  // All subpopulations can mutate as long as they have more than one cell
  // (so they won't go negative). Create table of their chances.
  double total = 0;
  for(uint32_t i=0; i<current_subpopulations.size(); i++) {
    fractional_chances_per_subpopulation[i] = (current_subpopulations[i].GetNumber() > 1) ? current_subpopulations[i].GetFitness() * current_subpopulations[i].GetNumber() : 0;
    total += fractional_chances_per_subpopulation[i];
  }
  
  for(uint32_t i=0; i<fractional_chances_per_subpopulation.size(); i++) {
    fractional_chances_per_subpopulation[i] /= total;
  }
  
  // Pick which one divided.
  // Gets a random number in range [0,1)
  double random_fraction = gsl_rng_uniform(rng);
  uint32_t chosen_subpopulation = 0;
  while (random_fraction > fractional_chances_per_subpopulation[chosen_subpopulation]) {
    random_fraction -= fractional_chances_per_subpopulation[chosen_subpopulation];
    chosen_subpopulation++;
  }
  
	cSubpopulation& ancestor = current_subpopulations[chosen_subpopulation];          
	cSubpopulation new_subpop;
  
  // There must be at least two cells for a mutation to have occurred...
  if (ancestor.GetNumber() <= 1) {
    cout << "Died: " << ancestor.GetNumber() << " " << chosen_subpopulation << " / "  << current_subpopulations.size() << " " << random_fraction << endl;
  }
  assert(ancestor.GetNumber() > 1.0);
	
  
  bool success_adding = new_subpop.CreateDescendant(
                              rng, 
                              ancestor, 
                              simulation_parameters,
                              genotype_tree,
                              unique_genotype_count++
                              );
  
  if (success_adding) {
  
    if (g_verbose) cout << "  Color: " << new_subpop.GetMarker() << endl;
    if (g_verbose) cout << "  New Fitness: " << new_subpop.GetFitness() << endl;
    
    AddSubpopulation(new_subpop);
    
    // Update the maximum subpopulation fitness
    maximum_subpopulation_fitness = max(maximum_subpopulation_fitness, new_subpop.GetFitness());
    
    if (g_verbose) cout << "  New Number: " << new_subpop.GetNumber() << endl;
    if (g_verbose) cout << "  Old Number: " << ancestor.GetNumber() << endl;
    
    replicate_statistics.total_mutations++;
  }
}

//////////  Methods for recording statistics //////////


void cPopulation::RecordStatisticsAtTransfer()
{	
  
  /////////////////////////////////////////
  /// Statistic: average_population_fitness
  /////////////////////////////////////////
  
  average_subpopulation_fitness = CalculateAverageSubpopulationFitness();
  
  // Only record if we are past a burn_in period
  if (num_completed_transfers <0) return;
  
  // Only record if we are on the coarse-graining interval
  if (num_completed_transfers % output_parameters.coarse_graining != 0)
    return;
  
  replicate_statistics.average_population_fitness.push_back(average_subpopulation_fitness);

  
  //  Debug: this checks to see if our population size was correct.
  if (debug) {
    double calculated_population_size = CalculatePopulationSize();
    assert(current_population_size == calculated_population_size);
  }
  
  /////////////////////////////////////////
  /// Statistic: mutation_counts
  /////////////////////////////////////////
  
  double average_total_mutation_count;  
  vector<double> average_mutation_counts(simulation_parameters.mutation_rates_per_division.size(), 0); 
  
  for (uint32_t it = 0; it<current_subpopulations.size(); ++it) {
    cSubpopulation& this_pop = current_subpopulations[it];
    for (uint32_t on_category = 0; on_category<average_mutation_counts.size(); ++on_category) {
      average_mutation_counts[on_category] += this_pop.GetGenotype().mutation_counts[on_category] * this_pop.GetNumber();
    }
    average_total_mutation_count += this_pop.GetGenotype().total_mutation_count * this_pop.GetNumber();
  }
  
  for (uint32_t on_category = 0; on_category<average_mutation_counts.size(); ++on_category) {
    average_mutation_counts[on_category] /= current_population_size;
  }
  average_total_mutation_count /= current_population_size;
  
  replicate_statistics.average_mutation_counts.push_back(average_mutation_counts);
  replicate_statistics.average_total_mutation_count.push_back(average_total_mutation_count);

  
  /////////////////////////////////////////
  /// Statistic: genotype_frequencies, 
  /// Statistic: clade_frequencies
  /////////////////////////////////////////

  // Initialize all entries
  GenotypeFrequencyMap current_genotype_frequency_map; // only of exactly this genotype 
  GenotypeFrequencyMap current_clade_frequency_map;    // AND including all descendants

	tree<cGenotype>::iterator update_location, genotype;
  for (uint32_t it = 0; it<current_subpopulations.size(); ++it) {
    
    cSubpopulation& this_pop = current_subpopulations[it];
    
    update_location = this_pop.GetGenotypeIter();
    assert(update_location != NULL);
    
    
    double& current_genotype_frequency = current_genotype_frequency_map[update_location->unique_node_id];
    current_genotype_frequency = this_pop.GetNumber() / current_population_size;
    
    // Add frequency to all parents to get clade numbers
		while(update_location != NULL) {
      
      double& current_clade_frequency = current_clade_frequency_map[update_location->unique_node_id];
      current_clade_frequency += this_pop.GetNumber() / current_population_size;
      
      //current_clade_frequencies[update_location->unique_node_id].unique_node_id = update_location->unique_node_id;
      //current_clade_frequencies[update_location->unique_node_id].m_frequency += floor(this_pop.GetNumber());
      
      update_location = genotype_tree.parent(update_location);
		}    
  }
  
  // Push to record
	replicate_statistics.genotype_frequencies.push_back(current_genotype_frequency_map);
	replicate_statistics.clade_frequencies.push_back(current_clade_frequency_map);
}

void cPopulation::RecordStatisticsAtEnd()
{	
  if (output_parameters.output_diverged_frequencies) {
  
    // Record mutations diverged by diverged_mutation_depth at each transfer
    uint32_t max_depth = output_parameters.diverged_mutation_depth;
    
    // Determine the dominant line of descent.
    //set<uint32_t> dominant_clade_set = GenotypesFromAncestorToFinalSweep();
    set<uint32_t> dominant_clade_set = GenotypesFromAncestorToFinalDominant();
    //cerr << "Number in dominant clade: " << dominant_clade_set.size() << endl;
    
    // @JEB - there are two ways to think about this
    //       1) Sum all genotypes that are this many steps diverged from the dominant LOD
    //             use 
    //       2) Count the size of the maximum clade that is diverged by this many mutations from the dominant LOD
    //  #2 is more like the LTEE results we are comparing to.
    
    
    // METHOD 1:
    /*
    // Step through each transfer (well, each recorded time point)
    for (uint32_t transfer=0; transfer < replicate_statistics.genotype_frequencies.size(); transfer++) {
      GenotypeFrequencyMap& this_genotype_frequencies = replicate_statistics.genotype_frequencies[transfer];
      vector<double> add_frequencies(max_depth, 0);    
      
      // Determine the number of steps back to the line of descent and add to this category and higher ones
      for (GenotypeFrequencyMap::iterator it=this_genotype_frequencies.begin(); it!= this_genotype_frequencies.end(); it++) {
        uint32_t this_node_id = it->first;
        
        tree<cGenotype>::iterator itt = FindGenotypeInTreeByID(this_node_id);

        int32_t this_depth = 0;
        while ((itt != NULL) && (!dominant_clade_set.count(itt->unique_node_id))) {
          

          
          if (this_depth < max_depth) {
            add_frequencies[this_depth] += this_genotype_frequencies[this_node_id];
          }
          
          itt = genotype_tree.parent(itt); 
          this_depth++;
        }
        
        assert(itt != NULL);
      } // end genotype loop
      
      replicate_statistics.diverged_frequencies_by_depth.push_back(add_frequencies);
    } // end transfer loop
     */
    
    // METHOD 2:
    // Step through each transfer (well, each recorded time point)
    for (uint32_t transfer=0; transfer < replicate_statistics.clade_frequencies.size(); transfer++) {
      GenotypeFrequencyMap& this_clade_frequencies = replicate_statistics.clade_frequencies[transfer];
      vector<double> add_frequencies(max_depth, 0);    
     
      // Determine the number of steps back to the line of descent and add to this category and higher ones
      for (GenotypeFrequencyMap::iterator it=this_clade_frequencies.begin(); it!= this_clade_frequencies.end(); it++) {
        uint32_t this_node_id = it->first;
     
        tree<cGenotype>::iterator itt = FindGenotypeInTreeByID(this_node_id);
     
        int32_t this_depth = 0;
        while ((itt != NULL) && (!dominant_clade_set.count(itt->unique_node_id))) {
     
          if (this_depth < max_depth) {
            add_frequencies[this_depth] = max(add_frequencies[this_depth],this_clade_frequencies[this_node_id]);
          }
     
        itt = genotype_tree.parent(itt); 
        this_depth++;
      }
     
       assert(itt != NULL);
     } // end genotype loop
     
     replicate_statistics.diverged_frequencies_by_depth.push_back(add_frequencies);
    } // end transfer loop
  }
  
  
  if (output_parameters.output_largest_sweep_size) {
    assert(false); // not implemented
  }
  
  if (output_parameters.output_longest_sweep_time) {
    bool debug = false;

    set<uint32_t> dominant_clade_set = GenotypesFromAncestorToFinalDominant();
    replicate_statistics.longest_sweep_time = 0;
    // for each genotype on the dominant line of descent
    for(set<uint32_t>::iterator it = dominant_clade_set.begin(); it != dominant_clade_set.end(); it++) {
      
      uint32_t on_genotype_id = *it;
      int32_t start_transfer = 0;
      int32_t end_transfer = 0;

      // determine how many transfers it took to go from ≤5% to ≥95%
      if (debug) cout << "== " << on_genotype_id << " ==" << endl;
      for (uint32_t transfer=0; transfer < replicate_statistics.clade_frequencies.size(); transfer++) {
        GenotypeFrequencyMap& this_clade_frequencies = replicate_statistics.clade_frequencies[transfer];
        if (this_clade_frequencies[on_genotype_id] <= 0.05) {
          start_transfer = transfer;
          if (debug) cout << this_clade_frequencies[on_genotype_id] << " (start)" << endl;
        }
        else if (this_clade_frequencies[on_genotype_id] >= 0.95) {
          end_transfer = transfer;
          if (debug) cout << this_clade_frequencies[on_genotype_id] << " (end)" << endl;
          break;
        } else {
          if (debug) cout << this_clade_frequencies[on_genotype_id] << endl;
        }
      }
      
      if (debug) cout << start_transfer << "-" << end_transfer << endl;
      
      int32_t this_sweep_time = end_transfer - start_transfer;
      replicate_statistics.longest_sweep_time  = max(replicate_statistics.longest_sweep_time, 
                                                     this_sweep_time * static_cast<int32_t>(output_parameters.coarse_graining));
    }    
  } // end output_longest_sweep_time
}

//////////  Utility Functions //////////


// Calculates the total population size by iterating over subpops. 
// This is really for debug checking of code only. Normally, use current_population_size.
double cPopulation::CalculatePopulationSize() 
{
  double calculated_population_size = 0;
  for (vector<cSubpopulation>::iterator it = current_subpopulations.begin(); it!=current_subpopulations.end(); ++it) {
    calculated_population_size += it->GetNumber();
  }  
  return calculated_population_size;
}

// Calculates the maximum population fitness by iterating over subpops. 
// This is really for debug checking of code only. Normally, use maximum_subpopulation_fitness.

double cPopulation::CalculateMaximumSubpopulationFitness() {
  double max_fitness = 0.0;
  for(vector<cSubpopulation>::iterator it=current_subpopulations.begin(); it!=current_subpopulations.end(); ++it) {
    max_fitness = max(max_fitness, it->GetFitness());
  }
  return max_fitness;
}

double cPopulation::CalculateAverageSubpopulationFitness() {
  double avg_fitness = 0.0;
  for(vector<cSubpopulation>::iterator it=current_subpopulations.begin(); it!=current_subpopulations.end(); ++it) {
    avg_fitness += it->GetNumber() * it->GetFitness();
  }
  return avg_fitness/current_population_size;
}



// Helper function for producing Muller plot matrix
// Assigns bounds to child genotype centered (and equally dividing) parent genotype slice.
double cPopulation::AssignChildFrequency(tree<cGenotype>::sibling_iterator this_node,
                                              double in_low,
                                              double in_high,
                                              vector<cFrequencySlice> * child_freqs,
                                              GenotypeFrequencyMap &frequencies, 
                                              int depth)  // current depth in tree, defaults to zero
{
  
  //kptree::print_tree_bracketed(*newtree);
  
  //The low for this mutation should be the low of the input interval
  double this_low = in_low;
  double this_high = this_low + frequencies.GetFrequency(this_node->unique_node_id);
  double size_depth1_children(0), half_size_parent_swath((this_high-this_low)/2);
  
  for (tree<cGenotype>::sibling_iterator it_node = genotype_tree.begin(this_node); it_node!=genotype_tree.end(this_node); ++it_node) {
    size_depth1_children += frequencies.GetFrequency(it_node->unique_node_id);
  }
  /*if(size_depth1_children != 0)
   cout << size_depth1_children << endl;*/
  
  double this_bottom_high(this_low + half_size_parent_swath - (size_depth1_children/2));
  double this_top_low(this_bottom_high);
  
  // The swath for this mutation may shrink on the bottom
  // due to its children taking a bite out of it.  
  double last_assigned_child_high = this_top_low;
  for (tree<cGenotype>::sibling_iterator it_node = genotype_tree.begin(this_node); it_node!=genotype_tree.end(this_node); ++it_node) {
    
    // is the frequency > 0? (It may have gone extinct, in which case it is a waste to keep going down the tree!)
    if( frequencies.GetFrequency(it_node->unique_node_id) > 0 ) {
      last_assigned_child_high = AssignChildFrequency(it_node, this_top_low, this_high, child_freqs, frequencies, depth+1);
    }
    
    // Our new low is the last high assigned to a child
    this_top_low = last_assigned_child_high;
  }
  
  // At this point we know the top and bottom of this node...
  // This is where we would want to paint into a bitmap between them!!
  // Draw between this_low and this_high ... rounding to only paint whole numbers ...
  
  (*child_freqs).push_back(cFrequencySlice(this_node->unique_node_id, this_low, this_bottom_high));
  (*child_freqs).push_back(cFrequencySlice(this_node->unique_node_id, this_top_low, this_high));
  
  // print indented version
  /*for (int i=0; i<depth; i++) {
   cout << " ";
   }*/
  //cout << " ID:" << this_node->unique_node_id << " Freq:" << (*frequencies)[this_node->unique_node_id].frequency;
  //cout << " [" << this_low << "," << this_high << "]" << endl;
  
  return this_high;  
}



tree<cGenotype>::iterator cPopulation::FindGenotypeInTreeByID(uint32_t id) 
{
  return genotype_id_map_into_tree[id];
}


// Collect a set of all genotypes from the ancestor to the final clade that achieved 100% of the population

set<uint32_t> cPopulation::GenotypesFromAncestorToFinalSweep() 
{
  
  // Go through current clades and find the one 
  set<uint32_t> genotype_set;
  
  // Find the most derived clade that is 100% of the population 
  tree<cGenotype>::iterator final_sweeping_clade = genotype_tree.end();
  
  uint32_t final_dominant_clade_id = 0;
  
  const double epsilon = 0.000001;
  for (GenotypeFrequencyMap::iterator it = replicate_statistics.clade_frequencies.back().begin(); 
       it != replicate_statistics.clade_frequencies.back().end(); it++) {
    
    // If the clade frequency is 1.0, then we keep the highest genotype ID (guaranteed to be latest)
    if ( abs(it->second - 1.0) < epsilon ) {
      final_dominant_clade_id = max(final_dominant_clade_id, it->first);
    }
  }
  
  // Find it in the tree
  tree<cGenotype>::iterator it = FindGenotypeInTreeByID(final_dominant_clade_id);
  
  // Now add it and all parents to the set
  
  while (it != NULL) {
    genotype_set.insert(it->unique_node_id);
    it = genotype_tree.parent(it);
  }
  
  return genotype_set;
}

// Start at ancestor and take the greatest clade at each step
// to give genotypes on the line of descent to the final dominant

set<uint32_t> cPopulation::GenotypesFromAncestorToFinalDominant() 
{
  // Go through all sub clades (mutants) starting at the ancestor and take the highest-frequency one at each step
  // = dominant line-of-descent
  set<uint32_t> genotype_set;
  
  GenotypeFrequencyMap& frequencies = replicate_statistics.clade_frequencies.back();
  
  tree<cGenotype>::iterator parent_node = genotype_tree.begin();
  genotype_set.insert(parent_node->unique_node_id);

  while (genotype_tree.number_of_children(parent_node) > 0) {
    
    double biggest_child_frequency = -1.0;
    tree<cGenotype>::iterator biggest_child = NULL;
    for (tree<cGenotype>::sibling_iterator it_node = genotype_tree.begin(parent_node); 
         it_node!=genotype_tree.end(parent_node); ++it_node) {
      
      double this_frequency = frequencies.GetFrequency(it_node->unique_node_id);
      if (this_frequency > biggest_child_frequency) {
        biggest_child = it_node;
        biggest_child_frequency = this_frequency;
      }
    }
    
    parent_node = biggest_child;
    genotype_set.insert(parent_node->unique_node_id);
  }
  
  //cerr << genotype_set.size() << endl;
  
  return genotype_set;
  
}

// Returns a list of all clades that were ever above the
// requested threshold at any transfer during the simulation.

vector<uint32_t> cPopulation::CladesAboveThreshold(float threshold) {
  
  set<uint32_t> genotype_set;
  
  for (uint32_t transfer=0; transfer < replicate_statistics.clade_frequencies.size(); transfer++)
  {
    GenotypeFrequencyMap& transfer_clade_frequencies = replicate_statistics.clade_frequencies[transfer];
    
    for (GenotypeFrequencyMap::iterator it = transfer_clade_frequencies.begin(); it != transfer_clade_frequencies.end(); ++it) {
      if (it->second >= threshold)
        genotype_set.insert(it->first);
    }
  }
  
  vector<uint32_t> genotype_list(genotype_set.begin(), genotype_set.end());
  sort(genotype_list.begin(), genotype_list.end());
  return genotype_list;
}

// Returns a list of all genotypes that were ever above the
// requested threshold at any transfer during the simulation.

vector<uint32_t> cPopulation::GenotypesAboveThreshold(float threshold) {
  
  set<uint32_t> genotype_set;
  
  for (uint32_t transfer=0; transfer < replicate_statistics.genotype_frequencies.size(); transfer++)
  {
    GenotypeFrequencyMap& transfer_genotype_frequencies = replicate_statistics.clade_frequencies[transfer];
    
    for (GenotypeFrequencyMap::iterator it = transfer_genotype_frequencies.begin(); it != transfer_genotype_frequencies.end(); ++it) {
      if (it->second >= threshold)
        genotype_set.insert(it->first);
    }
  }
  
  vector<uint32_t> genotype_list(genotype_set.begin(), genotype_set.end());
  sort(genotype_list.begin(), genotype_list.end());
  return genotype_list;
}


uint32_t cPopulation::Last_Sweep(float threshold) {
  
  for( vector<GenotypeFrequencyMap>::reverse_iterator time = replicate_statistics.clade_frequencies.rbegin(); 
      time != replicate_statistics.clade_frequencies.rend(); ++time ) {
    for( GenotypeFrequencyMap::reverse_iterator node = time->rbegin(); node != time->rend(); ++node) {
      if( node->second >= threshold )
        return node->first;
    }
  }
  
  cout << "Nothing met threshold." << std::endl;
  
  return 0;
}


//Used for MA?

void cPopulation::PrintSingleFitness(string output_folder) {
  string output_file;
  ofstream output_handle;
  
  float single_fitness;
  
  output_file = output_folder + "/SingleFitness.dat";
  output_handle.open(output_file.c_str(),ios_base::app);
  
  if(current_subpopulations.size() == 1) {
    for (vector<cSubpopulation>::iterator it = current_subpopulations.begin(); it!=current_subpopulations.end(); ++it) { 
      single_fitness = it->GetFitness();
    }
  }
  else {
    cerr << "There is more than one population remaining... there is a bug." << endl;
    abort();
  }
  
  output_handle << single_fitness << "\n";
  //if( single_fitness != 1 )
    //cout << "Fitness: " << single_fitness << "\n";
  output_handle.close();
}

void cPopulation::ConvertExternalData(const string &input_file) {
  
  /*
  fstream input_handle;
  input_handle.open(input_file.c_str(), fstream::in);
  assert(input_handle.is_open());
  
  char char_line[1000];
  
  uint32_t num_nodes = 1;
  
  while(!input_handle.getline(char_line,1000).eof()) {
    string line(char_line);
    vector<string> line_pieces(split(line, " ")), preformat_tree(split(line_pieces[0], ","));
    
    tree<cGenotype>::iterator top, next;
    
    if( num_nodes == 1 ) { 
      cGenotype origin;
      origin.unique_node_id = num_nodes;
      origin.name = preformat_tree[0];
      origin.fitness = 1;
      origin.mut_num = 0;
      
      top = genotype_tree.begin();
      genotype_tree.insert(top, origin);
      
      replicate_statistics.clade_frequencies.resize(line_pieces.size()-1);
      
      for(uint32_t this_time = 0; this_time < replicate_statistics.clade_frequencies.size(); ++this_time) {
        if( from_string<double>(line_pieces[this_time+1]) != 0 ) {
          cGenotypeFrequency temp_time;
          temp_time.unique_node_id = num_nodes;
          temp_time.name = *(preformat_tree.end()-1);
          temp_time.m_frequency = from_string<double>(line_pieces[this_time+1]);
          replicate_statistics.clade_frequencies[this_time].push_back(temp_time);
        }
      }
    }
    else {
      for(tree<cGenotype>::iterator this_node = genotype_tree.begin(); this_node != genotype_tree.end(); ++this_node) {
        if( *(preformat_tree.end()-2) == (*this_node).name ) {
          cSubpopulation new_child;
          cGenotype child;
          child.name = *(preformat_tree.end()-1);
          child.unique_node_id = num_nodes;
          child.fitness = 1;
          child.mut_num = num_nodes;
          new_child.AddToTree(genotype_tree, this_node, child);
          
          for(uint32_t this_time = 0; this_time < replicate_statistics.clade_frequencies.size(); ++this_time) {
            if( from_string<double>(line_pieces[this_time+1]) != 0 ) {
              cGenotypeFrequency temp_time;
              temp_time.unique_node_id = num_nodes;
              temp_time.name = *(preformat_tree.end()-1);
              temp_time.m_frequency = from_string<double>(line_pieces[this_time+1]);
              replicate_statistics.clade_frequencies[this_time].push_back(temp_time);
            }
          }
        }
      }
    }
    
    num_nodes++;
  }
  
  input_handle.close();
  */
  /*for(tree<cGenotype>::iterator this_node = m_tree.begin(); this_node != m_tree.end(); ++this_node) {
    cout << this_node->unique_node_id << " " << this_node->name << endl;
  }*/
  
  //PrintFreqsQuick();
  //PrintTree();
}




//Utilities Section

//The goal here is the find the fitness of the winning line of descent at every time point
//This allows comparison eventualy winners and losers style

void cPopulation::PrintWinningFitness(string output_folder, uint32_t run_number) {
  
/*
  vector<uint32_t> all_sweep_ids(GenotypesFromAncestorToFinalDominant(1.0));
  
  ofstream output_handle;
  string output_file;
  
  output_file.append(output_folder);
  output_file.append("/Winning_Fitness_");
  output_file.append(to_string(run_number));
  output_file.append(".dat");
	output_handle.open(output_file.c_str(),ios_base::app);
  
  uint32_t count(0);
  
  //First I iterate over the frequency vector per time
  for (vector< vector<cGenotypeFrequency> >::iterator this_time = m_frequencies.begin(); this_time < m_frequencies.end(); ++this_time) {
    if( count % m_coarse_graining == 0 ) {
      
      for (uint32_t a_node=0; a_node<all_sweep_ids.size(); a_node++) {
        
        //Here we grab the frequency of each of the nodes that are known to have swept by the END of the simulation
        double frequency( GenotypeFrequency( *this_time, all_sweep_ids[a_node]) );
        
        //On the first node that has zero frequency,
        //Find the fitness in the tree of the one that came immediately before it
        if( frequency == 0 ) {
          for (tree<cGenotype>::iterator it = m_tree.begin(); it!=m_tree.end(); it++) {
            if( it->unique_node_id == all_sweep_ids[a_node-1] )
               output_handle << it->fitness << "\t";
               
          }
          //After finding the youngest one at this time point in the winning line...
          //Go to the next time point
          break;
        }
      }
    }
    count++;
  }
  
  output_handle.close();
 */
}

void cPopulation::PrintExpectationValue(const string& output_folder) {
  /*
  //Print everything out
	ofstream output_handle;
  string output_file;
  
  output_file.append(output_folder);
  output_file.append("/ExpectationVal.dat");
	output_handle.open(output_file.c_str(),ios_base::app);
  
  double calculated_mutant_population_freq(0);

  for ( vector<cGenotypeFrequency>::iterator it = m_frequencies[m_frequencies.size()-1].begin(); it!=m_frequencies[m_frequencies.size()-1].end(); ++it) {
    if( it->unique_node_id != 0 ) {
      calculated_mutant_population_freq += it->m_frequency;
    }
  }
    
  cout << calculated_mutant_population_freq << endl;
  */
}

void cPopulation::PrintOut_RedWhiteOnly(const string& output_folder, 
                                        vector< vector<double> > * red_white_ratios,
                                        uint32_t transfer_interval_to_print) 
{  
	//Will print out only red and white
	ofstream output_handle;
  string output_file;
  
  output_file = output_folder + "_Gfreqs_RO.dat";
  
	output_handle.open(output_file.c_str(),ios_base::app);
  
  uint32_t largest_replicate(0);
  //find largest replicate
  for (uint32_t replicate = 0; replicate<(*red_white_ratios).size(); replicate++) {
    if( (*red_white_ratios)[replicate].size() > largest_replicate ) largest_replicate = replicate;
  }
  
  output_handle << "transfer";
  for (uint32_t replicate = 0; replicate<(*red_white_ratios).size(); replicate++) {
    output_handle << "\t" << replicate ;
  }

  output_handle << endl;
  
  for (uint32_t time = 0; time<(*red_white_ratios)[largest_replicate].size(); time++) {
    output_handle << time*transfer_interval_to_print;
    for (uint32_t replicate = 0; replicate<(*red_white_ratios).size(); replicate++) {
      if( time >= (*red_white_ratios)[replicate].size() )
        output_handle << "\t";
      else
        output_handle << "\t" << log((*red_white_ratios)[replicate][time]);
    }
    output_handle << endl;
	}
}

void cPopulation::PrintFrequenciesToScreen_RedWhiteOnly(string output_folder) {;}

//@agm This function determines the maximum difference in genotype frequency between a mutation
//     and its predecessor if they got above the threshold passed to the threshold function below

double cPopulation::CalculateSimilarity(string output_folder) {
  
  /*
  vector<uint32_t> all_sweep_ids(GenotypesFromAncestorToFinalDominant(.99));
  
  ofstream output_handle;
  string output_file;
  
  vector<float> max_diff(all_sweep_ids.size(), 0);
  float current_diff;
  double num_below_threshold(0);
  
  uint32_t i=0;
  
  for (uint32_t a_node=0; a_node<all_sweep_ids.size(); a_node++) {
    uint32_t time = 0;
    for ( vector<cGenotypeFrequency>::iterator it = run_statistics.clade_frequencies[i].begin(); it!=run_statistics.clade_frequencies[i].end(); ++it) {
      
      //This conditional is to coarse-grain the sampling to every 75 transfers (~500 generations)
      //and to make sure there is a value in the vector at that position
      //another way to think about it is this conditional is checking to make sure that
      //a particular mutation has arisen in time before comparing it to something.
      
      if ( time% output_parameters.coarse_graining == 0 ) {
        current_diff = fabs(GenotypeFrequency(*this_time, all_sweep_ids[a_node]) - GenotypeFrequency(*this_time, all_sweep_ids[a_node+1]));
        if( max_diff[i] < current_diff ) max_diff[i] = current_diff;
      }
      time++;
    }
    i++;
  }
  
  output_file = output_folder + "/" + "SweepClumpiness.dat";
  output_handle.open(output_file.c_str(), ios_base::app);
  
  uint32_t num_simultaneous(1);

  for (int i = 0; i<max_diff.size()-1; i++) {
    if( max_diff[i] <= .1 ) num_simultaneous++; 
    else {
      output_handle << num_simultaneous << endl;
      num_simultaneous = 1;
    }
  }
  
  output_handle.close();
  
  output_file = output_folder + "/" + "SignificantParallelMutations.dat";
	output_handle.open(output_file.c_str(), ios_base::app);
  
  for (int i = 0; i<max_diff.size()-1; i++) {
    cout << endl << i << " " << max_diff[i] << endl;
    output_handle << max_diff[i] << endl;
    if( max_diff[i] <= .1 ) num_below_threshold++;
  }
  
  output_handle.close();
  return num_below_threshold;
   */
}

void cPopulation::TimeToSweep(string output_folder) {
  
  /*
  vector<uint32_t> all_sweep_ids(GenotypesFromAncestorToFinalDominant(.98));
  
  ofstream output_handle;
  string output_file;
  
  vector<uint32_t> time_to_sweep(all_sweep_ids.size(), 0);
  
  uint32_t i=0;
  
  for (uint32_t a_node=0; a_node<all_sweep_ids.size(); a_node++) {
    uint32_t time = 0;
    for (vector< vector<cGenotypeFrequency> >::iterator this_time = m_frequencies.begin(); this_time < m_frequencies.end(); ++this_time) {
      
      //This conditional is to coarse-grain the sampling to every 75 transfers (~500 generations)
      //and to make sure there is a value in the vector at that position
      //another way to think about it is this conditional is checking to make sure that
      //a particular mutation has arisen in time before comparing it to something.
      
      double frequency( GenotypeFrequency(*this_time, all_sweep_ids[a_node]) );
      
      if ( time%m_coarse_graining == 0 && frequency != 1 && frequency > .01) {
        time_to_sweep[a_node]+=3;
      }
      time++;
    }
    i++;
  }
  
  //@agm this should output exactly the same thing to file as it does to screen
  //     the virtue of this format is that it is human readable
  output_file.append(output_folder);
  output_file.append("/");
  output_file.append("TimeToSweep.dat");
	output_handle.open(output_file.c_str(),ios_base::app);
  
  for (vector<uint32_t>::iterator it_node = time_to_sweep.begin(); it_node != time_to_sweep.end(); it_node++) {
    if( (*it_node) != 0 ) {
      cout << (*it_node) << endl;
      output_handle << (*it_node) << endl;
    }
  }
}

vector<uint32_t> cPopulation::GenotypesFromAncestorToFinalDominant(float threshold) {
  
  uint32_t youngest_sweep(Last_Sweep(threshold));
  
  vector<uint32_t> all_sweep_ids;
  
  for (tree<cGenotype>::iterator it = m_tree.begin(); it!=m_tree.end(); it++) {
    if( it->unique_node_id == youngest_sweep ) {
      while (it != NULL) {
        all_sweep_ids.push_back(it->unique_node_id);
        it = m_tree.parent(it);
      }
      break;
    }
  }
  
  sort(all_sweep_ids.begin(), all_sweep_ids.end());
  
  return all_sweep_ids;
}

uint32_t cPopulation::Last_Sweep(float threshold) {
  
  for( vector< vector<cGenotypeFrequency> >::reverse_iterator time = m_frequencies.rbegin(); time < m_frequencies.rend(); ++time ) {
    for( vector<cGenotypeFrequency>::reverse_iterator node = time->rbegin(); node < time->rend(); ++node) {
      if( node->m_frequency >= threshold )
        return node->unique_node_id;
    }
  }
  
  cout << "Nothing met threshold." << endl;
  
  return 0;
   */
}




////////// THINGS WE DONT REALLY NEED BELOW


void cPopulation::PrintFreqsQuick() {
  uint32_t time_keeper(0);
  
  for(vector<GenotypeFrequencyMap>::iterator this_time=replicate_statistics.clade_frequencies.begin(); 
      this_time!=replicate_statistics.clade_frequencies.end(); ++this_time) {
    time_keeper++;
    cout << "Time: " << time_keeper << endl;
    for(GenotypeFrequencyMap::iterator this_genotype=this_time->begin(); this_genotype!=this_time->end(); ++this_genotype) {
      cout << this_genotype->first << " " << this_genotype->second << endl;
    }
    cout << endl << endl;
  }
}
