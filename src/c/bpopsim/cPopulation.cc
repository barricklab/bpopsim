#include "common.h"
#include "cPopulation.h"

#include <algorithm>
#include <utility>

using namespace bpopsim;
using namespace std;

cPopulation::cPopulation(SimulationParameters in_simulation_parameters, const OutputParameters& in_output_parameters,
                         const MarkerDivergence& in_marker_divergence, RandomGenerator& in_rng, uint32_t in_replicate)
: simulation_parameters(std::move(in_simulation_parameters))
, output_parameters(in_output_parameters)
, marker_divergence(in_marker_divergence)
, debug(in_output_parameters.debug)
, replicate(in_replicate)
, maximum_subpopulation_fitness(simulation_parameters.initial_fitness)
, num_completed_transfers(-output_parameters.burn_in)  // This prevents recording
, rng(in_rng)
, log_2(log(2))
{}

void cPopulation::RunSimulation()
{
  // Initialize the population
  if (marker_divergence.marker_states) {
    SeedPopulationWithMarkedGenotypes(marker_divergence.marker_states);
  } else {
    SeedPopulationWithOneGenotype();
  }

  // Get an initial time point
  RecordStatisticsAtTransfer();

  cerr << "    Replicate: " << setw(3) << replicate << "  Transfer: " << setw(4) << num_completed_transfers;
  cerr << "  Fitness: " << setw(7) << fixed << setprecision(4) << average_subpopulation_fitness
       << "  Genotypes: " << setw(6) << existing_genotype_count;
  cerr << '\n';

  // Print the initial tree
  // if (output_parameters.verbose) population.PrintTree();

  // Get us started with a mutation
  num_divisions_until_next_mutation += CalculateDivisionsUntilNextBeneficialMutation();

  run_end_condition_met = false;
  while (!run_end_condition_met && (num_completed_transfers < simulation_parameters.maximum_number_of_transfers)) {

    ProcessCellDivisionTimeStepExactWithFractionalCells();

    // cerr << current_population_size << endl;

    // Debug check on mutate function
    double cells_before_mutate = 0;
    if (debug) cells_before_mutate = CalculatePopulationSize();

    // Slightly greater than one for floating point errors.
    if (num_divisions_until_next_mutation <= 0.00001) {
      MutateExactWithFractionalCells();
    }
    num_divisions_until_next_mutation += CalculateDivisionsUntilNextBeneficialMutation();

    if (debug) {
      double cells_after_mutate = CalculatePopulationSize();
      assert(abs(cells_before_mutate - cells_after_mutate) < 0.00001);
    }

    if (current_population_size >= simulation_parameters.final_population_size_at_transfer) {

      TransferResampleDilution();
      CullSubpopulationsThatDidNotEstablish();

      // Record statistics after the transfer
      RecordStatisticsAtTransfer();

      cerr << "    Replicate: " << setw(3) << replicate << "  Transfer: " << setw(4) << num_completed_transfers;
      cerr << "  Fitness: " << setw(7) << fixed << setprecision(4) << average_subpopulation_fitness
           << "  Genotypes: " << setw(6) << existing_genotype_count;
      cerr << '\n';

      if ((simulation_parameters.average_fitness_end_condition != 0) &&
          (average_subpopulation_fitness >= simulation_parameters.average_fitness_end_condition))
        run_end_condition_met = true;
    }
  }

  RecordStatisticsAtEnd();
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
  // vector<uint32_t> all_sweep_ids = GenotypesFromAncestorToFinalDominant(frequency_threshold);

  vector<uint32_t> all_sweep_ids = CladesAboveThreshold(output_parameters.minimum_output_frequency);

  set<uint32_t> lod = GenotypesFromAncestorToFinalDominant();

  // Print each sweeping genotype and its frequency per time
  ofstream output_file;
  string output_file_name = output_parameters.output_directory_name + "/clade_frequencies_" +
                            to_string(output_parameters.minimum_output_frequency) + "_" + to_string(replicate) + ".dat";
  output_file.open(output_file_name.c_str(), ios::out);  /// removed ",4" from the minimum_output_frequency to_string to
                                                         /// distinguish frequencies lower than 1E-4

  cerr << "Output: " << output_file_name << '\n';

  output_file << "#on_dominant_line_of_descent";
  for (unsigned int all_sweep_id : all_sweep_ids) {
    output_file << "\t" << (lod.count(all_sweep_id) ? 1 : 0);
  }
  output_file << '\n';

  output_file << "#fitness_effect";
  for (unsigned int all_sweep_id : all_sweep_ids) {
    output_file << "\t" << genotype_id_map_into_tree[all_sweep_id]->this_mutation_fitness_effect;
  }
  output_file << '\n';

  output_file << "transfer";
  for (unsigned int all_sweep_id : all_sweep_ids) {
    output_file << "\t" << "clade_" << all_sweep_id;
  }
  output_file << '\n';

  for (uint32_t transfer = 0; transfer < replicate_statistics.clade_frequencies.size(); transfer++) {
    output_file << (transfer * output_parameters.coarse_graining);
    GenotypeFrequencyMap& transfer_clade_frequencies = replicate_statistics.clade_frequencies[transfer];
    for (unsigned int all_sweep_id : all_sweep_ids) {
      output_file << "\t" << transfer_clade_frequencies.GetFrequency(all_sweep_id);
    }
    output_file << '\n';
  }

  output_file.close();
}

void cPopulation::OutputGenotypeFrequencies()
{
  vector<uint32_t> all_sweep_ids = GenotypesAboveThreshold(output_parameters.minimum_output_frequency);

  // Print each sweeping genotype and its frequency per time
  ofstream output_file;
  string output_file_name = output_parameters.output_directory_name + "/genotype_frequencies_" +
                            to_string(output_parameters.minimum_output_frequency, 4) + "_" + to_string(replicate) +
                            ".dat";
  output_file.open(output_file_name.c_str(), ios::out);

  cerr << "Output: " << output_file_name << '\n';

  output_file << "transfer";
  for (unsigned int all_sweep_id : all_sweep_ids) {
    output_file << "\t" << "genotype_" << all_sweep_id;
  }
  output_file << '\n';

  for (uint32_t transfer = 0; transfer < replicate_statistics.genotype_frequencies.size(); transfer++) {
    output_file << (transfer * output_parameters.coarse_graining);
    GenotypeFrequencyMap& transfer_genotype_frequencies = replicate_statistics.genotype_frequencies[transfer];
    for (unsigned int all_sweep_id : all_sweep_ids) {
      output_file << "\t" << transfer_genotype_frequencies.GetFrequency(all_sweep_id);
    }
    output_file << '\n';
  }

  output_file.close();
}

void cPopulation::OutputMullerMatrix(uint32_t frequency_resolution, bool line_of_descent_separate)
{
  ofstream output_file;
  string output_file_name = output_parameters.output_directory_name +
                            "/muller_matrix_"
                            /* + to_string(frequency_resolution) + "_" */
                            + to_string(replicate) + ".dat";
  output_file.open(output_file_name.c_str(), ios_base::out);

  cerr << "Output: " << output_file_name << '\n';

  map<uint32_t, int32_t> renumber;
  int32_t renumber_value(0);
  int32_t renumber_lod_value(-1);

  // algorithm below is the same no matter what, just depends on whether lod_genotypes is populated
  set<uint32_t> lod_genotypes;
  if (line_of_descent_separate) {
    lod_genotypes = GenotypesFromAncestorToFinalDominant();
  }

  uint32_t time = 0;

  // step through simulation time
  for (auto& clade_frequencie : replicate_statistics.clade_frequencies) {
    cerr << " Writing Transfer: " << setw(4) << (time * output_parameters.coarse_graining) << '\n';
    time++;
    vector<cFrequencySlice> child_freqs;

    tree<cGenotype>::sibling_iterator location;
    location = genotype_tree.begin();

    AssignChildFrequency(location, 0, 1, &child_freqs, clade_frequencie);
    sort(child_freqs.begin(), child_freqs.end(), cSortByLow());

    // Note: we look at frequencies that are !!centered!! on the locations of the pixels.
    // if you look at frequencies on the bottom or top of pixels then you get a row that
    // has essentailly random genotypes in it on the edge
    double pixel_step = 1.0 / frequency_resolution;

    // create the pixels
    vector<int32_t> output_pixels(frequency_resolution, 0);

    for (auto& child_freq : child_freqs) {

      // We look at frequencies that are !!centered!! on the locations of the pixels.
      // this is why we are subtracting half here.
      int32_t low = ceil((child_freq.low / pixel_step) - 0.5);
      int32_t high = ceil((child_freq.high / pixel_step) - 0.5);

      for (int32_t pixel = low; pixel < high; pixel++) {

        if (renumber.count(child_freq.unique_node_id) == 0) {

          if (lod_genotypes.count(child_freq.unique_node_id) != 0) {
            renumber[child_freq.unique_node_id] = renumber_lod_value--;
          } else {
            renumber[child_freq.unique_node_id] = renumber_value++;
          }
        }
        output_pixels[pixel] = renumber[child_freq.unique_node_id];
      }
    }

    for (uint32_t pixel = 0; pixel < output_pixels.size(); pixel++) {
      if (pixel != 0) output_file << " ";
      output_file << left << setw(6) << output_pixels[pixel];
    }

    output_file << '\n';
  }
}

void cStatistics::OutputAveragePopulationFitness()
{

  ofstream output_file;
  string output_file_name = output_directory_name + "/average_fitness.tab";
  output_file.open(output_file_name.c_str(), ios::out);

  cerr << "Output: " << output_file_name << '\n';

  uint32_t max_transfers_printed = 0;
  for (auto& replicate_statistics : *this) {
    max_transfers_printed = std::max<std::vector<double>::size_type>(
        replicate_statistics.average_population_fitness.size(), max_transfers_printed);
  }

  output_file << "replicate";
  for (uint32_t transfer = 0; transfer < max_transfers_printed; ++transfer) {
    output_file << "\t" << (transfer * coarse_graining);
  }
  output_file << '\n';

  for (uint32_t replicate = 0; replicate < this->size(); ++replicate) {

    cReplicateStatistics& replicate_statistics = (*this)[replicate];
    output_file << (replicate + 1);

    for (uint32_t transfer = 0; transfer < max_transfers_printed; ++transfer) {

      if (transfer < replicate_statistics.average_population_fitness.size())
        output_file << "\t" << replicate_statistics.average_population_fitness[transfer];
      else
        output_file << "\t" << "NA";
    }
    output_file << '\n';
  }
  output_file.close();
}

void cStatistics::OutputAveragePopulationMutationCounts()
{

  //
  // Main output file of overall totals
  //

  ofstream output_file;
  string output_file_name = output_directory_name + "/average_mutations.tab";
  output_file.open(output_file_name.c_str(), ios::out);

  cerr << "Output: " << output_file_name << '\n';

  uint32_t num_entries = (*this)[0].average_total_mutation_count.size();

  uint32_t max_transfers_printed = 0;
  for (auto& replicate_statistics : *this) {
    max_transfers_printed = std::max<std::vector<double>::size_type>(
        replicate_statistics.average_total_mutation_count.size(), max_transfers_printed);
  }

  output_file << "replicate";
  for (uint32_t transfer = 0; transfer < num_entries; ++transfer) {
    output_file << "\t" << (transfer * coarse_graining);
  }
  output_file << '\n';

  for (uint32_t replicate = 0; replicate < this->size(); ++replicate) {

    cReplicateStatistics& replicate_statistics = (*this)[replicate];
    output_file << (replicate + 1);

    for (uint32_t transfer = 0; transfer < num_entries; ++transfer) {

      if (transfer < replicate_statistics.average_total_mutation_count.size())
        output_file << "\t" << replicate_statistics.average_total_mutation_count[transfer];
      else
        output_file << "\t" << "NA" << '\n';
    }
    output_file << '\n';
  }
  output_file.close();

  //
  // Output files for each category of mutation
  //

  for (uint32_t mutation_category = 0; mutation_category < (*this)[0].average_mutation_counts[0].size();
       mutation_category++) {

    ofstream output_file;
    string output_file_name =
        output_directory_name + "/average_mutations_category_" + to_string(mutation_category + 1) + ".tab";
    output_file.open(output_file_name.c_str(), ios::out);

    cerr << "Output: " << output_file_name << '\n';

    output_file << "transfer";
    for (uint32_t transfer = 0; transfer < num_entries; ++transfer) {
      output_file << "\t" << (transfer * coarse_graining);
    }
    output_file << '\n';

    for (uint32_t replicate = 0; replicate < this->size(); ++replicate) {

      cReplicateStatistics& replicate_statistics = (*this)[replicate];
      output_file << (replicate + 1);

      for (uint32_t transfer = 0; transfer < num_entries; ++transfer) {

        if (transfer < replicate_statistics.average_mutation_counts.size())
          output_file << "\t" << replicate_statistics.average_mutation_counts[transfer][mutation_category];
        else
          output_file << "\t" << "NA" << '\n';
      }
      output_file << '\n';
    }
    output_file.close();
  }
}

void cStatistics::OutputDivergedFrequenciesByDepth()
{
  //
  // Output files for each depth of mutational divergence
  //

  uint32_t num_entries = (*this)[0].diverged_frequencies_by_depth.size();

  uint32_t max_transfers_printed = 0;
  for (auto& replicate_statistics : *this) {
    max_transfers_printed = std::max<std::vector<class std::vector<double>>::size_type>(
        replicate_statistics.diverged_frequencies_by_depth.size(), max_transfers_printed);
  }

  for (uint32_t divergence_depth = 0; divergence_depth < (*this)[0].diverged_frequencies_by_depth[0].size();
       divergence_depth++) {

    ofstream output_file;
    string output_file_name =
        output_directory_name + "/mutation_divergence_depth_frequency_" + to_string(divergence_depth + 1) + ".tab";
    output_file.open(output_file_name.c_str(), ios::out);

    cerr << "Output: " << output_file_name << '\n';

    output_file << "replicate";
    for (uint32_t transfer = 0; transfer < num_entries; ++transfer) {
      output_file << "\t" << (transfer * coarse_graining);
    }
    output_file << '\n';

    for (uint32_t replicate = 0; replicate < this->size(); ++replicate) {

      cReplicateStatistics& replicate_statistics = (*this)[replicate];
      output_file << (replicate + 1);

      for (uint32_t transfer = 0; transfer < num_entries; ++transfer) {

        if (transfer < replicate_statistics.average_mutation_counts.size())
          output_file << "\t" << replicate_statistics.diverged_frequencies_by_depth[transfer][divergence_depth];
        else
          output_file << "\t" << "NA" << '\n';
      }
      output_file << '\n';
    }
    output_file.close();
  }
}

// Not implemented
void cStatistics::OutputLargestSweepSize() const
{
  ofstream output_file;
  string output_file_name = output_directory_name + "/largest_sweep_size.tab";
  output_file.open(output_file_name.c_str(), ios::out);
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
  output_file.open(output_file_name.c_str(), ios::out);

  for (auto& replicate_statistics : *this) {
    output_file << replicate_statistics.longest_sweep_time << '\n';
  }

  output_file.close();
}

void cPopulation::UpdateSubpopulationsForGrowthExactWithFractionalCells(double update_time)
{
  current_population_size = 0;  // Update the population size.
  for (auto& current_subpopulation : current_subpopulations) {
    // N = No * exp(log(2) * relative_growth_rate * t)
    double new_number =
        current_subpopulation.GetNumber() * exp(log_2 * current_subpopulation.GetFitness() * update_time);
    current_subpopulation.SetNumber(new_number);
    current_population_size += new_number;
  }
}

double cPopulation::TimeToNextWholeCell()
{
  double time_to_next_whole_cell = -1;
  for (auto& current_subpopulation : current_subpopulations) {

    if (current_subpopulation.GetNumber() == 0) continue;

    // what is the time to get to the next whole number of cells?
    double current_cells(current_subpopulation.GetNumber());
    double next_whole_cells(floor(current_subpopulation.GetNumber()) + 1);

    // cout << "Current cells: " << current_cells << " Next whole cells: " << next_whole_cells << endl;
    //  N = No * exp(log(2) * growth_rate * t)
    double this_time_to_next_whole_cell = log2(next_whole_cells / current_cells) / (current_subpopulation.GetFitness());

    if (time_to_next_whole_cell == -1 || (this_time_to_next_whole_cell < time_to_next_whole_cell)) {
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

  if (output_parameters.verbose) {
    for (auto& current_subpopulation : current_subpopulations) {
      cout << "Genotype: " << current_subpopulation.GetNode_id() << " Fitness: " << current_subpopulation.GetFitness()
           << " Number: " << current_subpopulation.GetNumber() << '\n';
    }
  }

  double desired_divisions = num_divisions_until_next_mutation;

  // Adjust divisions to be to the end of the transfer if it is larger
  if (desired_divisions + current_population_size > simulation_parameters.final_population_size_at_transfer) {
    desired_divisions = simulation_parameters.final_population_size_at_transfer - current_population_size;
  }

  // Note: we underestimate by a few divisions so that we can step forward by single division increments
  // as we get close to the one where the mutation happened (or right before a transfer).

  desired_divisions = std::max<double>(desired_divisions, 1);

  if (output_parameters.verbose) {
    cout << "Total pop size: " << current_population_size << '\n';
    cout << "Desired divisions: " << desired_divisions << '\n';
  }

  // How much time would we like to pass to achieve the desired number of divisions?
  // (Assuming the entire population has the maximum fitness, makes us underestimate by a few)

  double update_time(0);

  // This is the exact amount of time to be sure that this many fractional cells were produced
  update_time = log2((desired_divisions + current_population_size) / current_population_size) /
                CalculateAverageSubpopulationFitness();

  if (output_parameters.verbose) {
    cout << "Update time: " << update_time << '\n';
  }

  // Now update all lineages by the time that actually passed
  double previous_population_size = current_population_size;
  UpdateSubpopulationsForGrowthExactWithFractionalCells(update_time);
  double completed_divisions = current_population_size - previous_population_size;
  // assert(completed_divisions != 0);

  // Check approximation -- it's very good.
  // cerr << completed_divisions << " " << desired_divisions << endl;

  if (output_parameters.verbose) {
    cerr << "Completed divisions: " << completed_divisions << '\n';
    for (auto& current_subpopulation : current_subpopulations) {
      cerr << "Genotype: " << current_subpopulation.GetNode_id() << " Frequency: " << current_subpopulation.GetNumber()
           << " Fitness: " << current_subpopulation.GetFitness() << '\n';
    }
  }

  num_divisions_until_next_mutation -= completed_divisions;

  /*
   if (num_divisions_until_next_mutation <= 0)
   if (completed_divisions > 3)
   cerr << completed_divisions << endl;
   */ // We should really backtrack so that all pops divided simultaneously, this is slightly off.
  // We could move time forward by the number of mutations and then weight which subpop it occurs in by the number of
  // divisions in each one

  if (output_parameters.verbose) {
    cerr << "Divisions until mutation: " << num_divisions_until_next_mutation << '\n';
  }
}

void cPopulation::TransferResampleDilution()
{
  // When it is time for a transfer, resample population

  // @JEB: Move to main simulation loop
  num_completed_transfers++;

  if (output_parameters.verbose) cerr << "  Population size before transfer: " << current_population_size << '\n';

  current_population_size = 0;  // recalculate population size
  for (auto it = current_subpopulations.begin(); it != current_subpopulations.end(); ++it) {

    // Perform binomial sampling
    it->Transfer(simulation_parameters.binomial_sampling_transfer_probability, rng);

    //@agm This section is to delete new mutations from the tree that do not get passed,
    //     it also deletes subpopulations from the list that have zero population

    current_population_size += it->GetNumber();

    if (it->GetNumber() == 0) {
      if (new_genotype_ids_since_last_transfer.count(it->GetNode_id())) {
        it->GetGenotype().marked_for_deletion = true;
      }
      it = current_subpopulations.erase(it);
      --it;
      replicate_statistics.total_subpopulations_lost++;
      existing_genotype_count--;
    }
  }

  if (output_parameters.verbose) cerr << "  Population size after transfer: " << current_population_size << '\n';

  // Check our calculated population size
  // Check how we have updated the population size
  if (debug) {
    double calculated_population_size = CalculatePopulationSize();
    assert(calculated_population_size == current_population_size);
  }
}

void cPopulation::SeedPopulationWithMarkedGenotypes(uint32_t genetic_marker_states)
{
  // For now assumes...
  assert(genetic_marker_states == 2);

  cGenotype r(simulation_parameters, unique_genotype_count++);
  cGenotype w(simulation_parameters, unique_genotype_count++);
  tree<cGenotype>::iterator top, red_side, white_side;
  double starting_fitness = simulation_parameters.initial_fitness;

  // initialize object of cSubpopulation type
  cSubpopulation red, white;

  /*@agm The unique code and fitness is set. */

  /*@jeb note that we increment genotype count here.
   It should really happen in AddSubpopulation?
   */

  r.fitness = starting_fitness;
  w.fitness = starting_fitness;

  // Start building the tree
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

// update: Victory conditions are now a command line option.  If not set to false the simulation
//         will run to the max number of iterations.

void cPopulation::SeedPopulationWithOneGenotype()
{
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

  if (output_parameters.verbose) {
    cout << "Completed divisions: " << num_completed_divisions << '\n';
    for (auto& current_subpopulation : current_subpopulations) {
      cout << "Genotype1: " << current_subpopulation.GetNode_id()
           << " Frequency1: " << current_subpopulation.GetNumber() << '\n';
    }
  }
}

// Note that if you call add subpopulation, you may need to update the population size!
void cPopulation::AddSubpopulation(cSubpopulation& subpop)
{
  genotype_id_map_into_tree[subpop.GetGenotype().unique_node_id] = subpop.GetGenotypeIter();
  current_subpopulations.push_back(subpop);
  existing_genotype_count++;

  new_genotype_ids_since_last_transfer.insert(subpop.GetNode_id());  // Keep track of all mutations since last transfer
  // cout << subpop.GetNode_id() << " " << subpop.GetFitness() << endl;
}

// Don't cull all populations that disappeared because we need the final tree of all genotypes
void cPopulation::CullSubpopulationsThatDidNotEstablish()
{

  // Be sure to use leaf_begin() and leaf_end()

  for (tree<cGenotype>::leaf_iterator it = genotype_tree.begin_leaf(); it != genotype_tree.end_leaf();) {
    if (it->marked_for_deletion) {
      tree<cGenotype>::leaf_iterator ite = it;
      it++;
      genotype_id_map_into_tree.erase(ite->unique_node_id);
      genotype_tree.erase(ite);
    } else {
      it++;
    }
  }
  new_genotype_ids_since_last_transfer.clear();
}

// Generates new mutant and adds it to the tree
void cPopulation::MutateExactWithFractionalCells()
{
  if (output_parameters.verbose) cout << "* Mutating!" << '\n';
  if (output_parameters.verbose) cout << "Total population: " << CalculatePopulationSize() << '\n';

  // Mutation depends on instantaneuous Fitness * population size of subpopulations.
  vector<double> fractional_chances_per_subpopulation(current_subpopulations.size());

  // All subpopulations can mutate as long as they have more than one cell
  // (so they won't go negative). Create table of their chances.
  double total = 0;
  for (uint32_t i = 0; i < current_subpopulations.size(); i++) {
    fractional_chances_per_subpopulation[i] =
        (current_subpopulations[i].GetNumber() > 1)
            ? current_subpopulations[i].GetFitness() * current_subpopulations[i].GetNumber()
            : 0;
    total += fractional_chances_per_subpopulation[i];
  }

  for (double& i : fractional_chances_per_subpopulation) {
    i /= total;
  }

  // Pick which one divided.
  // Gets a random number in range [0,1)
  double random_fraction = rng.Uniform();
  uint32_t chosen_subpopulation = 0;
  while (random_fraction > fractional_chances_per_subpopulation[chosen_subpopulation]) {
    random_fraction -= fractional_chances_per_subpopulation[chosen_subpopulation];
    chosen_subpopulation++;
  }

  cSubpopulation& ancestor = current_subpopulations[chosen_subpopulation];
  cSubpopulation new_subpop;

  // There must be at least two cells for a mutation to have occurred...
  if (ancestor.GetNumber() <= 1) {
    cout << "Died: " << ancestor.GetNumber() << " " << chosen_subpopulation << " / " << current_subpopulations.size()
         << " " << random_fraction << '\n';
  }
  assert(ancestor.GetNumber() > 1.0);

  bool success_adding = new_subpop.CreateDescendant(rng, ancestor, simulation_parameters, genotype_tree,
                                                    unique_genotype_count++, output_parameters.verbose);

  if (success_adding) {

    if (output_parameters.verbose) cout << "  Color: " << new_subpop.GetMarker() << '\n';
    if (output_parameters.verbose) cout << "  New Fitness: " << new_subpop.GetFitness() << '\n';

    AddSubpopulation(new_subpop);

    // Update the maximum subpopulation fitness
    maximum_subpopulation_fitness = max(maximum_subpopulation_fitness, new_subpop.GetFitness());

    if (output_parameters.verbose) cout << "  New Number: " << new_subpop.GetNumber() << '\n';
    if (output_parameters.verbose) cout << "  Old Number: " << ancestor.GetNumber() << '\n';

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
  if (num_completed_transfers < 0) return;

  // Only record if we are on the coarse-graining interval
  if (num_completed_transfers % output_parameters.coarse_graining != 0) return;

  replicate_statistics.average_population_fitness.push_back(average_subpopulation_fitness);

  //  Debug: this checks to see if our population size was correct.
  if (debug) {
    double calculated_population_size = CalculatePopulationSize();
    assert(current_population_size == calculated_population_size);
  }

  /////////////////////////////////////////
  /// Statistic: mutation_counts
  /////////////////////////////////////////

  double average_total_mutation_count = 0;
  vector<double> average_mutation_counts(simulation_parameters.mutation_rates_per_division.size(), 0);

  for (auto& this_pop : current_subpopulations) {
    for (uint32_t on_category = 0; on_category < average_mutation_counts.size(); ++on_category) {
      average_mutation_counts[on_category] +=
          this_pop.GetGenotype().mutation_counts[on_category] * this_pop.GetNumber();
    }
    average_total_mutation_count += this_pop.GetGenotype().total_mutation_count * this_pop.GetNumber();
  }

  for (double& average_mutation_count : average_mutation_counts) {
    average_mutation_count /= current_population_size;
  }
  average_total_mutation_count /= current_population_size;

  replicate_statistics.average_mutation_counts.push_back(average_mutation_counts);
  replicate_statistics.average_total_mutation_count.push_back(average_total_mutation_count);

  /////////////////////////////////////////
  /// Statistic: genotype_frequencies,
  /// Statistic: clade_frequencies
  /////////////////////////////////////////

  // Initialize all entries
  GenotypeFrequencyMap current_genotype_frequency_map;  // only of exactly this genotype
  GenotypeFrequencyMap current_clade_frequency_map;     // AND including all descendants

  tree<cGenotype>::iterator update_location, genotype;
  for (auto& this_pop : current_subpopulations) {

    update_location = this_pop.GetGenotypeIter();
    assert(update_location != nullptr);

    double& current_genotype_frequency = current_genotype_frequency_map[update_location->unique_node_id];
    current_genotype_frequency = this_pop.GetNumber() / current_population_size;

    // Add frequency to all parents to get clade numbers
    while (update_location != nullptr) {

      double& current_clade_frequency = current_clade_frequency_map[update_location->unique_node_id];
      current_clade_frequency += this_pop.GetNumber() / current_population_size;

      // current_clade_frequencies[update_location->unique_node_id].unique_node_id = update_location->unique_node_id;
      // current_clade_frequencies[update_location->unique_node_id].m_frequency += floor(this_pop.GetNumber());

      update_location =
          tree<bpopsim::cGenotype, std::allocator<tree_node_<bpopsim::cGenotype>>>::parent(update_location);
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
    // set<uint32_t> dominant_clade_set = GenotypesFromAncestorToFinalSweep();
    set<uint32_t> dominant_clade_set = GenotypesFromAncestorToFinalDominant();
    // cerr << "Number in dominant clade: " << dominant_clade_set.size() << endl;

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
      for (GenotypeFrequencyMap::iterator it=this_genotype_frequencies.begin(); it!= this_genotype_frequencies.end();
    it++) { uint32_t this_node_id = it->first;

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
    for (auto& this_clade_frequencies : replicate_statistics.clade_frequencies) {
      vector<double> add_frequencies(max_depth, 0);

      // Determine the number of steps back to the line of descent and add to this category and higher ones
      for (auto it = this_clade_frequencies.begin(); it != this_clade_frequencies.end(); it++) {
        uint32_t this_node_id = it->first;

        tree<cGenotype>::iterator itt = FindGenotypeInTreeByID(this_node_id);

        uint32_t this_depth = 0;
        while ((itt != nullptr) && (!dominant_clade_set.count(itt->unique_node_id))) {

          if (this_depth < max_depth) {
            add_frequencies[this_depth] = max(add_frequencies[this_depth], this_clade_frequencies[this_node_id]);
          }

          itt = tree<bpopsim::cGenotype, std::allocator<tree_node_<bpopsim::cGenotype>>>::parent(itt);
          this_depth++;
        }

        assert(itt != nullptr);
      }  // end genotype loop

      replicate_statistics.diverged_frequencies_by_depth.push_back(add_frequencies);
    }  // end transfer loop
  }

  if (output_parameters.output_largest_sweep_size) {
    assert(false);  // not implemented
  }

  if (output_parameters.output_longest_sweep_time) {
    bool debug = false;

    set<uint32_t> dominant_clade_set = GenotypesFromAncestorToFinalDominant();
    replicate_statistics.longest_sweep_time = 0;
    // for each genotype on the dominant line of descent
    for (unsigned int on_genotype_id : dominant_clade_set) {

      int32_t start_transfer = 0;
      int32_t end_transfer = 0;

      // determine how many transfers it took to go from ≤5% to ≥95%
      if (debug) cout << "== " << on_genotype_id << " ==" << '\n';
      for (uint32_t transfer = 0; transfer < replicate_statistics.clade_frequencies.size(); transfer++) {
        GenotypeFrequencyMap& this_clade_frequencies = replicate_statistics.clade_frequencies[transfer];
        if (this_clade_frequencies[on_genotype_id] <= 0.05) {
          start_transfer = transfer;
          if (debug) cout << this_clade_frequencies[on_genotype_id] << " (start)" << '\n';
        } else if (this_clade_frequencies[on_genotype_id] >= 0.95) {
          end_transfer = transfer;
          if (debug) cout << this_clade_frequencies[on_genotype_id] << " (end)" << '\n';
          break;
        } else {
          if (debug) cout << this_clade_frequencies[on_genotype_id] << '\n';
        }
      }

      if (debug) cout << start_transfer << "-" << end_transfer << '\n';

      int32_t this_sweep_time = end_transfer - start_transfer;
      replicate_statistics.longest_sweep_time =
          max(replicate_statistics.longest_sweep_time,
              this_sweep_time * static_cast<int32_t>(output_parameters.coarse_graining));
    }
  }  // end output_longest_sweep_time
}

//////////  Utility Functions //////////

// Calculates the total population size by iterating over subpops.
// This is really for debug checking of code only. Normally, use current_population_size.
double cPopulation::CalculatePopulationSize()
{
  double calculated_population_size = 0;
  for (auto& current_subpopulation : current_subpopulations) {
    calculated_population_size += current_subpopulation.GetNumber();
  }
  return calculated_population_size;
}

// Calculates the maximum population fitness by iterating over subpops.
// This is really for debug checking of code only. Normally, use maximum_subpopulation_fitness.

double cPopulation::CalculateMaximumSubpopulationFitness()
{
  double max_fitness = 0.0;
  for (auto& current_subpopulation : current_subpopulations) {
    max_fitness = max(max_fitness, current_subpopulation.GetFitness());
  }
  return max_fitness;
}

double cPopulation::CalculateAverageSubpopulationFitness()
{
  double avg_fitness = 0.0;
  for (auto& current_subpopulation : current_subpopulations) {
    avg_fitness += current_subpopulation.GetNumber() * current_subpopulation.GetFitness();
  }
  return avg_fitness / current_population_size;
}

// Helper function for producing Muller plot matrix
// Assigns bounds to child genotype centered (and equally dividing) parent genotype slice.
double cPopulation::AssignChildFrequency(const tree<cGenotype>::sibling_iterator& this_node, double in_low,
                                         double /*in_high*/, vector<cFrequencySlice>* child_freqs,
                                         GenotypeFrequencyMap& frequencies,
                                         int depth)  // current depth in tree, defaults to zero
{

  // kptree::print_tree_bracketed(*newtree);

  // The low for this mutation should be the low of the input interval
  double this_low = in_low;
  double this_high = this_low + frequencies.GetFrequency(this_node->unique_node_id);
  double size_depth1_children(0), half_size_parent_swath((this_high - this_low) / 2);

  for (tree<cGenotype>::sibling_iterator it_node = genotype_tree.begin(this_node);
       it_node != genotype_tree.end(this_node); ++it_node) {
    size_depth1_children += frequencies.GetFrequency(it_node->unique_node_id);
  }
  /*if(size_depth1_children != 0)
   cout << size_depth1_children << endl;*/

  double this_bottom_high(this_low + half_size_parent_swath - (size_depth1_children / 2));
  double this_top_low(this_bottom_high);

  // The swath for this mutation may shrink on the bottom
  // due to its children taking a bite out of it.
  double last_assigned_child_high = this_top_low;
  for (tree<cGenotype>::sibling_iterator it_node = genotype_tree.begin(this_node);
       it_node != genotype_tree.end(this_node); ++it_node) {

    // is the frequency > 0? (It may have gone extinct, in which case it is a waste to keep going down the tree!)
    if (frequencies.GetFrequency(it_node->unique_node_id) > 0) {
      last_assigned_child_high =
          AssignChildFrequency(it_node, this_top_low, this_high, child_freqs, frequencies, depth + 1);
    }

    // Our new low is the last high assigned to a child
    this_top_low = last_assigned_child_high;
  }

  // At this point we know the top and bottom of this node...
  // This is where we would want to paint into a bitmap between them!!
  // Draw between this_low and this_high ... rounding to only paint whole numbers ...

  (*child_freqs).emplace_back(this_node->unique_node_id, this_low, this_bottom_high);
  (*child_freqs).emplace_back(this_node->unique_node_id, this_top_low, this_high);

  // print indented version
  /*for (int i=0; i<depth; i++) {
   cout << " ";
   }*/
  // cout << " ID:" << this_node->unique_node_id << " Freq:" << (*frequencies)[this_node->unique_node_id].frequency;
  // cout << " [" << this_low << "," << this_high << "]" << endl;

  return this_high;
}

tree<cGenotype>::iterator cPopulation::FindGenotypeInTreeByID(uint32_t id)
{
  return genotype_id_map_into_tree[id];
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

  while (tree<bpopsim::cGenotype, std::allocator<tree_node_<bpopsim::cGenotype>>>::number_of_children(parent_node) >
         0) {

    double biggest_child_frequency = -1.0;
    tree<cGenotype>::iterator biggest_child = nullptr;
    for (tree<cGenotype>::sibling_iterator it_node = genotype_tree.begin(parent_node);
         it_node != genotype_tree.end(parent_node); ++it_node) {

      double this_frequency = frequencies.GetFrequency(it_node->unique_node_id);
      if (this_frequency > biggest_child_frequency) {
        biggest_child = it_node;
        biggest_child_frequency = this_frequency;
      }
    }

    parent_node = biggest_child;
    genotype_set.insert(parent_node->unique_node_id);
  }

  // cerr << genotype_set.size() << endl;

  return genotype_set;
}

// Returns a list of all clades that were ever above the
// requested threshold at any transfer during the simulation.

vector<uint32_t> cPopulation::CladesAboveThreshold(float threshold)
{

  set<uint32_t> genotype_set;

  for (auto& transfer_clade_frequencies : replicate_statistics.clade_frequencies) {
    for (auto& transfer_clade_frequencie : transfer_clade_frequencies) {
      if (transfer_clade_frequencie.second >= threshold) genotype_set.insert(transfer_clade_frequencie.first);
    }
  }

  vector<uint32_t> genotype_list(genotype_set.begin(), genotype_set.end());
  sort(genotype_list.begin(), genotype_list.end());
  return genotype_list;
}

// Returns a list of all genotypes that were ever above the
// requested threshold at any transfer during the simulation.

vector<uint32_t> cPopulation::GenotypesAboveThreshold(float threshold)
{

  set<uint32_t> genotype_set;

  for (uint32_t transfer = 0; transfer < replicate_statistics.genotype_frequencies.size(); transfer++) {
    GenotypeFrequencyMap& transfer_genotype_frequencies = replicate_statistics.clade_frequencies[transfer];

    for (auto& transfer_genotype_frequencie : transfer_genotype_frequencies) {
      if (transfer_genotype_frequencie.second >= threshold) genotype_set.insert(transfer_genotype_frequencie.first);
    }
  }

  vector<uint32_t> genotype_list(genotype_set.begin(), genotype_set.end());
  sort(genotype_list.begin(), genotype_list.end());
  return genotype_list;
}
