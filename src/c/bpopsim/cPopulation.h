#pragma once

#include "common.h"
#include "cSubpopulation.h"

namespace bpopsim {

  class GenotypeFrequencyMap : public map<uint32_t, double>
  {

  public:
    // Returns zero if not found
    double GetFrequency(uint32_t in_genotype_id)
    {
      if (!count(in_genotype_id)) return 0.0;

      return (*this)[in_genotype_id];
    }
  };

  // Statistics for a single run
  struct cReplicateStatistics
  {

    vector<double> average_population_fitness;
    vector<double> average_total_mutation_count;     // average number of total mutations
    vector<vector<double>> average_mutation_counts;  // average number of mutations by category

    vector<GenotypeFrequencyMap> genotype_frequencies;  // Frequencies of each genotype, map by genotype_id
    vector<GenotypeFrequencyMap>
        clade_frequencies;  // Frequencies of each clade with the given genotype as its ancestor, map by genotype_id

    vector<vector<double>> diverged_frequencies_by_depth;  // percent of population diverged from dominant lineage by
                                                           // certain depth (number of mutations) or more
    int32_t longest_sweep_time = 0;                        // for any mutation in the replicate; measured in transfers

    vector<vector<double>> marker_frequencies;  // [recorded transfer][marker - 1]: fraction of the population
                                                // carrying each neutral marker state (marker mode only)

    uint32_t total_mutations = 0;            // Number of mutations so far
    uint32_t total_subpopulations_lost = 0;  // Number of subpopulations that no longer exist
  };

  // Options controlling what is recorded and written, plus run-time switches.
  // Plain data with defaults matching the command-line defaults in options.cc;
  // filled from the command line by ParseOutputParameters().
  struct OutputParameters
  {
    string output_directory_name = ".";
    uint32_t coarse_graining = 1;
    uint32_t burn_in = 0;  // Number of transfers to perform before recording output
    uint32_t diverged_mutation_depth = 10;
    bool output_diverged_frequencies = false;
    bool output_longest_sweep_time = false;
    bool output_marker_frequencies = false;
    double minimum_output_frequency = 0.01;
    bool verbose = false;  // Print step-by-step details of the run (-v)
    bool debug = false;    // Perform extra consistency checks (--debug)
  };

  // Options for marker-divergence mode (-k). Filled by ParseMarkerDivergence().
  struct MarkerDivergence
  {
    uint32_t marker_states = 0;          // Number of neutral marker states the population starts with (0 = off)
    double divergence_end_fraction = 0;  // End a replicate once one marker reaches this fraction (0 = off)
  };

  // Statistics across all runs
  struct cStatistics : vector<cReplicateStatistics>
  {

    explicit cStatistics(const OutputParameters& output_parameters)
    : output_directory_name(output_parameters.output_directory_name)
    , coarse_graining(output_parameters.coarse_graining)
    {}

    void OutputAveragePopulationFitness();
    void OutputAveragePopulationMutationCounts();
    void OutputDivergedFrequenciesByDepth();
    void OutputLongestSweepTime();
    void OutputMarkerFrequencies();

    string output_directory_name;
    uint32_t coarse_graining;  // Only used for calculating statistics and ignoring intermediate time points
  };

  // Parameters of the evolutionary model. Plain data with defaults matching the
  // command-line defaults in options.cc. Set the input fields, then call
  // Finalize() to compute the derived fields (ParseSimulationParameters() does
  // both from a parsed command line).
  class SimulationParameters
  {

  public:
    SimulationParameters() { Finalize(); }

    //! Throws std::invalid_argument if the inputs are inconsistent.
    void Validate() const
    {
      const size_t categories = mutation_rates_per_division.size();
      if (categories == 0) throw std::invalid_argument("At least one mutation category (-u) is required.");
      if (mutation_fitness_effects.size() != categories || mutation_fitness_effect_models.size() != categories) {
        throw std::invalid_argument(
            "The -u, -s and -f lists must have the same shape (same '|' groups and ',' items).");
      }
      size_t grouped = 0;
      for (uint32_t size : mutation_category_group_sizes) grouped += size;
      if (grouped != categories) {
        throw std::invalid_argument("Mutation category group sizes do not add up to the number of categories.");
      }
      for (double rate : mutation_rates_per_division) {
        if (rate < 0) throw std::invalid_argument("Mutation rates (-u) must be 0 or greater.");
      }
      for (char model : mutation_fitness_effect_models) {
        if (model != 'u' && model != 'e' && model != 'o') {
          throw std::invalid_argument(std::string("Unknown fitness effect model (-f) '") + model +
                                      "'. Choices are u, e, o.");
        }
      }
    }

    //! Compute the derived fields from the input fields. Call after changing any input.
    void Finalize()
    {
      // Hand-built parameters may set only the rates and effects: broadcast a single
      // model to every category (as the command line does) and repair the display
      // grouping if it no longer matches.
      const size_t categories = mutation_rates_per_division.size();
      if (mutation_fitness_effect_models.size() == 1 && categories > 1) {
        const char model = mutation_fitness_effect_models[0];  // copy: assign() must not alias an element
        mutation_fitness_effect_models.assign(categories, model);
      }
      size_t grouped = 0;
      for (uint32_t size : mutation_category_group_sizes) grouped += size;
      if (grouped != categories) mutation_category_group_sizes = {static_cast<uint32_t>(categories)};

      transfer_dilution_factor = exp(log(2) * generations_per_transfer);
      final_population_size_at_transfer = transfer_dilution_factor * initial_population_size_after_transfer;
      binomial_sampling_transfer_probability = 1.0 / transfer_dilution_factor;

      // Sum mutation rates
      total_mutation_rate_per_division = 0;
      for (double& it : mutation_rates_per_division) total_mutation_rate_per_division += it;

      total_mutation_rate_exponential_mean = 1.0 / total_mutation_rate_per_division;

      // Fill in the fractional table for deciding which category mutations belong in.
      fractional_chances_of_mutation_categories.clear();
      for (double& it : mutation_rates_per_division) {
        fractional_chances_of_mutation_categories.push_back(it / total_mutation_rate_per_division);
      }
    }

    void Print() const
    {
      std::cerr << "***** Simulation Settings *****" << '\n';
      std::cerr << "  Initial population size                  (-p) = " << initial_population_size << '\n';
      std::cerr << "  Generations of growth per transfer       (-T) = " << generations_per_transfer << '\n';
      std::cerr << "  Transfer dilution                             = " << transfer_dilution_factor << '\n';
      std::cerr << "  Population size after transfer           (-N) = " << initial_population_size_after_transfer
                << '\n';
      std::cerr << "  Population size at transfer                   = " << final_population_size_at_transfer << '\n';
      std::cerr << "  Mutation rates per division              (-u) = " << GroupedString(mutation_rates_per_division)
                << '\n';
      std::cerr << "  Mutation fitness effects                 (-s) = " << GroupedString(mutation_fitness_effects)
                << '\n';
      std::cerr << "  Mutation fitness effect models           (-f) = " << GroupedString(mutation_fitness_effect_models)
                << '\n';
      std::cerr << "  Initial mutation fitness effects         (-1) = "
                << substitute(first_mutation_fitness_effects_string, "\n", ",") << '\n';
      std::cerr << "  Average fitness end condition (0=OFF)         = " << average_fitness_end_condition << '\n';
    }

    uint32_t DetermineMutationCategory(RandomGenerator& rng)
    {
      // Don't draw a random number if there is only one category
      if (fractional_chances_of_mutation_categories.empty()) return 0;

      double random_fraction = rng.Uniform();
      const uint32_t last_category = fractional_chances_of_mutation_categories.size() - 1;
      uint32_t category = 0;
      while (category < last_category && random_fraction > fractional_chances_of_mutation_categories[category]) {
        random_fraction -= fractional_chances_of_mutation_categories[category];
        category++;
      }
      // The fractions sum to 1 up to rounding; the loop bound guarantees the
      // last category absorbs any residue instead of indexing past the end.
      assert(category < fractional_chances_of_mutation_categories.size());
      return category;
    }

    //! Render a per-category list in the command-line notation: ',' within a
    //! '|' group, '|' between groups (see mutation_category_group_sizes).
    template <class T> std::string GroupedString(const vector<T>& per_category) const
    {
      std::string out;
      size_t index = 0;
      for (size_t group = 0; group < mutation_category_group_sizes.size(); ++group) {
        if (group > 0) out += '|';
        for (uint32_t i = 0; i < mutation_category_group_sizes[group] && index < per_category.size(); ++i, ++index) {
          if (i > 0) out += ',';
          out += to_string(per_category[index]);
        }
      }
      return out;
    }

    // Input
    double maximum_number_of_transfers =
        50;  // Maximum number of transfers before ending simulation. May end earlier for other reasons.
    double initial_population_size = 5E6;                 // Initial population size at the beginning of the simulation
    double initial_population_size_after_transfer = 5E6;  // Initial population size after transfer dilution
    double generations_per_transfer = 6.64;               // Number of binary cell divisions per transfer

    double initial_fitness = 1.0;                 // Fitness of cells at the beginning of the simulation
    bool multiplicative_fitness_effects = false;  // Multiply rather than add selection coefficients *(1+s) versus +s
    // Mutation categories are flat lists indexed by category. On the command line
    // they are written in '|'-separated groups of ','-separated items;
    // mutation_category_group_sizes records that grouping for display only.
    vector<double> mutation_rates_per_division = {1E-7};   // Rate per cell division of mutations in each category
    vector<double> mutation_fitness_effects = {0.1};       // Effect (or mean effect) of a mutation in each category
    vector<char> mutation_fitness_effect_models = {'u'};   // Per category: 'u' fixed, 'e' exponential, 'o' one-time
    vector<uint32_t> mutation_category_group_sizes = {1};  // Number of categories in each '|' group
    vector<double> first_mutation_fitness_effects;  // Fitness effects of initial mutations (overrides other options)
    string first_mutation_fitness_effects_string;   //   string version of previous for printing
    double average_fitness_end_condition = 0;       // End replicate when average fitness is this or greater (0 = off)

    // Calculated from Input by Finalize()
    double transfer_dilution_factor;                           // The transfer dilution fraction
    double final_population_size_at_transfer;                  // Final population size that triggers a transfer
    double binomial_sampling_transfer_probability;             // Probability for binomial sampling
    double total_mutation_rate_exponential_mean;               // Mean for Poisson draws of divisions until next mut
    double total_mutation_rate_per_division;                   // Sum of all mutation rates
    vector<double> fractional_chances_of_mutation_categories;  // Fraction of mutations in each category
  };

  class cPopulation
  {

  private:
    SimulationParameters simulation_parameters;
    OutputParameters output_parameters;
    MarkerDivergence marker_divergence;

  public:
    cReplicateStatistics replicate_statistics;

  private:
    // Internal state variables
    bool debug;
    uint32_t replicate;

    // Tracks characteristics of all subpopulations (including location in Tree)
    vector<cSubpopulation> current_subpopulations;

    // Tree that tracks ancestry of the subpopulations
    tree<cGenotype> genotype_tree;
    uint32_t unique_genotype_count{
        0};  // used to assign node ids in tree, should equal number of nodes that ever existed
    map<uint32_t, tree<cGenotype>::iterator> genotype_id_map_into_tree;
    uint32_t existing_genotype_count{0};  // only the current number of genotypes in existence
    double maximum_subpopulation_fitness;

    double current_population_size{0};  // Current population size in terms of whole cells

    // New subpopulations that arose this transfer. Allows culling subpopulations did not survive a transfer.
    set<uint32_t> new_genotype_ids_since_last_transfer;

    // All of the following should be initialized in the constructor
    int32_t
        num_completed_transfers;  // Number of transfers that have been completed thus far, can be negative for burn-in
    double average_subpopulation_fitness{0};
    bool run_end_condition_met{false};

    // @JEB: An double rather than a uint because this can go negative by a few cells
    //       when cells (usually the ancestors) divide simultaneously.
    double num_divisions_until_next_mutation{0};
    double num_completed_divisions{0};

    RandomGenerator& rng;
    const double log_2;

  public:
    // CONSTRUCTOR
    cPopulation(SimulationParameters in_simulation_parameters, const OutputParameters& in_output_parameters,
                const MarkerDivergence& in_marker_divergence, RandomGenerator& in_rng, uint32_t in_replicate);

    // DESTRUCTOR
    virtual ~cPopulation() = default;

    void RunSimulation();

    void OutputCladeFrequencies();
    void OutputGenotypeFrequencies();
    void OutputMullerMatrix(uint32_t frequency_resolution, bool line_of_descent_separate);

  private:
    //!!! Private Simulation Helper Methods

    //! Move time forward by this increment, growing all subpopulations
    void UpdateSubpopulationsForGrowthExactWithFractionalCells(double update_time);

    //! Calculate the time until the next subpopulation divides (passes a whole number of cells)
    double TimeToNextWholeCell();

    //! Calculate an amount of time to divide that won't take us over the time until the next mutation
    double CalculateDivisionsUntilNextBeneficialMutation()
    {
      return round(rng.Exponential(simulation_parameters.total_mutation_rate_exponential_mean));
    }
    void ProcessCellDivisionTimeStepExactWithFractionalCells();

    //! Population transfer: binomial dilution of every subpopulation
    void TransferResampleDilution();

    //! Methods for initializing the population
    void SeedPopulationWithMarkedGenotypes(uint32_t genetic_markers);
    void SeedPopulationWithOneGenotype();

    //! Methods for introducing mutants into the population
    void MutateExactWithFractionalCells();

    // General function used whenever a new subpopulation is added
    void AddSubpopulation(cSubpopulation& subpop);

    // Used to remove genotypes that died before a transfer to save memory
    void CullSubpopulationsThatDidNotEstablish();

    // Methods for recording statistics during a run
    void RecordStatisticsAtTransfer();
    void RecordStatisticsAtEnd();

    //!!! Utility calculations - mainly used in debug mode

    double CalculatePopulationSize();
    double CalculateMaximumSubpopulationFitness();
    double CalculateAverageSubpopulationFitness();
    //! Fraction of the population carrying each marker state, index = marker - 1 (marker mode only)
    vector<double> CalculateMarkerFrequencies();

    double AssignChildFrequency(const tree<cGenotype>::sibling_iterator& child_node, double parent_low,
                                double parent_high, vector<cFrequencySlice>* child_freqs,
                                GenotypeFrequencyMap& frequencies, int depth = 0);

    //!!! Tree queries used by the output methods

    tree<cGenotype>::iterator FindGenotypeInTreeByID(uint32_t id);
    //! Genotype ids from the ancestor down the line of descent to the final dominant genotype
    set<uint32_t> GenotypesFromAncestorToFinalDominant();
    //! Ids of genotypes / clades that were ever at or above the threshold frequency
    vector<uint32_t> GenotypesAboveThreshold(float threshold);
    vector<uint32_t> CladesAboveThreshold(float threshold);
  };

}  // namespace bpopsim
