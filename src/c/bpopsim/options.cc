#include "options.h"

namespace bpopsim {

  void DefineBpopsimOptions(OptionParser& options)
  {
    // @JEB: doubles must be put in quotes if in scientific notation
    options("help,h", "Produce this help message", NoArgument);
    options("verbose,v", "Output extra details of runs", NoArgument);
    options("debug", "Perform extra checks on simulation", NoArgument);

    options.AddUsage("");
    options.AddUsage("==== Simulation Options ====");
    options.AddUsage("");
    options("rng-seed,d", "Seed for random number generator. [0 = choose a random seed and report it]", 0);
    options("replicates,r", "Replicates", 1);
    options("initial-population-size,p", "Initial population size.", "5E6");
    options("generations-per-transfer,T", "Generations per transfer", 6.64);
    options("population-size-after-transfer,N", "Population size after transfer", "5E6");
    options("number-of-transfers,n", "Max number of transfer to replicate", 50);
    options("initial-fitness,z", "Initial fitness.", 1.0);
    options("multiplicative-fitness,m", "Use multiplicative fitness model [f * (1 + s)]. Default = additive [f + s].",
            NoArgument);
    options("mutation-rates,u",
            "Beneficial mutation rate per cell division for each category of mutations. Categories are written as "
            "','-separated items within '|'-separated groups, e.g. \"1E-7,1E-8|1E-9\"; repeating the option adds "
            "another group. A mutation occurs at the total rate and is assigned to a category in proportion to its "
            "rate.",
            "1E-7");
    options("fitness-effects,s",
            "Fitness effect (or mean effect) of a mutation in each category. Same '|' and ',' shape as -u.", 0.1);
    options("fitness-effect-model,f",
            "How the fitness effect is drawn in each category: 'u' every mutation has exactly the -s value, 'e' "
            "exponentially distributed with mean -s, 'o' exactly the -s value but only one mutation per category per "
            "lineage. Same '|' and ',' shape as -u, or a single letter applied to every category.",
            'u');
    options("first-mutation-fitness-effects,1",
            "Use these fitness effects for the first mutations. Supply option multiple times to define the fitness "
            "effects of the 1st, 2nd, 3rd, etc. mutations. These values will override the -s option until the number "
            "of mutations exceeds the number of these values.");
    options("average-fitness-end-condition",
            "End a replicate when the average fitness of the population is this or greater (0=OFF)", 0);

    options.AddUsage("");
    options.AddUsage("==== Marker Divergence Options ====");
    options.AddUsage("");
    options("marker-states,k",
            "Begin with the initial population divided equally into this many genotypes, each carrying a different "
            "state of a neutral genetic marker that its descendants inherit. Each marker starts with "
            "initial-population-size / N cells. Must be 0 (OFF) or at least 2.",
            0);
    options("marker-divergence-end-fraction,x",
            "End a replicate once the most common marker state makes up at least this fraction of the population "
            "(0=OFF). Requires --marker-states.",
            0);
    options("output-marker-frequencies",
            "Output the fraction of the population carrying each marker state at each recorded transfer to "
            "'marker_frequencies.tab'. Requires --marker-states.",
            NoArgument);
    options.AddUsage("");
    options.AddUsage("==== Statistical Output Files and Options ====");
    options.AddUsage("");
    options("output-folder,o", "Base output folder. All output files will be created here.", ".");
    options("burn-in,b",
            "Perform this number of transfers before beginning to output statistics (before transfer numbered 0).", 0);
    options("coarse-graining,c", "Only print stats every this many transfers to all output files.", 1);
    options("output-average-fitness", "Output average fitness to file 'average_fitness.tab'", NoArgument);
    options("output-average-mutations", "Output mutation numbers to file 'average_mutations_count.tab'", NoArgument);
    options("output-clade-frequencies",
            "Output clade (mutation) frequencies (above minimum) for replicate X to file 'clade_frequencies_X.tab'",
            NoArgument);
    options("output-genotype-frequencies",
            "Output genotype frequencies (above minimum) for replicate X to file 'genotype_frequencies_X.tab'.",
            NoArgument);
    options("minimum-output-frequency",
            "Only output clade/genotype frequencies that are above this frequency at some time point.", 0.01);
    options("output-diverged-frequencies",
            "Output frequencies of clades that differ by at least X mutations from the line of descent to the final "
            "dominant to 'diverged_frequencies_X.tab'",
            NoArgument);
    options("diverged-mutation-depth", "Maximum depth of diverged mutations to output.", 10);
    options("output-longest-sweep-time",
            "Output the longest time that a mutation along the line of descent to the final dominant took to sweep "
            "from ≥5% to ≥95% frequency in the population to 'longest_sweep_time.tab'",
            NoArgument);
    options.AddUsage("");
    options.AddUsage("==== Muller Plot Output Files and Options ====");
    options.AddUsage("");
    options("output-muller",
            "Output Muller matrix for replicate X to file 'muller_matrix_X.dat'. Cell values unique genotypes. This "
            "file can be used to make plots of population dynamics.",
            NoArgument);
    options("muller-resolution", "Muller matrix vertical resolution.", 200);
    options("muller-lod",
            "In the Muller matrix, give genotypes on the line of descent to the final dominant negative numbers, so "
            "that they can be plotted in separate colors.",
            NoArgument);

    // Not (re)implemented
    //  options("output-dominant-genotypes", "Output dominant genotype frequencies to file 'dominant_genotypes.csv'",
    //  NoArgument);
    //  ("output-time-to-sweep", "Print time to sweep for each mutation", NoArgument)
    //  ("output-max-diff", "Print max difference of sweeping mutations", NoArgument)
    //  ("output-single-fitness", "Print the fitness of the single cell.", NoArgument)
    //  ("output-sweeping-descent-fitness", "Output the average fitness of sweeping descent per time.",
    //  NoArgument)

    // options.AddUsage("");
    // options.AddUsage("==== Other Options ====");
    // options.AddUsage("");
    // options("exact-mutations-per-transfer,y", "Cause exactly this number of mutations each transfer
    // --beneficial-mutation-rate will not be used. (0=OFF)", 0); options("new-simulator", "Use the new simulator
    // type.", NoArgument);
  }

}  // namespace bpopsim

namespace bpopsim {

  vector<vector<string>> ParseGroupedList(const string& value)
  {
    // Repeated options arrive joined with "\n"; treat that like '|'.
    const string groups_text = substitute(value, "\n", "|");

    vector<vector<string>> groups;
    for (const string& group_text : split(groups_text, "|")) {
      vector<string> items;
      for (const string& raw_item : split(group_text, ",")) {
        const string item = chomp(raw_item.substr(
            raw_item.find_first_not_of(" \t") == string::npos ? raw_item.size() : raw_item.find_first_not_of(" \t")));
        if (item.empty()) {
          throw std::invalid_argument("Empty item in list \"" + value + "\" (check for stray ',' or '|').");
        }
        items.push_back(item);
      }
      if (items.empty()) {
        throw std::invalid_argument("Empty group in list \"" + value + "\" (check for stray '|').");
      }
      groups.push_back(items);
    }
    if (groups.empty()) throw std::invalid_argument("Empty list where at least one value is required.");
    return groups;
  }

  namespace {

    bool SameShape(const vector<vector<string>>& a, const vector<vector<string>>& b)
    {
      if (a.size() != b.size()) return false;
      for (size_t g = 0; g < a.size(); ++g) {
        if (a[g].size() != b[g].size()) return false;
      }
      return true;
    }

  }  // namespace

  SimulationParameters ParseSimulationParameters(OptionParser& options)
  {
    SimulationParameters p;
    p.maximum_number_of_transfers = from_string<uint32_t>(options["number-of-transfers"]);
    p.initial_population_size = from_string<double>(options["initial-population-size"]);
    p.initial_population_size_after_transfer = from_string<double>(options["population-size-after-transfer"]);
    p.generations_per_transfer = from_string<double>(options["generations-per-transfer"]);

    p.initial_fitness = from_string<double>(options["initial-fitness"]);
    p.multiplicative_fitness_effects = options.count("multiplicative-fitness");
    // Mutation categories: '|' groups of ',' items, identical shape for -u, -s
    // and -f (a single -f letter is broadcast to every category).
    const vector<vector<string>> rates = ParseGroupedList(options["mutation-rates"]);
    const vector<vector<string>> effects = ParseGroupedList(options["fitness-effects"]);
    vector<vector<string>> models = ParseGroupedList(options["fitness-effect-model"]);

    if (!SameShape(rates, effects)) {
      throw std::invalid_argument("The -u and -s lists must have the same shape (same '|' groups and ',' items).");
    }
    if (models.size() == 1 && models[0].size() == 1) {
      // Broadcast one model to the shape of -u.
      const string model = models[0][0];
      models.clear();
      for (const auto& group : rates) models.emplace_back(group.size(), model);
    } else if (!SameShape(rates, models)) {
      throw std::invalid_argument(
          "The -f list must be a single model or have the same shape as -u (same '|' groups and ',' items).");
    }

    p.mutation_rates_per_division.clear();
    p.mutation_fitness_effects.clear();
    p.mutation_fitness_effect_models.clear();
    p.mutation_category_group_sizes.clear();
    for (size_t g = 0; g < rates.size(); ++g) {
      p.mutation_category_group_sizes.push_back(rates[g].size());
      for (size_t i = 0; i < rates[g].size(); ++i) {
        p.mutation_rates_per_division.push_back(from_string<double>(rates[g][i]));
        p.mutation_fitness_effects.push_back(from_string<double>(effects[g][i]));
        if (models[g][i].size() != 1) {
          throw std::invalid_argument("Fitness effect model (-f) entries must be single letters: '" + models[g][i] +
                                      "'.");
        }
        p.mutation_fitness_effect_models.push_back(models[g][i][0]);
      }
    }

    p.average_fitness_end_condition = from_string<double>(options["average-fitness-end-condition"]);

    if (options.count("first-mutation-fitness-effects")) {
      p.first_mutation_fitness_effects = from_string<vector<double>>(options["first-mutation-fitness-effects"]);
      p.first_mutation_fitness_effects_string = options["first-mutation-fitness-effects"];
    }

    p.Finalize();
    return p;
  }

  OutputParameters ParseOutputParameters(OptionParser& options)
  {
    OutputParameters p;
    p.output_directory_name = options["output-folder"];
    p.coarse_graining = from_string<uint32_t>(options["coarse-graining"]);
    auto burn_in = from_string<int64_t>(options["burn-in"]);
    if (burn_in < 0) throw std::invalid_argument("--burn-in must be 0 or greater.");
    p.burn_in = static_cast<uint32_t>(burn_in);
    p.diverged_mutation_depth = from_string<uint32_t>(options["diverged-mutation-depth"]);
    p.output_diverged_frequencies = options.count("output-diverged-frequencies");
    p.minimum_output_frequency = from_string<double>(options["minimum-output-frequency"]);
    p.output_longest_sweep_time = options.count("output-longest-sweep-time");
    p.output_marker_frequencies = options.count("output-marker-frequencies");
    p.verbose = options.count("verbose");
    p.debug = options.count("debug");
    return p;
  }

  MarkerDivergence ParseMarkerDivergence(OptionParser& options)
  {
    MarkerDivergence p;
    p.marker_states = from_string<uint32_t>(options["marker-states"]);
    p.divergence_end_fraction = from_string<double>(options["marker-divergence-end-fraction"]);
    return p;
  }

}  // namespace bpopsim
