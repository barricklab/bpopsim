#include "options.h"

namespace bpopsim {

  void DefineBpopsimOptions(AnyOption& options)
  {
    // @JEB: doubles must be put in quotes if in scientific notation
    options("help,h", "Produce this help message", TAKES_NO_ARGUMENT);
    options("verbose,v", "Output extra details of runs", TAKES_NO_ARGUMENT);
    options("debug", "Perform extra checks on simulation", TAKES_NO_ARGUMENT);

    options.addUsage("");
    options.addUsage("==== Simulation Options ====");
    options.addUsage("");
    options("rng-seed,d", "Seed for random number generator. [0 = choose a random seed and report it]", 0);
    options("replicates,r", "Replicates", 1);
    options("initial-population-size,p", "Initial population size.", "5E6");
    options("generations-per-transfer,T", "Generations per transfer", 6.64);
    options("population-size-after-transfer,N", "Population size after transfer", "5E6");
    options("number-of-transfers,n", "Max number of transfer to replicate", 50);
    options("initial-fitness,z", "Initial fitness.", 1.0);
    options("multiplicative-fitness,m", "Use multiplicative fitness model [f * (1 + s)]. Default = additive [f + s].",
            TAKES_NO_ARGUMENT);
    options("mutation-rates,u",
            "Mutation rate per cell division. Supply option multiple times to define categories of mutations.", "1E-7");
    options("fitness-effects,s",
            "Mean fitness increment per mutation. Supply option multiple times to define categories of mutations.",
            0.1);
    options("fitness-effect-model,f",
            "Distribution of mutation fitness effects. Choices are 'u' for uniform, 'e' for exponential distributions, "
            "'o' for only one beneficial mutation possible in that category.",
            'u');
    options("first-mutation-fitness-effects,1",
            "Use these fitness effects for the first mutations. Supply option multiple times to define the fitness "
            "effects of the 1st, 2nd, 3rd, etc. mutations. These values will override the -s option until the number "
            "of mutations exceeds the number of these values.");
    options("average-fitness-end-condition",
            "End a replicate when the average fitness of the population is this or greater (0=OFF)", 0);

    options.addUsage("");
    options.addUsage("==== Marker Divergence Options ====");
    options.addUsage("");
    options("marker-states,k",
            "Begin with the initial population divided equally into this many subpopulations each with a different "
            "state of a neutral genetic marker. Do not track genotypes. Only track the neutral marker state of each "
            "subpopulation. Default (0=OFF).",
            0);
    options("max-marker-divergence-ratio,x",
            "Stop if divergence factor between the subpopulations with the most abundant and next-most abundant marker "
            "state exceeds this factor.",
            100);
    options("transfer-interval-to-print,t", "Marker divergence printing intervals.", 1);
    options.addUsage("");
    options.addUsage("==== Statistical Output Files and Options ====");
    options.addUsage("");
    options("output-folder,o", "Base output folder. All output files will be created here.", ".");
    options("burn-in,b",
            "Perform this number of transfers before beginning to output statistics (before transfer numbered 0).", 0);
    options("coarse-graining,c", "Only print stats every this many transfers to all output files.", 1);
    options("output-average-fitness", "Output average fitness to file 'average_fitness.tab'", TAKES_NO_ARGUMENT);
    options("output-average-mutations", "Output mutation numbers to file 'average_mutations_count.tab'",
            TAKES_NO_ARGUMENT);
    options("output-clade-frequencies",
            "Output clade (mutation) frequencies (above minimum) for replicate X to file 'clade_frequencies_X.tab'",
            TAKES_NO_ARGUMENT);
    options("output-genotype-frequencies",
            "Output genotype frequencies (above minimum) for replicate X to file 'genotype_frequencies_X.tab'.",
            TAKES_NO_ARGUMENT);
    options("minimum-output-frequency",
            "Only output clade/genotype frequencies that are above this frequency at some time point.", 0.01);
    options("output-diverged-frequencies",
            "Output frequencies of clades that differ by at least X mutations from the line of descent to the final "
            "dominant to 'diverged_frequencies_X.tab'",
            TAKES_NO_ARGUMENT);
    options("diverged-mutation-depth", "Maximum depth of diverged mutations to output.", 10);
    options("output-largest-sweep-size",
            "Output the largest number of mutations sweeping simultaneously along the line of descent to the final "
            "dominant to 'largest_sweep_size.tab'",
            TAKES_NO_ARGUMENT);
    options("minimum-sweep-frequency-difference",
            "Adjacent mutations along the line of descent must have frequncies that are different by at least this "
            "threshold at one recorded time point to be considered part of different sweeps.",
            0.05);
    options("output-longest-sweep-time",
            "Output the longest time that a mutation along the line of descent to the final dominant took to sweep "
            "from ≥5% to ≥95% frequency in the population to 'longest_sweep_time.tab'",
            TAKES_NO_ARGUMENT);
    options.addUsage("");
    options.addUsage("==== Muller Plot Output Files and Options ====");
    options.addUsage("");
    options("output-muller",
            "Output Muller matrix for replicate X to file 'muller_matrix_X.dat'. Cell values unique genotypes. This "
            "file can be used to make plots of population dynamics.",
            TAKES_NO_ARGUMENT);
    options("muller-resolution", "Muller matrix vertical resolution.", 200);
    options("muller-lod",
            "In the Muller matrix, give genotypes on the line of descent to the final dominant negative numbers, so "
            "that they can be plotted in separate colors.",
            TAKES_NO_ARGUMENT);

    // Not (re)implemented
    //  options("output-dominant-genotypes", "Output dominant genotype frequencies to file 'dominant_genotypes.csv'",
    //  TAKES_NO_ARGUMENT);
    //  ("output-time-to-sweep", "Print time to sweep for each mutation", TAKES_NO_ARGUMENT)
    //  ("output-max-diff", "Print max difference of sweeping mutations", TAKES_NO_ARGUMENT)
    //  ("output-single-fitness", "Print the fitness of the single cell.", TAKES_NO_ARGUMENT)
    //  ("output-sweeping-descent-fitness", "Output the average fitness of sweeping descent per time.",
    //  TAKES_NO_ARGUMENT)

    // options.addUsage("");
    // options.addUsage("==== Other Options ====");
    // options.addUsage("");
    // options("exact-mutations-per-transfer,y", "Cause exactly this number of mutations each transfer
    // --beneficial-mutation-rate will not be used. (0=OFF)", 0); options("new-simulator", "Use the new simulator
    // type.", TAKES_NO_ARGUMENT);
  }

}  // namespace bpopsim

namespace bpopsim {

  SimulationParameters ParseSimulationParameters(AnyOption& options)
  {
    SimulationParameters p;
    p.maximum_number_of_transfers = from_string<uint32_t>(options["number-of-transfers"]);
    p.initial_population_size = from_string<double>(options["initial-population-size"]);
    p.initial_population_size_after_transfer = from_string<double>(options["population-size-after-transfer"]);
    p.generations_per_transfer = from_string<double>(options["generations-per-transfer"]);

    p.initial_fitness = from_string<double>(options["initial-fitness"]);
    p.multiplicative_fitness_effects = options.count("multiplicative-fitness");
    p.mutation_rates_per_division = from_string<vector<double>>(options["mutation-rates"]);
    p.mutation_fitness_effects = from_string<vector<double>>(options["fitness-effects"]);
    p.mutation_fitness_effect_model = options["fitness-effect-model"];

    p.average_fitness_end_condition = from_string<double>(options["average-fitness-end-condition"]);

    if (options.count("first-mutation-fitness-effects")) {
      p.first_mutation_fitness_effects = from_string<vector<double>>(options["first-mutation-fitness-effects"]);
      p.first_mutation_fitness_effects_string = options["first-mutation-fitness-effects"];
    }

    p.Finalize();
    return p;
  }

  OutputParameters ParseOutputParameters(AnyOption& options)
  {
    OutputParameters p;
    p.output_directory_name = options["output-folder"];
    p.coarse_graining = from_string<uint32_t>(options["coarse-graining"]);
    p.burn_in = from_string<uint32_t>(options["burn-in"]);
    p.diverged_mutation_depth = from_string<uint32_t>(options["diverged-mutation-depth"]);
    p.output_diverged_frequencies = options.count("output-diverged-frequencies");
    p.minimum_output_frequency = from_string<double>(options["minimum-output-frequency"]);
    p.output_largest_sweep_size = options.count("output-largest-sweep-size");
    p.output_longest_sweep_time = options.count("output-longest-sweep-time");
    p.verbose = options.count("verbose");
    p.debug = options.count("debug");
    return p;
  }

  MarkerDivergence ParseMarkerDivergence(AnyOption& options)
  {
    MarkerDivergence p;
    p.marker_states = from_string<uint32_t>(options["marker-states"]);
    p.max_divergence_factor = from_string<double>(options["max-marker-divergence-ratio"]);
    p.transfer_interval_to_print = from_string<double>(options["transfer-interval-to-print"]);
    return p;
  }

}  // namespace bpopsim
