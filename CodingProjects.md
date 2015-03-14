# Software Engineering #

  * ~~Remove Boost and change to any\_option.cpp/any\_option.h from breseq.~~
  * Add subcommands.
  * Audit memory usage.

# Muller Plot #

  * ~~Read in new input files and produce Muller plot.~~
  * Example format:

```

anc 1.00 1.00 1.00 1.00 1.00
anc,pykF 0.00 0.10 0.40 0.85 1.00
anc,pykF,nadR 0.00 0.00 0.10 0.40 0.90
anc,pykF,rbs1 0.00 0.00 0.20 0.40 0.10```

  * Lines are space-delimited.
  * Genotypes are comma-delimited.
  * Assume parents always come first and they are only one generation away.
  * Each line is a unique genotype.
  * The parent's genotype frequency must always be greater than or equal to the sum of its children.

# Mutation Accumulation #

  * Throw down user-specified number of mutations in a transfer
  * Record frequency of mutation in final population.

# Veracode Genotyping #

Questions about what is possible:
  * ~~Log plot of mutations that fix over time.~~
  * ~~How low a frequency would you have to detect to tell the true order of mutations?~~
  * Graph of "derivative" of mutation trajectories. Take current time point and previous time point to calculate slope. Graph slope over time to compare.
    * What is the apparent selection coefficient for a mutation?
    * Selection coefficient = log( (number at time T with mutation)**(dilution factor) / (number at time T-1 with mutation)) / log ((number at time T without mutation)**(dilution factor) / (number at time T-1 without mutation))
    * What if we knew linkage between mutations? Could we improve things?

Making the model more realistic:
  * Add neutral mutations in at a fixed frequency relative to beneficial mutations (10%? 50%?).
  * Draw mutations from an exponential fitness distribution (after finding what lambda = mean parameter is characteristic of the population).