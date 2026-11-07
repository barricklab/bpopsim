# bpopsim

Simulator of asexual (bacterial) population evolution under serial transfer.

A population grows by binary cell division for a fixed number of generations,
is diluted by binomial sampling to start the next transfer, and repeats.
Beneficial mutations arise at a per-division rate with fitness effects drawn
from a chosen model. Every genotype that ever arose is a node in an ancestry
tree, which lets the program report clade and genotype frequencies over time,
line-of-descent statistics such as the longest sweep, and a Muller matrix for
plotting population dynamics.

bpopsim was written in the Barrick lab to model experiments like the *E. coli*
long-term evolution experiment, where populations are propagated by daily
1:100 dilution (6.64 generations per transfer).

## Building

All build and run prerequisites come from a conda environment described in
`environment.yml`: a C++17 compiler, CMake, and the
[GNU Scientific Library](https://www.gnu.org/software/gsl/) (GSL), which
supplies the random number generator and distribution samplers.

```sh
conda env create -f environment.yml --prefix=$PWD/env

conda run -p $PWD/env cmake -S . -B build
conda run -p $PWD/env cmake --build build --parallel
conda run -p $PWD/env ctest --test-dir build --output-on-failure
```

The binary is `build/bpopsim`. It carries an rpath into the environment's
`lib` directory, so it runs without activating the environment.

`./binarydist.sh` builds a redistributable tarball with GSL linked statically
(a universal arm64 + x86_64 binary on macOS).

Without conda, any C++17 compiler, CMake 3.15+, and a GSL development package
(`libgsl-dev` on Debian/Ubuntu, `gsl` from Homebrew) will do.

## Running

At least one `--output-*` option is required. For example, three replicate
populations of 5 x 10^6 cells, a beneficial mutation rate of 10^-7 per cell
division, selection coefficient 0.1, 500 transfers:

```sh
build/bpopsim -r 3 -p 5E6 -N 5E6 -T 6.64 -u 1E-7 -s 0.1 -n 500 \
    -o output --output-average-fitness --output-muller
```

`bpopsim --help` lists every option. The main ones:

| Option | Meaning |
|---|---|
| `-r, --replicates` | Number of replicate populations |
| `-d, --rng-seed` | Random seed. 0 (default) picks one and prints it, so any run can be reproduced |
| `-p, --initial-population-size` | Cells at the start of the simulation |
| `-N, --population-size-after-transfer` | Cells after each dilution |
| `-T, --generations-per-transfer` | Doublings per transfer; the dilution factor is 2^T |
| `-n, --number-of-transfers` | Transfers to simulate |
| `-u, --mutation-rates` | Beneficial mutation rate per cell division. Repeat to define categories |
| `-s, --fitness-effects` | Fitness effect per mutation. Repeat to match `-u` categories |
| `-f, --fitness-effect-model` | `u` fixed effect, `e` exponential with mean `-s`, `o` one mutation per category |
| `-m, --multiplicative-fitness` | Fitness = f (1 + s) instead of f + s |
| `-1, --first-mutation-fitness-effects` | Override the effects of the 1st, 2nd, ... mutations |
| `-b, --burn-in` | Transfers to run before recording begins (in addition to `-n`) |
| `-c, --coarse-graining` | Record every this many transfers |
| `-o, --output-folder` | Where output files go |

Quote values in scientific notation (`-p "5E6"`) if your shell would otherwise
mangle them.

## Output files

All files are tab-separated text written to the output folder.

| Option | File(s) | Contents |
|---|---|---|
| `--output-average-fitness` | `average_fitness.tab` | Mean fitness per transfer, one row per replicate |
| `--output-average-mutations` | `average_mutations.tab`, `average_mutations_category_N.tab` | Mean mutations per cell per transfer, total and per category |
| `--output-genotype-frequencies` | `genotype_frequencies_<min>_<rep>.dat` | Frequency of each genotype that ever exceeded `--minimum-output-frequency` |
| `--output-clade-frequencies` | `clade_frequencies_<min>_<rep>.dat` | Same for clades (a genotype plus all its descendants); header rows flag the line of descent and each mutation's fitness effect |
| `--output-diverged-frequencies` | `mutation_divergence_depth_frequency_<d>.tab` | Largest clade frequency diverged by at least `d` mutations from the final dominant line of descent, up to `--diverged-mutation-depth` |
| `--output-longest-sweep-time` | `longest_sweep_time.tab` | Longest time (transfers) for a mutation on the line of descent to go from at most 5% to at least 95% |
| `--output-muller` | `muller_matrix_<rep>.dat` | Muller matrix, `--muller-resolution` rows per transfer, cell values are genotype ids (`--muller-lod` makes the line of descent negative) |

Scripts for plotting these in R, Python, and MATLAB live under `src/R`,
`src/python`, and `src/matlab`.

## Reproducibility

Given the same seed, GSL version, and command line, bpopsim produces the same
output on every platform at the precision it writes. The test suite pins the
random stream for one seed and compares a full run against stored reference
output (`tests/golden/`), and continuous integration runs that check on both
Linux and macOS. Keep the GSL version pinned in `environment.yml` if you need
to reproduce old runs exactly.

One caveat: GSL's uniform, integer, and binomial draws are integer algorithms
and are bit-identical everywhere, but the exponential sampler calls the C
library's `log()`, whose last binary digit can differ between glibc and
Apple's libm. Written output is identical on Linux and macOS at the precision
bpopsim uses, but bit-for-bit identity of very long runs across operating
systems is not guaranteed.

## Development

See `CLAUDE.md` for the code layout, conventions, and known issues. Tests use
[doctest](https://github.com/doctest/doctest) (vendored) under CTest; run
them with `ctest --test-dir build`.

## License

GPL v3. See `COPYING`.
