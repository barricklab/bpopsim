# bpopsim

Simulator of asexual (bacterial) population evolution under serial transfer.
A population grows by binary division for a fixed number of generations per
transfer, is diluted by binomial sampling, and repeats. Beneficial mutations
arise at a per-division rate with fitness effects drawn from a chosen model.
Every genotype that ever arose is a node in an ancestry tree, which lets the
program report clade frequencies, line-of-descent statistics, sweep sizes and
times, and a Muller matrix for plotting population dynamics.

Written ~2010 (Barrick lab, originally alongside breseq). Pre-C++11 style
throughout. A staged modernization roadmap lives in `docs/MODERNIZATION.md`;
work through it one phase at a time.

## Repository layout

- `src/c/bpopsim/` C++ core. Our code: `main.cpp`, `cPopulation.{h,cc}`,
  `cSubpopulation.{h,cc}`, `common.h`.
- Vendored third-party code in the same directory, not ours to restyle:
  `tree.hh` / `tree_util.hh` (Kasper Peeters' STL-like n-ary tree, holds the
  genotype ancestry), `anyoption.{h,cpp}` (command-line parser, locally
  modified), `icsilog.{h,c}` (fast log lookup; compiled in but never called).
- `src/c/bpopsim.1` man page (still the Xcode template, not filled in).
- `CMakeLists.txt` build; `tests/` holds the CTest golden regression test
  (`tests/CMakeLists.txt`, `tests/golden_test.cmake`, `tests/golden/`).
- `binarydist.sh` builds a release tarball with GSL linked statically. On
  macOS it uses the system Apple clang and builds GSL from source as a static
  universal (arm64 + x86_64) archive into `binarydist_build/` (gitignored,
  cached across runs), mirroring breseq's release script.
- `environment.yml` conda environment with every build and run prerequisite.
- Analysis helpers, out of scope for C++ work: `src/R/` (Muller plots,
  average fitness), `src/python/muller_plot.py`, `src/matlab/`,
  `src/mathematica/`, `src/perl/MarkerDivergence/` (marker divergence
  statistics, its own Module::Build project).
- An Xcode project is not checked in; generate one with
  `cmake -S . -B build-xcode -G Xcode` when wanted.

## Build and run

All build and run prerequisites come from a local conda environment defined
in `environment.yml` (compiler, cmake, make, GSL). Create it once at
`./env` (gitignored), then run commands inside it with `conda run` so nothing
needs activating:

```
conda env create -f environment.yml --prefix=$PWD/env

conda run -p $PWD/env cmake -S . -B build
conda run -p $PWD/env cmake --build build --parallel
conda run -p $PWD/env ctest --test-dir build --output-on-failure
build/bpopsim -u 1E-7 -s 0.1 -n 50 -r 3 -o out \
    --output-average-fitness --output-muller
```

CMake picks up the env's compiler shims (`cc`, `c++`) from PATH and finds
the env's GSL through `find_package(GSL)`. The binary carries an rpath into
`env/lib`, so it runs without activating the env. Build type defaults to
Release. `-DBPOPSIM_STATIC_GSL=ON` links GSL statically (used by
`binarydist.sh`).

Unlike breseq, there is no separate runtime environment: GSL is the only
runtime dependency and it is also a build dependency.

At least one `--output-*` flag is required or the program prints usage and
exits. `bpopsim --help` lists all options.

Regression oracle: `ctest` runs `golden_regression`, which executes the
fixed-seed command recorded in `tests/golden/COMMAND` and requires every
output file to be byte-identical to `tests/golden/`. Run it after any change.
Output must stay identical unless the commit says otherwise and regenerates
the golden files. If running the command by hand in zsh, pass the arguments
as an array; an unquoted `$ARGS` string is passed as a single word and the
program falls back to printing usage.

## Architecture

- `main.cpp`: defines options with `AnyOption`, builds one `gsl_rng` seeded
  once for the whole execution, loops over replicates constructing a fresh
  `cPopulation` each time, then writes per-replicate and summary outputs.
- `cPopulation` (`cPopulation.h`): owns `SimulationParameters`,
  `OutputParameters`, `MarkerDivergence`, a `vector<cSubpopulation>` of
  currently living genotypes with their cell counts, and a `tree<cGenotype>`
  of every genotype ever created. `RunSimulation()` loops transfers;
  within a transfer it alternates
  `ProcessCellDivisionTimeStepExactWithFractionalCells()` (grow until the
  next whole-cell event or the next mutation), `MutateExactWithFractionalCells()`
  (pick a parent subpopulation weighted by growth, spawn a child genotype),
  and `TransferResampleDilution()` (binomial dilution, cull extinct
  genotypes). `RecordStatisticsAtTransfer()` fills `replicate_statistics`.
- `cSubpopulation` / `cGenotype` (`cSubpopulation.h`): a subpopulation is a
  cell count, a marker state, and an iterator into the genotype tree.
  `cGenotype::AddOneMutation()` applies the fitness-effect model
  (`u` uniform, `e` exponential, `o` one-time per category).
- `cStatistics` (`cPopulation.h`): a vector of `cReplicateStatistics`;
  the `Output*` methods average across replicates and write the summary
  `.tab` files.
- `common.h`: string/number helpers (`split`, `join`, `from_string<>`,
  `to_string`), `create_path`, and the global `g_verbose`.

## Working rules

- **Do not commit without explicit OK from the user.** Make the changes,
  report what changed, and wait. The user often edits the result or asks for
  changes before anything is committed. Approval of a plan is not approval to
  commit; ask again when the work is ready.

## Conventions and gotchas

- GSL is the single permitted external dependency. It stays, because GSL
  gives identical random streams on every platform for a given seed and
  version. Do not replace it with `<random>`; its distributions are not
  bit-for-bit specified across libc++ and libstdc++.
- One `gsl_rng` (taus2) per execution. Never reseed per replicate.
- Option values in scientific notation must be quoted on the command line
  (`-p "5E6"`).
- `-u` and `-s` may be repeated to define mutation categories; the counts
  must match or the program exits with usage.
- `-1` repeated fixes the fitness effects of the first N mutations.
- Output precision is set with `setprecision` in the output methods; changing
  float formatting changes output files.
- No new third-party dependencies without discussion. Vendored single
  headers (like `tree.hh`, and doctest for tests) are acceptable.
- Seeded output should stay byte-identical across refactors. If a change is
  intentionally behavioral, say so in the commit message.
- Run `ctest` before asking to commit.

## Known issues (fix via roadmap, not ad hoc)

- `--output-largest-sweep-size` is `assert(false)`, not implemented.
- `GenotypesAboveThreshold()` reads `clade_frequencies` instead of
  `genotype_frequencies` (copy-paste), so the genotype frequency file
  selects columns by clade frequency. Fixing this changes output.
- `TransferResampleDilution()` declares `double by_color[0]` and writes to
  it when markers are `'r'`/`'w'` (marker-divergence mode only).
- Seed is read with `from_string<int16_t>` into a `uint16_t`; the default
  `time(NULL) * getpid()` is truncated to 16 bits.
- `chomp()` in `common.h` drops one character too many.
- `CalculateSimilarity()` in `cPopulation.cc` is a non-void function with no
  return statement (unreferenced dead code, `-Wreturn-type`).
- `main()` declares `char* argv_new[argc]`, a variable-length array, which is
  a Clang extension in C++ (`-Wvla-cxx-extension`).
- The build currently emits about 60 warnings under `-Wall -Wextra`, mostly
  `-Wignored-qualifiers` (`const` on by-value returns), `-Wsign-compare`,
  and `-Wunused-parameter`. Phase 5 of the roadmap clears them.
- `OutputParameters::burn_in` is `int32_t` but parsed as `uint32_t`.
- Roughly a dozen `cPopulation` methods below the comment
  "Need to check output below here" are unreferenced, as is the
  `CONVERT-TREE` subcommand in `main.cpp`.
- `SimulationParameters`, `OutputParameters`, `cStatistics` take `AnyOption&`
  in their constructors, so nothing is constructible in a test without a
  parsed command line.

## Uncommitted work

As of 2026-11-07 there is a stash (or working-tree diff) that changes
`mutation_rates_per_division` and `mutation_fitness_effects` to
`vector<vector<double>>` with a `|` separator between groups and adds a
`mutation_path` field to `cGenotype`. It is incomplete and does not compile.
Do not build on it until it is finished or dropped.
