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
- `src/c/bpopsim.xcodeproj/` hand-maintained Xcode project; superseded once
  CMake lands (`cmake -G Xcode` regenerates one).
- Analysis helpers, out of scope for C++ work: `src/R/` (Muller plots,
  average fitness), `src/python/muller_plot.py`, `src/matlab/`,
  `src/mathematica/`, `src/perl/MarkerDivergence/` (marker divergence
  statistics, its own Module::Build project).
- `bootstrap.sh`, `configure.ac`, `Makefile.am`, `aux_build/m4/` autotools
  build. `binarydist.sh` builds a static universal macOS binary tarball.

## Build and run

All build and run prerequisites come from a local conda environment defined
in `environment.yml` (compiler, make, cmake, autotools, GSL). Create it
once at `./env` (gitignored), then run commands inside it with `conda run`
so nothing needs activating:

```
conda env create -f environment.yml --prefix=$PWD/env

conda run -p $PWD/env ./bootstrap.sh
conda run -p $PWD/env sh -c './configure --with-gsl-prefix="$CONDA_PREFIX" \
    CC="$CONDA_PREFIX/bin/clang" CXX="$CONDA_PREFIX/bin/clang++"'
conda run -p $PWD/env make -j
src/c/bpopsim/bpopsim -u 1E-7 -s 0.1 -n 50 -r 3 -o out \
    --output-average-fitness --output-muller
```

Why the explicit flags: conda-forge's compiler packages no longer export
`CC`/`CXX` on activation, so without them autoconf picks up Apple's
`/usr/bin/gcc` and the GSL link check fails. `--with-gsl-prefix` is needed
because this machine also has a stale x86_64 GSL 1.13 in `/usr/local` that
`gsl.m4` finds first. The `sh -c` wrapper lets `$CONDA_PREFIX` expand inside
`conda run`. The resulting binary carries an rpath into `env/lib`, so it
runs without activating the env.

Unlike breseq, there is no separate runtime environment: GSL is the only
runtime dependency and it is also a build dependency.

If `./configure` or `Makefile` is missing, the bootstrap and configure steps
have not been run in this checkout. After switching between toolchains, run
`make distclean` first; otherwise `make` re-runs the old configure with its
old arguments.

At least one `--output-*` flag is required or the program prints usage and
exits. `bpopsim --help` lists all options. There is no test target yet
(`make test` is an empty rule in `Makefile.am`).

Regression oracle: `tests/golden/` holds the output of the fixed-seed run
recorded in `tests/golden/COMMAND`. Rerun that command into a scratch
directory and `diff -r` against `tests/golden/` after any change. Output must
stay byte-identical unless the commit says otherwise and regenerates it.
In zsh, pass the arguments as an array; an unquoted `$ARGS` string is passed
as a single word and the program falls back to printing usage.

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
- Once tests exist (roadmap Phase 2), run `ctest` before committing.

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
