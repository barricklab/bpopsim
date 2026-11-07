# bpopsim

Simulator of asexual (bacterial) population evolution under serial transfer.
A population grows by binary division for a fixed number of generations per
transfer, is diluted by binomial sampling, and repeats. Beneficial mutations
arise at a per-division rate with fitness effects drawn from a chosen model.
Every genotype that ever arose is a node in an ancestry tree, which lets the
program report clade frequencies, line-of-descent statistics, sweep sizes and
times, and a Muller matrix for plotting population dynamics.

C++17, built with CMake inside a conda environment, tested with doctest and
a golden-output regression under CTest. `docs/MODERNIZATION.md` is the
development log; nothing in it is needed to work on the code.

## Repository layout

- `src/c/bpopsim/` C++ core, all ours: `main.cpp`, `cPopulation.{h,cc}`,
  `cSubpopulation.{h,cc}`, `option_parser.{h,cc}` (small GNU/POSIX style
  command-line parser), `options.{h,cc}` (the bpopsim option table and the
  `Parse*` functions, shared by `main` and the tests), `random.{h,cc}` (the
  GSL wrapper), `common.h` (string helpers).
- `src/c/third_party/` vendored code, built as a `SYSTEM` include directory
  and excluded from clang-format/clang-tidy: `tree.hh` / `tree_util.hh`
  (Kasper Peeters' STL-like n-ary tree, holds the genotype ancestry;
  `tree.hh` is locally patched to drop its `using namespace std;` and
  qualify its std uses).
- `src/c/bpopsim.1` man page, covering every option and output file. Update
  it together with `options.cc` and the option table in `README.md`.
- `.github/workflows/build-and-test.yml` builds and runs `ctest` plus a
  clang-format check on Ubuntu and macOS from `environment.yml`.
- `CMakeLists.txt` build. `tests/` holds the doctest unit tests
  (`test_*.cc`, one `TEST_CASE` per CTest entry), the vendored doctest in
  `tests/third_party/doctest/`, and the golden regression test
  (`golden_test.cmake`, `golden/`).
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

Tests: `ctest` runs every doctest `TEST_CASE` in `tests/test_*.cc` as its own
entry, then `golden_regression`. To run only the unit tests with doctest's own
reporter: `build/tests/bpopsim_tests` (add `-tc="<name>*"` to filter). New
tests go in an existing `test_*.cc` or a new one listed in
`tests/CMakeLists.txt`.

Regression oracle: `golden_regression` executes the fixed-seed command
recorded in `tests/golden/COMMAND` and requires every output file to be
byte-identical to `tests/golden/`. Run it after any change.
Output must stay identical unless the commit says otherwise and regenerates
the golden files. If running the command by hand in zsh, pass the arguments
as an array; an unquoted `$ARGS` string is passed as a single word and the
program falls back to printing usage.

## Architecture

- `main.cpp`: parses the command line with `OptionParser` (option table in
  `options.cc`), converts it into `SimulationParameters`, `OutputParameters`
  and `MarkerDivergence` via the `Parse*` functions in `options.cc`, calls
  `Validate()`, builds one `RandomGenerator` seeded once for the whole
  execution, loops over replicates constructing a fresh `cPopulation` from
  the structs each time, then writes per-replicate and summary outputs.
  Parse errors and validation errors are `std::invalid_argument`; `main`
  prints the message with the usage text and returns -1. `--help` returns 0.
- `OptionParser` (`option_parser.h`): accepts `--long value`, `--long=value`,
  `-s value`, `-svalue`, and combined flags `-abc`. Repeated options join
  their values with `"\n"`, which `from_string<vector<T>>` reads; that is how
  `-u`/`-s` categories reach `SimulationParameters`. Flags report `"true"`.
  Unknown options, missing values, and positional arguments are errors.
- Parameter structs (`cPopulation.h`): plain data with defaults equal to the
  command-line defaults. `SimulationParameters::Finalize()` computes derived
  values (dilution factor, binomial probability, category fractions); call it
  after changing inputs. Only `options.cc` touches `OptionParser`, so tests
  construct the structs directly.
- `cPopulation` (`cPopulation.h`): owns copies of the three parameter
  structs, a `vector<cSubpopulation>` of currently living genotypes with
  their cell counts, and a `tree<cGenotype>` of every genotype ever created.
  `RunSimulation()` loops transfers;
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
  `to_string`) and `create_path`. Verbose and debug switches live on
  `OutputParameters`, not in globals.

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
- All randomness goes through `RandomGenerator` (`random.h`), which wraps a
  taus2 `gsl_rng`. `random.cc` is the only file that includes GSL headers;
  keep it that way. One generator per execution, seeded once. Never reseed
  per replicate. `tests/test_random.cc` pins the stream for seed 12345.
- Cross-platform limits of that reproducibility: `Uniform`, `UniformInt` and
  `Binomial` are integer algorithms and bit-identical on Linux and macOS.
  `Exponential` goes through libm `log()`, whose last bit can differ between
  glibc and Apple libm. Written output (6 significant digits) is identical on
  both; do not write tests that expect bit-exact exponential draws across
  platforms. Use `Approx` with a tight epsilon as `test_random.cc` does.
- `assert()` is active in every build type (CMake strips `-DNDEBUG`). The
  simulation uses asserts as internal consistency checks and the tests rely
  on them.
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
- Run `ctest` before asking to commit. The build is warning-free under
  `-Wall -Wextra` and must stay that way.
- Style is enforced by `.clang-format` (run
  `conda run -p $PWD/env clang-format -i <files>` on our sources only; vendored
  files are listed in `.clang-format-ignore`). `.clang-tidy` holds the
  project's check set; run it per file with
  `conda run -p $PWD/env clang-tidy -p build src/c/bpopsim/<file>`.
- No header has a `using namespace` directive. `common.h` helpers are at
  global scope and fully `std::`-qualified. Inside `namespace bpopsim`,
  `cSubpopulation.h` carries using-declarations for `map`, `pair`, `set`,
  `string` and `vector`, which is what lets the other bpopsim headers use
  those names unqualified; anything else from `std` is written qualified.
  `.cc` files and tests may `using namespace std;` at file scope. The
  vendored `tree.hh` is patched to qualify its own `std` uses.
- `OptionParser::PrintUsage` probes the terminal width with `ioctl` only when
  stdout is a terminal and otherwise uses 80 columns. Keep the `winsize`
  struct zero-initialized and the return value checked.

## Known issues

- `OutputParameters::burn_in` is `int32_t` but parsed as `uint32_t`.
- Marker-divergence mode (`-k 2`) seeds two genotypes carrying markers
  `'r'` and `'w'` and then runs an ordinary simulation. Nothing
  marker-specific is reported; a divergence stopping rule and a marker-ratio
  output would have to be written.
- Burn-in transfers (`-b`) are performed in addition to `-n`, not within it.

## Stashed work

A git stash named `WIP vector<vector> mutation rates` holds an incomplete
change that makes `mutation_rates_per_division` and `mutation_fitness_effects`
`vector<vector<double>>` with a `|` separator between groups and adds a
`mutation_path` field to `cGenotype`. It does not compile. Do not build on it
until it is finished or dropped.
