# bpopsim modernization roadmap

Staged plan to bring the C++ core to modern C++17 with a test suite while
adding no dependency beyond GSL. Work one phase at a time. Every phase ends
with a green build and, from Phase 2 on, green `ctest` and a byte-identical
regression output for a fixed seed.

## Decisions

| Topic | Decision | Why |
|---|---|---|
| GSL | Keep, behind one small `RandomGenerator` class | Identical seeded streams on every platform. `<random>` distributions differ between libc++ and libstdc++. |
| Build | CMake replaces autotools | Standard; `find_package(GSL)` is built in; `cmake -G Xcode` replaces the hand-maintained project. |
| Tests | Vendored doctest single header + CTest | One file in-repo, fast compile, `doctest::Approx` for numerics. |
| C++ standard | C++17 | Apple clang 17 supports it fully. `std::filesystem`, `std::optional`, structured bindings. |
| WIP edits | Stash before starting | The `vector<vector<double>>` mutation-rate edits do not compile. |

Toolchain at time of writing: Apple clang 17, cmake 3.11 (MacPorts), GSL 1.13.
`cmake_minimum_required(VERSION 3.11)` works today. `brew install cmake` for a
current release is recommended but not required.

## Phase 0: Safety net (done 2026-11-07)

Everything later is measured against this.

- [x] `git stash push -m "WIP vector<vector> mutation rates"` of the four
      modified files (`stash@{0}`).
- [x] Build HEAD. Two fixes were required, both committed separately:
  - `common.h` `file_exists()` returned an `ifstream` as `bool`, rejected
    since C++11 (`explicit operator bool`). Now an explicit cast.
  - `cPopulation.cc` `CullSubpopulationsThatDidNotEstablish()` erased a
    tree node and then dereferenced it to erase its map entry. With a
    modern allocator this left a dangling iterator in
    `genotype_id_map_into_tree`, crashing `--output-clade-frequencies`
    (segfault) and `--output-diverged-frequencies` (assert). Lines swapped.
  - Environment: the only GSL present was x86_64 1.13 from 2009; installed
    GSL 2.8 via Homebrew and configure with
    `--with-gsl-prefix=/opt/homebrew/opt/gsl`.
- [x] Golden output committed under `tests/golden/`, command in
      `tests/golden/COMMAND`. Runs in well under a second, two replicates,
      30 transfers, every working output option. Verified byte-identical
      across two runs. `--output-largest-sweep-size` is excluded because it
      is `assert(false)`.

## Phase 1: CMake build, remove autotools

- [ ] Top-level `CMakeLists.txt`:
  - `project(bpopsim VERSION 0.1 LANGUAGES CXX)`, `CMAKE_CXX_STANDARD 17`,
    extensions off.
  - `find_package(GSL REQUIRED)`.
  - `add_library(bpopsim_core STATIC ...)` from `cPopulation.cc`,
    `cSubpopulation.cc`, `anyoption.cpp`; link `GSL::gsl GSL::gslcblas`.
  - `add_executable(bpopsim src/c/bpopsim/main.cpp)` linking `bpopsim_core`.
  - `-Wall -Wextra -Wpedantic` on our targets. Mark the vendored include
    directory `SYSTEM` or add targeted `-Wno-*` so `tree.hh` and `anyoption`
    do not drown real warnings.
  - `install(TARGETS bpopsim)` and `install(FILES src/c/bpopsim.1 ...)`.
  - `enable_testing()`, `add_subdirectory(tests)` (empty until Phase 2).
- [ ] Delete: `configure.ac`, `Makefile.am`, `src/c/bpopsim/Makefile.am`,
      `bootstrap.sh`, `aux_build/`, `src/c/bpopsim/compile`,
      `src/c/bpopsim.xcodeproj/`.
- [ ] Rewrite `binarydist.sh` around `cmake --build` + `cmake --install`
      (or CPack), keeping the static universal macOS flags.
- [ ] Trim autotools entries from `.gitignore`; add `build*/`.
- [ ] Update the Build section of `CLAUDE.md` and `README`.
- [ ] Verify: `cmake -B build && cmake --build build`, rerun the Phase 0
      command, diff against `tests/golden/`.

## Phase 2: Test harness

- [ ] Vendor `tests/doctest.h` (pin the version in a header comment) and
      doctest's `doctest.cmake` for `doctest_discover_tests`.
- [ ] `tests/CMakeLists.txt`: one `bpopsim_tests` executable from
      `tests/main.cpp` (`#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN`) plus
      `test_*.cc` files, linked to `bpopsim_core`.
- [ ] `tests/test_common.cc`: `split`, `join`, `from_string<double>`,
      `from_string<vector<double>>` (newline-separated), `to_string(double,
      precision)`, `substitute`, `chomp` (write the test for the correct
      behavior; it fails until the off-by-one is fixed).
- [ ] `tests/test_parameters.cc`: `SimulationParameters` derived values:
      `transfer_dilution_factor == 2^T`, `binomial_sampling_transfer_probability
      == 1 / that`, `final_population_size_at_transfer`, category fractions
      sum to 1, mismatch of `-u`/`-s` counts is rejected. Until Phase 3 this
      needs a helper that builds an `AnyOption` from a `vector<string>` argv.
- [ ] Regression test: `add_test` that runs `bpopsim` with the `COMMAND`
      arguments into `${CMAKE_CURRENT_BINARY_DIR}/golden_out` and compares
      each file with `cmake -E compare_files` against `tests/golden/`
      (every file except `COMMAND`).
- [ ] Add a second golden case in marker-divergence mode (`-k 2`) once the
      `by_color[0]` zero-length array in `TransferResampleDilution` is fixed.
- [ ] `ctest --test-dir build --output-on-failure` is green.

## Phase 3: Decouple parameters from option parsing

Biggest structural win. Makes the simulator constructible in tests.

- [ ] Give `SimulationParameters`, `OutputParameters`, `MarkerDivergence`,
      `cStatistics` default member initializers and no `AnyOption&` ctor.
- [ ] Move derived-value computation into `SimulationParameters::finalize()`
      (dilution factor, binomial p, total rate, category fractions).
- [ ] Add `options.h/.cc` with `parse_simulation_parameters(const AnyOption&)`,
      `parse_output_parameters(...)`, `parse_marker_divergence(...)`. Move the
      `-u`/`-s` count check into `validate()` throwing `std::invalid_argument`;
      `main` catches, prints the message and usage, returns nonzero.
- [ ] Replace global `g_verbose` with `OutputParameters::verbose` passed to
      what needs it. Keep output identical.
- [ ] `tests/test_population.cc`: construct `cPopulation` from parameter
      structs and a seeded rng; run N transfers; assert population size stays
      within `[N_after_transfer, N_at_transfer]`, genotype counts are
      consistent, and `CalculatePopulationSize()` equals the tracked
      `current_population_size` (what the `--debug` asserts check today).
- [ ] Golden test still byte-identical.

## Phase 4: Isolate GSL behind `RandomGenerator`

- [ ] `src/c/bpopsim/random.h/.cc`: class owning `gsl_rng*` via
      `std::unique_ptr<gsl_rng, decltype(&gsl_rng_free)>`. Methods:
      `uniform()`, `uniform_int(unsigned long n)`, `binomial(double p,
      unsigned n)`, `exponential(double mean)`, `seed(unsigned long)`.
      Constructed with `gsl_rng_taus2` so streams do not change.
- [ ] Thread `RandomGenerator&` through `cPopulation`, `cSubpopulation`,
      `cGenotype::AddOneMutation`, `SimulationParameters::DetermineMutationCategory`
      in place of `gsl_rng*`. `random.cc` becomes the only file including
      GSL headers. Drop `gsl_matrix.h`.
- [ ] Fix the seed: store as `unsigned long` (what `gsl_rng_set` takes),
      parse with the matching `from_string<>`, default from
      `std::random_device` or `time ^ pid` without 16-bit truncation.
      Fixed-seed output is unchanged; the unseeded default changes, say so in
      the commit.
- [ ] `tests/test_random.cc`: a seeded generator reproduces a short recorded
      sequence of `uniform()`, `uniform_int()`, `binomial()`, `exponential()`
      values (this pins cross-platform reproducibility); binomial and
      exponential sample mean/variance over 1e5 draws match within `Approx`.

## Phase 5: Mechanical modern C++ cleanup

Mostly `clang-tidy` driven. One check per commit so each is reviewable.

- [ ] Add `.clang-format` matching existing style (2-space indent) and
      `.clang-tidy` enabling `modernize-*`, `readability-*`, `performance-*`,
      `bugprone-*`. Exclude vendored files.
- [ ] Range-for over iterator loops; `auto` where the type is obvious;
      `nullptr`; `override`; `using` aliases; `enum class` for `e_colors`;
      `constexpr` for `log_2`; default member initializers; `#pragma once`;
      `[[nodiscard]]` on getters; remove `const` on by-value return types.
- [ ] Delete `cSubpopulation`'s hand-written copy constructor (the defaulted
      one is correct).
- [ ] `std::filesystem::create_directories` replaces `create_path` /
      `mkdir` / `umask`.
- [ ] Remove `using namespace std` from headers; qualify or add
      `using std::vector;` etc. inside `namespace bpopsim`.
- [ ] Replace `make_map` / `make_list` helpers with brace initialization.
- [ ] Keep `setprecision` behavior exactly; the golden test is the guard.

## Phase 6: Dead code and vendored code

- [ ] Delete `icsilog.c/.h`, `fix_flags`, the `CONVERT-TREE` subcommand, and
      the unreferenced `cPopulation` methods: `PrintSingleFitness`,
      `ConvertExternalData`, `PrintWinningFitness`, `PrintExpectationValue`,
      `PrintOut_RedWhiteOnly`, `PrintFrequenciesToScreen*`,
      `CalculateSimilarity`, `TimeToSweep`, `PrintFreqsQuick`,
      `PrintUniqueGenotypes`, and the duplicate `Last_Sweep` /
      `GenotypesFromAncestorToFinalDominant(float)` overloads. Remove the
      large commented-out blocks.
- [ ] Move `tree.hh`, `tree_util.hh`, `anyoption.*` to
      `src/c/third_party/`; exclude from clang-tidy and clang-format. Patch
      the single `using namespace std;` line in `tree.hh` so it stops leaking.
- [ ] Optional: replace `AnyOption` (1,200 lines, `char*` buffers, macro
      constants) with a ~200-line `OptionParser` supporting long/short names,
      flags, repeatable options, and usage text. Only after Phase 3, so the
      swap touches `options.cc` and `main.cpp` alone.

## Phase 7: Quality of life

- [ ] CI (GitHub Actions or a local script): build and `ctest` on macOS and
      Ubuntu with GSL from the package manager. Passing the golden test on
      both confirms cross-platform reproducibility.
- [ ] Write `README` from the CLAUDE.md overview and build sections.
- [ ] Replace the template `src/c/bpopsim.1` with real content, or generate
      it from the usage text.
- [ ] Remove the "Known issues" and "Uncommitted work" sections of
      `CLAUDE.md` as they are resolved.

## Verification at every phase

```
cmake -B build && cmake --build build && ctest --test-dir build --output-on-failure
```

The `tests/golden/` regression test must pass byte-identical unless the
commit says it changes behavior and regenerates the golden files.
