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

Toolchain comes from the local conda env in `environment.yml` (see
CLAUDE.md): conda-forge clang 21, cmake 4.4, GSL 2.8 pinned. Phase 1 will use
`cmake_minimum_required(VERSION 3.15)` so the project also builds on Ubuntu LTS
for Phase 7 CI. System copies (Apple clang 17, MacPorts cmake 3.11, a 2009
x86_64 GSL 1.13 in `/usr/local`) are not used.

Cross-toolchain check already done: the golden output is byte-identical when
built with Apple clang 17 + Homebrew GSL 2.8 and with conda clang 21 + conda
GSL 2.8.

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
  - `icsilog.h` used the `register` storage class, a hard error under the
    C++17 default of conda's clang 21. Removed (it is a no-op hint).
  - Environment: the only GSL present was x86_64 1.13 from 2009. Added
    `environment.yml` (conda-forge: cxx-compiler, cmake, autotools,
    gsl=2.8) and build inside `./env` as breseq does.
- [x] Golden output committed under `tests/golden/`, command in
      `tests/golden/COMMAND`. Runs in well under a second, two replicates,
      30 transfers, every working output option. Verified byte-identical
      across two runs. `--output-largest-sweep-size` is excluded because it
      is `assert(false)`.

## Phase 1: CMake build, remove autotools (done 2026-11-07)

- [x] Top-level `CMakeLists.txt`: C++17, extensions off, Release default,
      `find_package(GSL REQUIRED)`, `bpopsim_core` static library
      (`cPopulation.cc`, `cSubpopulation.cc`, vendored `anyoption.cpp` and
      `icsilog.c`), `bpopsim` executable, `-Wall -Wextra` on our code with
      vendored sources compiled `-w`, install rules for the binary and man
      page, `include(CTest)`. `BPOPSIM_STATIC_GSL` option for release builds.
- [x] `tests/CMakeLists.txt` + `tests/golden_test.cmake`: CTest
      `golden_regression` runs the command in `tests/golden/COMMAND` and
      compares every file byte for byte, also failing on missing or extra
      output files. (Pulled forward from Phase 2 because it is Phase 1's
      verification.)
- [x] Deleted `configure.ac`, `Makefile.am`, `src/c/bpopsim/Makefile.am`,
      `bootstrap.sh`, `aux_build/`, `src/c/bpopsim/compile`,
      `src/c/bpopsim.xcodeproj/`. Generate an Xcode project with
      `cmake -G Xcode` when wanted.
- [x] `binarydist.sh` rewritten around `cmake --build`, `ctest`,
      `cmake --install` with static GSL. On macOS it follows breseq: system
      Apple clang, GSL built from source (version from the `environment.yml`
      pin) as a static universal arm64 + x86_64 archive in
      `binarydist_build/`, deployment target 10.13. The tarball binary links
      only `/usr/lib/libc++` and `libSystem`. Verified that both the arm64
      and the x86_64 (under Rosetta) slices reproduce `tests/golden/`.
- [x] `.gitignore` reduced to CMake, dist, Xcode-generated, Perl, editor
      entries. `environment.yml` autotools block removed.
- [x] `CLAUDE.md` build section updated. `README` is still empty (Phase 7).
- [x] Verified: `ctest` golden regression passes with the CMake build.

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
- [x] Regression test: done in Phase 1 (`tests/golden_test.cmake`).
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
