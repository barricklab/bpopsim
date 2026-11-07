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

## Phase 2: Test harness (done 2026-11-07)

- [x] doctest 2.5.3 vendored unmodified in `tests/third_party/doctest/`
      (`doctest.h`, `doctest.cmake`, `doctestAddTests.cmake`, `LICENSE.txt`,
      plus a `README` recording the source tag).
- [x] `tests/CMakeLists.txt`: `bpopsim_tests` executable from `tests/main.cpp`
      (`DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN`, defines `g_verbose`) plus
      `test_*.cc`, linked to `bpopsim_core`, registered per test case with
      `doctest_discover_tests`.
- [x] `src/c/bpopsim/options.{h,cc}`: the option table moved out of `main()`
      into `DefineBpopsimOptions(AnyOption&)` so tests can build the same
      parser programmatically. Behavior unchanged; golden still passes.
- [x] `tests/test_common.cc` (15 cases): `split`, `join`, `from_string`,
      `from_string<vector<double>>` (newline-separated, as AnyOption joins
      repeated options), `to_string(double, precision)`, `substitute`,
      `to_upper`/`to_lower`, `chomp`, `dirname`, `regex_m`.
- [x] `chomp()` off-by-one fixed (`-1` dropped a real character; now handles
      an all-whitespace string too). It had no callers, so no output change.
- [x] `tests/test_parameters.cc` (7 cases): defaults, dilution factor and
      binomial probability from `-T`, category rates/fractions from repeated
      `-u`/`-s`, `DetermineMutationCategory` sampling proportions over 1e5
      draws, `-m`/`-1`/`-f` parsing, and `cGenotype::AddOneMutation` under the
      additive, multiplicative, one-time (`o`), and first-mutation-override
      models. The `-u`/`-s` count mismatch calls `exit(-1)` so it is not
      testable until Phase 3 turns it into an exception.
- [x] Regression test: done in Phase 1 (`tests/golden_test.cmake`).
- [ ] Add a second golden case in marker-divergence mode (`-k 2`) once the
      `by_color[0]` zero-length array in `TransferResampleDilution` is fixed.
- [x] `ctest --test-dir build --output-on-failure`: 23/23 green.

## Phase 3: Decouple parameters from option parsing (done 2026-11-07)

Biggest structural win. Makes the simulator constructible in tests.

- [x] `SimulationParameters`, `OutputParameters`, `MarkerDivergence` are plain
      structs with default member initializers equal to the command-line
      defaults, no `AnyOption&` constructors. `OutputParameters` and
      `MarkerDivergence` moved from private nested types of `cPopulation` to
      namespace scope. `cStatistics` takes `const OutputParameters&`.
- [x] `SimulationParameters::Finalize()` computes the derived fields
      (dilution factor, binomial p, total rate, category fractions) and can be
      re-run after changing inputs. The default constructor calls it.
- [x] `SimulationParameters::Validate()` throws `std::invalid_argument` on a
      `-u`/`-s` count mismatch. `main` catches it, appends the message to the
      usage text, prints usage, returns -1 (same exit code as before).
- [x] `options.{h,cc}`: `ParseSimulationParameters`, `ParseOutputParameters`,
      `ParseMarkerDivergence` are the only code that reads `AnyOption` values.
      `cPopulation` is constructed from the three structs plus rng and
      replicate number. `main` prints settings with
      `simulation_parameters.Print()` instead of constructing a throwaway
      `cPopulation`.
- [x] Global `g_verbose` removed. `OutputParameters::verbose` and `::debug`
      replace it; `cSubpopulation::CreateDescendant` takes a `verbose`
      argument. Output identical.
- [x] `tests/test_parameters.cc` rewritten (11 cases): direct construction,
      `Finalize` re-run, `Validate` throwing, parse functions for all three
      structs including `-v`/`--debug`/`-k`/`-x`/`-t`, category sampling,
      `AddOneMutation` models, child inheritance.
- [x] `tests/test_population.cc` (8 cases) runs `cPopulation` end to end with
      `debug` on (internal population-size asserts active): record counts,
      fitness monotone and starting at initial fitness, genotype frequencies
      summing to one, ancestral clade at frequency one, mutations accumulate,
      same seed reproduces and different seed differs, fitness end condition,
      burn-in and coarse-graining record counts.
- [x] Golden test still byte-identical. 36/36 tests green.

Learned while testing: burn-in transfers are performed in addition to `-n`;
the run continues until the transfer counter (which starts at `-burn_in`)
reaches `-n`, so `-n 6 -b 3` performs 9 transfers and records 0..6.

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
