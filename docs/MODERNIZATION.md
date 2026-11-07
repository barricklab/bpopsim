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

## Phase 4: Isolate GSL behind `RandomGenerator` (done 2026-11-07)

- [x] `src/c/bpopsim/random.{h,cc}`: `RandomGenerator` owns a taus2
      `gsl_rng` behind a pimpl (GSL's `gsl_rng` is an anonymous-struct typedef
      and cannot be forward-declared). Methods `Uniform()`, `UniformInt(n)`,
      `Binomial(p, n)`, `Exponential(mean)`, `Seed()`, `seed()`, and static
      `RandomSeed()` (from `std::random_device`, never 0). Argument and return
      types match the GSL calls exactly so streams are unchanged.
- [x] `RandomGenerator&` threaded through `cPopulation`, `cSubpopulation`,
      `cGenotype::AddOneMutation`, `DetermineMutationCategory`. `random.cc` is
      the only file that includes GSL headers; `gsl_matrix.h` and the unused
      `gsl_rng_env_setup()` call are gone. `cSubpopulation::MutantFitness`,
      declared but never defined, removed.
- [x] Seed is `unsigned long` parsed with the matching `from_string<>`;
      0 draws a random seed from the OS. Fixed-seed output unchanged (golden
      passes). The unseeded default changed from a 16-bit truncation of
      `time * pid` to a full-width OS random value; the help text now says
      "0 = choose a random seed and report it".
- [x] `tests/test_random.cc` (7 cases): recorded reference sequence for seed
      12345 (identical from conda-forge and Homebrew GSL 2.8), same seed
      reproduces / different differs, `Seed()` restarts, `UniformInt` range
      and coverage, binomial mean/variance and edge cases, exponential mean,
      `RandomSeed` nonzero and varying.
- [x] Found while doing this: CMake's Release flags define `NDEBUG`, which
      silently removed every `assert` in the simulation, including the
      population-size consistency checks the tests rely on. `CMakeLists.txt`
      now strips `-DNDEBUG` from all build types. Golden still passes, so
      the asserts hold on the reference run.
- [x] 43/43 tests green.

## Phase 5: Mechanical modern C++ cleanup (done 2026-11-07)

`clang-tidy` driven, applied in three groups (modernize, readability,
performance/bugprone) with a full build and test run after each, then
`clang-format`, then a few manual items.

- [x] `clang-tools` added to `environment.yml`. `.clang-tidy` (modernize,
      readability, performance, bugprone with noisy checks disabled),
      `.clang-format` (2-space, 120 columns, braces after functions/classes),
      `.clang-format-ignore` for vendored files. Run on our sources with
      `conda run -p $PWD/env clang-tidy -p build -fix src/c/bpopsim/<file>`
      (`run-clang-tidy` matched zero files here; the per-file form works).
- [x] ~230 fixes applied: range-for, `auto`, `nullptr`, `= default`, default
      member initializers, `using` aliases, `emplace_back`, `const` on
      by-value returns removed, `empty()` over `size() == 0`, `std::max`,
      `'\n'` over `endl`, const& parameters, const member functions,
      consistent parameter names, parentheses on mixed `/` `-` precedence.
- [x] Manual: `#pragma once` in our five headers; hand-written
      `cSubpopulation` copy constructor and empty destructor replaced with
      `= default`; `main()` VLA replaced with a NULL-terminated
      `vector<char*>`; `create_path` reimplemented on
      `std::filesystem::create_directories` and no longer calls `umask(0)`
      (output directories were world-writable before); unused helpers removed
      from `common.h` (`file_exists`, `file_empty`, `fix_flags`, `make_map`,
      `make_list`, `get_keys`, `split_on_any`, `reverse_string`,
      `repeat_char`, `grep`, `shift`, `from_string(bool)`); C headers switched
      to `<c...>` forms; `using namespace std` removed from `cPopulation.h`
      and `common.h`.
- [x] `using namespace std` still reaches everything through vendored
      `tree.hh` (global scope). Patching that line is Phase 6 and will require
      qualifying names in `common.h` and the headers.
- [x] Live code builds warning-free under `-Wall -Wextra`. What remains: a
      dozen `-Wunused-parameter`/`-Wreturn-type`/`-Wunused-variable` warnings
      inside the unreferenced `cPopulation` methods that Phase 6 deletes, and
      `-Wdeprecated-copy-with-user-provided-copy` / `-Wdeprecated-declarations`
      from vendored `tree.hh`, which goes quiet once Phase 6 moves it to a
      `SYSTEM` include directory.
- [x] Golden still byte-identical; 44/44 tests; suite run repeatedly.

Bugs found and fixed along the way (none affect golden output):
- `cReplicateStatistics::total_mutations`, `total_subpopulations_lost`,
  `longest_sweep_time` were never initialized. The sanitizer build showed two
  identical-seed runs reporting 24683 and 76 total mutations. Now `= 0`.
- `AnyOption::init()` left `struct winsize ws` uninitialized and ignored the
  `ioctl` result. With stdin not a terminal (tests, batch jobs, pipelines)
  `terminal_width` was stack garbage, and a small value sent `word_wrap` into
  an infinite loop while the option table was being built. Reproduced at ~7%
  of runs with stdin redirected, on the pre-cleanup code too. Fixed by
  zero-initializing and checking the return value, and `word_wrap` no longer
  jumps backwards when a word is longer than the width. This could hang the
  real program in any non-interactive run.
- `e_colors` left as a plain enum: it indexes the (zero-length, marker-mode
  only) `by_color` array, which Phase 6 removes with the rest of marker mode
  cleanup.

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
