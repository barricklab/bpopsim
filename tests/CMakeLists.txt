# ---- Unit tests (doctest, vendored in third_party/doctest) ----
include("${CMAKE_CURRENT_SOURCE_DIR}/third_party/doctest/doctest.cmake")

add_executable(bpopsim_tests
  main.cpp
  test_common.cc
  test_option_parser.cc
  test_parameters.cc
  test_population.cc
  test_random.cc
)
target_include_directories(bpopsim_tests SYSTEM PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/third_party/doctest")
target_link_libraries(bpopsim_tests PRIVATE bpopsim_core)
target_compile_options(bpopsim_tests PRIVATE
  $<$<COMPILE_LANG_AND_ID:CXX,Clang,AppleClang,GNU>:-Wall -Wextra>
)

# Registers each TEST_CASE as its own CTest test at build time.
doctest_discover_tests(bpopsim_tests)

# ---- Golden regression test ----
# Run the fixed-seed command recorded in golden/COMMAND and require every
# output file to be byte-identical to golden/.
add_test(NAME golden_regression
  COMMAND "${CMAKE_COMMAND}"
    -DBPOPSIM=$<TARGET_FILE:bpopsim>
    -DGOLDEN_DIR=${CMAKE_CURRENT_SOURCE_DIR}/golden
    -DOUT_DIR=${CMAKE_CURRENT_BINARY_DIR}/golden_out
    -P "${CMAKE_CURRENT_SOURCE_DIR}/golden_test.cmake"
)
