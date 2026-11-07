# Regression test: run the fixed-seed command recorded in golden/COMMAND and
# require every output file to be byte-identical to golden/.
add_test(NAME golden_regression
  COMMAND "${CMAKE_COMMAND}"
    -DBPOPSIM=$<TARGET_FILE:bpopsim>
    -DGOLDEN_DIR=${CMAKE_CURRENT_SOURCE_DIR}/golden
    -DOUT_DIR=${CMAKE_CURRENT_BINARY_DIR}/golden_out
    -P "${CMAKE_CURRENT_SOURCE_DIR}/golden_test.cmake"
)
