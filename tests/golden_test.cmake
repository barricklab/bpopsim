# Run bpopsim with the arguments in ${GOLDEN_DIR}/COMMAND, writing into ${OUT_DIR},
# then compare every file in ${GOLDEN_DIR} (except COMMAND) byte for byte.
#
# Invoked by CTest as:
#   cmake -DBPOPSIM=<exe> -DGOLDEN_DIR=<dir> -DOUT_DIR=<dir> -P golden_test.cmake

foreach(var BPOPSIM GOLDEN_DIR OUT_DIR)
  if(NOT DEFINED ${var})
    message(FATAL_ERROR "golden_test.cmake: -D${var}=... is required")
  endif()
endforeach()

# The last non-comment, non-blank line of COMMAND is the invocation.
file(STRINGS "${GOLDEN_DIR}/COMMAND" command_lines)
set(command_line "")
foreach(line IN LISTS command_lines)
  if(NOT line MATCHES "^[ \t]*#" AND NOT line MATCHES "^[ \t]*$")
    set(command_line "${line}")
  endif()
endforeach()
if(command_line STREQUAL "")
  message(FATAL_ERROR "No command line found in ${GOLDEN_DIR}/COMMAND")
endif()

# Drop the leading program name; we supply the real executable path.
string(REGEX REPLACE "^[ \t]*bpopsim[ \t]+" "" args "${command_line}")
separate_arguments(arg_list UNIX_COMMAND "${args}")

file(REMOVE_RECURSE "${OUT_DIR}")
file(MAKE_DIRECTORY "${OUT_DIR}")

execute_process(
  COMMAND "${BPOPSIM}" ${arg_list} -o "${OUT_DIR}"
  RESULT_VARIABLE run_result
  OUTPUT_VARIABLE run_stdout
  ERROR_VARIABLE run_stderr
)
if(NOT run_result EQUAL 0)
  message(FATAL_ERROR
    "bpopsim exited with ${run_result}\n"
    "command: ${BPOPSIM} ${args} -o ${OUT_DIR}\n"
    "stdout:\n${run_stdout}\nstderr:\n${run_stderr}")
endif()

file(GLOB golden_files RELATIVE "${GOLDEN_DIR}" "${GOLDEN_DIR}/*")
list(REMOVE_ITEM golden_files COMMAND)
list(SORT golden_files)

file(GLOB out_files RELATIVE "${OUT_DIR}" "${OUT_DIR}/*")
list(SORT out_files)

set(failures "")

foreach(f IN LISTS golden_files)
  if(NOT EXISTS "${OUT_DIR}/${f}")
    list(APPEND failures "missing output file: ${f}")
    continue()
  endif()
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E compare_files "${GOLDEN_DIR}/${f}" "${OUT_DIR}/${f}"
    RESULT_VARIABLE cmp_result
    OUTPUT_QUIET ERROR_QUIET
  )
  if(NOT cmp_result EQUAL 0)
    list(APPEND failures "differs from golden: ${f}")
  endif()
endforeach()

foreach(f IN LISTS out_files)
  if(NOT f IN_LIST golden_files)
    list(APPEND failures "unexpected output file: ${f}")
  endif()
endforeach()

if(failures)
  list(JOIN failures "\n  " failure_text)
  message(FATAL_ERROR
    "Golden regression failed:\n  ${failure_text}\n"
    "golden: ${GOLDEN_DIR}\noutput: ${OUT_DIR}\n"
    "If the change is intentional, regenerate golden/ and say so in the commit.")
endif()

list(LENGTH golden_files n)
message(STATUS "Golden regression passed: ${n} files identical")
