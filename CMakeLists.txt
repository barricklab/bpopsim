cmake_minimum_required(VERSION 3.15)

project(bpopsim
  VERSION 0.1
  DESCRIPTION "Simulator of asexual population evolution under serial transfer"
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# For clangd and other editor tooling: build/compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

# Keep assert() active in every build type. The simulation uses asserts as
# cheap internal consistency checks (population size bookkeeping, tree
# invariants) and the unit tests rely on them firing. CMake's default Release
# flags would define NDEBUG and silently remove them.
foreach(cfg RELEASE RELWITHDEBINFO MINSIZEREL)
  string(REPLACE "-DNDEBUG" "" CMAKE_CXX_FLAGS_${cfg} "${CMAKE_CXX_FLAGS_${cfg}}")
  string(REPLACE "-DNDEBUG" "" CMAKE_C_FLAGS_${cfg} "${CMAKE_C_FLAGS_${cfg}}")
endforeach()

# GSL is the single external dependency. It supplies every random draw and is
# kept (rather than <random>) because its streams are identical on every platform
# for a given seed and version. See CLAUDE.md.
option(BPOPSIM_STATIC_GSL "Link GSL statically (used by binarydist.sh)" OFF)
if(BPOPSIM_STATIC_GSL)
  set(_saved_suffixes "${CMAKE_FIND_LIBRARY_SUFFIXES}")
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif()
find_package(GSL REQUIRED)
if(BPOPSIM_STATIC_GSL)
  set(CMAKE_FIND_LIBRARY_SUFFIXES "${_saved_suffixes}")
endif()

set(BPOPSIM_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/c/bpopsim")
set(BPOPSIM_THIRD_PARTY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/c/third_party")

# Core simulation code as a library so tests can link against it without main().
add_library(bpopsim_core STATIC
  "${BPOPSIM_SRC_DIR}/cPopulation.cc"
  "${BPOPSIM_SRC_DIR}/cSubpopulation.cc"
  "${BPOPSIM_SRC_DIR}/options.cc"
  "${BPOPSIM_SRC_DIR}/random.cc"
  # Vendored third-party source (tree.hh and tree_util.hh are header-only)
  "${BPOPSIM_THIRD_PARTY_DIR}/anyoption.cpp"
)
target_include_directories(bpopsim_core PUBLIC "${BPOPSIM_SRC_DIR}")
# SYSTEM so warnings inside the vendored headers are not reported against our code.
target_include_directories(bpopsim_core SYSTEM PUBLIC "${BPOPSIM_THIRD_PARTY_DIR}")
target_link_libraries(bpopsim_core PUBLIC GSL::gsl GSL::gslcblas)
target_compile_options(bpopsim_core PRIVATE
  $<$<COMPILE_LANG_AND_ID:CXX,Clang,AppleClang,GNU>:-Wall -Wextra>
)
# Silence warnings inside the vendored source file only.
set_source_files_properties("${BPOPSIM_THIRD_PARTY_DIR}/anyoption.cpp" PROPERTIES COMPILE_OPTIONS "-w")

add_executable(bpopsim "${BPOPSIM_SRC_DIR}/main.cpp")
target_link_libraries(bpopsim PRIVATE bpopsim_core)
target_compile_options(bpopsim PRIVATE
  $<$<COMPILE_LANG_AND_ID:CXX,Clang,AppleClang,GNU>:-Wall -Wextra>
)

include(GNUInstallDirs)
install(TARGETS bpopsim RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/c/bpopsim.1"
  DESTINATION "${CMAKE_INSTALL_MANDIR}/man1")

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
