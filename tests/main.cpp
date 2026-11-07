#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

// Defined in src/c/bpopsim/main.cpp for the real program; the test binary
// needs its own definition because bpopsim_core references it.
bool g_verbose = false;
