#ifndef bpopsim_options_h
#define bpopsim_options_h

#include "common.h"

namespace bpopsim {

  //! Register every bpopsim command-line option (names, help text, defaults)
  //! on an AnyOption parser. Shared by main() and the test suite so tests can
  //! build the same option set programmatically.
  void DefineBpopsimOptions(AnyOption& options);

}

#endif
