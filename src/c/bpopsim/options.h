#pragma once

#include "common.h"
#include "cPopulation.h"

namespace bpopsim {

  //! Register every bpopsim command-line option (names, help text, defaults)
  //! on an AnyOption parser. Shared by main() and the test suite so tests can
  //! build the same option set programmatically.
  void DefineBpopsimOptions(AnyOption& options);

  //! Convert parsed command-line options into the parameter structs the
  //! simulation runs on. These are the only functions that read AnyOption
  //! values; everything downstream works from the structs.
  //! ParseSimulationParameters() calls Finalize() but not Validate().
  SimulationParameters ParseSimulationParameters(AnyOption& options);
  OutputParameters ParseOutputParameters(AnyOption& options);
  MarkerDivergence ParseMarkerDivergence(AnyOption& options);

}  // namespace bpopsim
