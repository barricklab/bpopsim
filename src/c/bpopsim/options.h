#pragma once

#include "common.h"
#include "cPopulation.h"
#include "option_parser.h"

namespace bpopsim {

  //! Register every bpopsim command-line option (names, help text, defaults)
  //! on an OptionParser parser. Shared by main() and the test suite so tests can
  //! build the same option set programmatically.
  void DefineBpopsimOptions(OptionParser& options);

  //! Convert parsed command-line options into the parameter structs the
  //! simulation runs on. These are the only functions that read OptionParser
  //! values; everything downstream works from the structs.
  //! ParseSimulationParameters() calls Finalize() but not Validate().
  //! ParseOutputParameters() throws std::invalid_argument for a negative --burn-in.
  SimulationParameters ParseSimulationParameters(OptionParser& options);
  OutputParameters ParseOutputParameters(OptionParser& options);
  MarkerDivergence ParseMarkerDivergence(OptionParser& options);

}  // namespace bpopsim
