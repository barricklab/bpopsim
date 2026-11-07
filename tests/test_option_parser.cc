#include "doctest.h"
#include "option_parser.h"

#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;
using namespace bpopsim;

namespace {

  // Build a parser with a representative option set and parse the given arguments.
  struct Parsed
  {
    explicit Parsed(const vector<string>& args)
    : parser("Usage: test")
    {
      parser("help,h", "Print help", NoArgument);
      parser("verbose,v", "Be chatty", NoArgument);
      parser("debug", "Extra checks", NoArgument);
      parser("rate,u", "Mutation rate", "1E-7");
      parser("effect,s", "Fitness effect", 0.1);
      parser("transfers,n", "Number of transfers", 50);
      parser("model,f", "Effect model", 'u');
      parser("generations,T", "Generations per transfer", 6.64);
      parser("first,1", "First mutation effects");
      parser("output-folder,o", "Output folder", ".");

      storage.emplace_back("prog");
      for (const string& a : args) storage.push_back(a);
      for (string& a : storage) argv.push_back(const_cast<char*>(a.c_str()));
      argv.push_back(nullptr);
      parser.Parse(static_cast<int>(storage.size()), argv.data());
    }
    OptionParser parser;
    vector<string> storage;
    vector<char*> argv;
  };

}  // namespace

TEST_CASE("defaults are formatted the way the simulation expects")
{
  Parsed p({});
  CHECK(p.parser["rate"] == "1E-7");
  CHECK(p.parser["effect"] == "0.1");
  CHECK(p.parser["transfers"] == "50");
  CHECK(p.parser["model"] == "u");
  CHECK(p.parser["generations"] == "6.64");
  CHECK(p.parser["output-folder"] == ".");
  CHECK(p.parser["first"] == "");  // no default
  CHECK_FALSE(p.parser.count("rate"));
  CHECK_FALSE(p.parser.count("help"));
  CHECK(p.parser["help"] == "");
}

TEST_CASE("long options with a separate value or with =")
{
  Parsed p({"--rate", "2E-7", "--transfers=10", "--output-folder", "out dir"});
  CHECK(p.parser.count("rate"));
  CHECK(p.parser["rate"] == "2E-7");
  CHECK(p.parser["transfers"] == "10");
  CHECK(p.parser["output-folder"] == "out dir");
}

TEST_CASE("short options attached or detached")
{
  Parsed p({"-u", "3E-7", "-n20", "-f", "e"});
  CHECK(p.parser["rate"] == "3E-7");
  CHECK(p.parser["transfers"] == "20");
  CHECK(p.parser["model"] == "e");
}

TEST_CASE("flags, alone and combined")
{
  Parsed p({"-hv", "--debug"});
  CHECK(p.parser.count("help"));
  CHECK(p.parser.count("verbose"));
  CHECK(p.parser.count("debug"));
  CHECK(p.parser["help"] == "true");
  CHECK(p.parser["verbose"] == "true");
}

TEST_CASE("a combined run of flags may end in an option that takes the rest")
{
  Parsed p({"-vn7", "-hT", "3"});
  CHECK(p.parser.count("verbose"));
  CHECK(p.parser["transfers"] == "7");
  CHECK(p.parser.count("help"));
  CHECK(p.parser["generations"] == "3");
}

TEST_CASE("repeated options accumulate newline-separated values")
{
  Parsed p({"-u", "1E-7", "-u", "2E-7", "--rate=3E-7", "-1", "0.05", "-1", "0.02"});
  CHECK(p.parser["rate"] == "1E-7\n2E-7\n3E-7");
  CHECK(p.parser["first"] == "0.05\n0.02");
}

TEST_CASE("negative numbers are accepted as values")
{
  Parsed p({"-s", "-0.5", "--effect=-0.25"});
  CHECK(p.parser["effect"] == "-0.5\n-0.25");
}

TEST_CASE("errors are reported as std::invalid_argument")
{
  CHECK_THROWS_AS(Parsed({"--bogus"}), std::invalid_argument);
  CHECK_THROWS_AS(Parsed({"-z"}), std::invalid_argument);
  CHECK_THROWS_AS(Parsed({"--rate"}), std::invalid_argument);      // missing value
  CHECK_THROWS_AS(Parsed({"-u"}), std::invalid_argument);          // missing value
  CHECK_THROWS_AS(Parsed({"--help=yes"}), std::invalid_argument);  // flag with value
  CHECK_THROWS_AS(Parsed({"stray"}), std::invalid_argument);       // positional
  CHECK_THROWS_WITH(Parsed({"--bogus"}), "Unknown option: --bogus");
  CHECK_THROWS_WITH(Parsed({"--rate"}), "Option --rate requires a value");
}

TEST_CASE("asking for an unregistered option is a programming error")
{
  Parsed p({});
  CHECK_THROWS_AS(p.parser.count("nope"), std::logic_error);
  CHECK_THROWS_AS(p.parser["nope"], std::logic_error);
}

TEST_CASE("registering duplicate or malformed names is a programming error")
{
  OptionParser parser;
  parser("alpha,a", "first", NoArgument);
  CHECK_THROWS_AS(parser("alpha", "again", NoArgument), std::logic_error);
  CHECK_THROWS_AS(parser("beta,a", "same short", NoArgument), std::logic_error);
  CHECK_THROWS_AS(parser("gamma,gg", "two-char short", NoArgument), std::logic_error);
  CHECK_THROWS_AS(parser("", "empty", NoArgument), std::logic_error);
}

TEST_CASE("usage text lists options in order with sections, defaults and wrapping")
{
  OptionParser parser("Usage: prog [options]");
  parser("help,h", "Print this help", NoArgument);
  parser.AddUsage("");
  parser.AddUsage("==== Section ====");
  parser("rate,u",
         "A description long enough that it certainly has to wrap onto a second line when the "
         "terminal is eighty columns wide, which is the width used when output is not a terminal.",
         "1E-7");
  parser("a-very-long-option-name-that-overflows", "Short text", NoArgument);
  parser.AddUsage("Trailing note.");

  ostringstream out;
  parser.PrintUsage(out);
  string text = out.str();

  CHECK(text.find("Usage: prog [options]") != string::npos);
  CHECK(text.find("==== Section ====") != string::npos);
  CHECK(text.find("-h,--help") != string::npos);
  CHECK(text.find("-u,--rate <arg>") != string::npos);
  CHECK(text.find("(DEFAULT=1E-7)") != string::npos);
  CHECK(text.find("--a-very-long-option-name-that-overflows") != string::npos);
  CHECK(text.find("Trailing note.") != string::npos);
  // Order: section heading before --rate, trailing note last.
  CHECK(text.find("==== Section ====") < text.find("--rate"));
  CHECK(text.find("Trailing note.") > text.find("--a-very-long"));
  // Every line fits in 80 columns.
  istringstream lines(text);
  string line;
  while (getline(lines, line)) CHECK(line.size() <= 80);
}
