#include "doctest.h"
#include "common.h"

#include <cmath>
#include <string>
#include <vector>

using namespace std;

TEST_CASE("split on a single-character delimiter")
{
  vector<string> parts = split("a,b,c", ",");
  REQUIRE(parts.size() == 3);
  CHECK(parts[0] == "a");
  CHECK(parts[1] == "b");
  CHECK(parts[2] == "c");
}

TEST_CASE("split on a multi-character delimiter")
{
  vector<string> parts = split("one::two::three", "::");
  REQUIRE(parts.size() == 3);
  CHECK(parts[1] == "two");
}

TEST_CASE("split of an empty string is empty")
{
  CHECK(split("", ",").empty());
}

TEST_CASE("split keeps empty fields")
{
  vector<string> parts = split("a,,b", ",");
  REQUIRE(parts.size() == 3);
  CHECK(parts[1] == "");
}

TEST_CASE("join")
{
  vector<string> v = {"x", "y", "z"};
  CHECK(join(v, "\t") == "x\ty\tz");
  CHECK(join(vector<string>(), ",") == "");
  CHECK(join(vector<string>{"only"}, ",") == "only");
}

TEST_CASE("split and join round-trip")
{
  string s = "1E-7\n2E-7\n3E-7";
  CHECK(join(split(s, "\n"), "\n") == s);
}

TEST_CASE("from_string parses numbers including scientific notation")
{
  CHECK(from_string<double>("1E5") == doctest::Approx(100000.0));
  CHECK(from_string<double>("6.64") == doctest::Approx(6.64));
  CHECK(from_string<int32_t>("42") == 42);
  CHECK(from_string<uint32_t>("0") == 0);
}

TEST_CASE("from_string<vector<double>> reads newline-separated values")
{
  // AnyOption joins repeated options with "\n"; this is how -u/-s categories
  // reach SimulationParameters.
  vector<double> v = from_string<vector<double> >("1E-7\n2E-7");
  REQUIRE(v.size() == 2);
  CHECK(v[0] == doctest::Approx(1e-7));
  CHECK(v[1] == doctest::Approx(2e-7));

  vector<double> single = from_string<vector<double> >("0.1");
  REQUIRE(single.size() == 1);
  CHECK(single[0] == doctest::Approx(0.1));
}

TEST_CASE("to_string with precision")
{
  CHECK(to_string(3.14159, 2) == "3.14");
  CHECK(to_string(2.0, 3) == "2.000");
  CHECK(to_string(numeric_limits<double>::quiet_NaN(), 2) == "NA");
}

TEST_CASE("to_string of a vector is comma-joined")
{
  vector<int> v = {1, 2, 3};
  CHECK(to_string(v) == "1,2,3");
}

TEST_CASE("substitute replaces every occurrence")
{
  CHECK(substitute("a-b-c", "-", "+") == "a+b+c");
  CHECK(substitute("abc", "x", "y") == "abc");
  // Replacement containing the pattern must not loop forever.
  CHECK(substitute("a.b", ".", "..") == "a..b");
}

TEST_CASE("to_upper and to_lower")
{
  CHECK(to_upper("convert-tree") == "CONVERT-TREE");
  CHECK(to_lower("ABC") == "abc");
}

TEST_CASE("chomp strips trailing whitespace and newlines")
{
  CHECK(chomp("abc \n") == "abc");
  CHECK(chomp("abc") == "abc");
  CHECK(chomp("abc\t\n\n") == "abc");
  CHECK(chomp("   ") == "");
  CHECK(chomp("") == "");
}

TEST_CASE("dirname")
{
  CHECK(dirname("a/b/c.txt") == "a/b");
  CHECK(dirname("c.txt") == "");
  CHECK(dirname("/abs/path/") == "/abs/path");
}

TEST_CASE("regex_m is a plain substring match")
{
  CHECK(regex_m("needle", "hay needle stack"));
  CHECK_FALSE(regex_m("needle", "haystack"));
}
