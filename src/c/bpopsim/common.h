#pragma once

// Shared includes and small string/number helpers used throughout bpopsim.

#include <cassert>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// Vendored third-party headers (src/c/third_party, a SYSTEM include directory).
// tree.hh has `using namespace std;` at global scope, which is what lets this
// file and the headers that include it use unqualified std names.
#include "tree.hh"
#include "tree_util.hh"
#include "random.h"

//!< Split a string on a delimiter into a vector
inline vector<string> split(const string& theString, const string& theDelimiter)
{
  assert(!theDelimiter.empty());  // My own ASSERT macro.
  assert(theString.size() != 1);  //@GRC check that programmer hasn't switched variables in declaration.

  size_t start = 0, end = 0;
  vector<string> theStringVector;

  if (theString.empty()) return theStringVector;
  // return empty list if the string is empty

  while (end != string::npos) {
    end = theString.find(theDelimiter, start);

    // If at end, use length=maxLength.  Else use length=end-start.
    theStringVector.push_back(theString.substr(start, (end == string::npos) ? string::npos : end - start));

    // If at end, use start=maxSize.  Else use start=end+delimiter.
    start = (end > (string::npos - theDelimiter.size())) ? string::npos : end + theDelimiter.size();
  }
  return theStringVector;
}

inline string join(const vector<string>& values, const string& separator)
{
  if (values.empty()) return "";

  string::size_type size = separator.length() * values.size();
  for (const auto& value : values) size += value.size();

  string retval;
  retval.reserve(size);
  retval = values[0];
  for (uint32_t i = 1; i < values.size(); i++) retval += separator + values[i];

  return retval;
}

inline string chomp(const string& str)
{
  size_t last = str.find_last_not_of("\n \t");
  return (last == string::npos) ? "" : str.substr(0, last + 1);
}

inline ostream& operator<<(ostream& stream, const vector<string>& lhs)
{
  stream << join(lhs, ",");
  return stream;
}

template <typename T> inline istream& operator>>(istream& stream, vector<T>& rhs)
{
  rhs.clear();
  string value;
  // the different values are separated by newlines "\n" in the string
  // so they are read out correctly this way! @JEB
  while (!stream.eof()) {
    T t;
    stream >> boolalpha >> t;
    rhs.push_back(t);
  }
  return stream;
}

template <typename T> inline string to_string(const T& t)
{
  stringstream ss;
  ss << t;
  return ss.str();
}

inline string to_string(const pair<int, int>& t)
{
  return to_string(t.first) + '/' + to_string(t.second);
}
inline string to_string(const double& t, const uint32_t precision)
{
  if (isnan(t)) {
    return "NA";
  } else {
    ostringstream interpreter;
    interpreter << fixed << setprecision(precision) << t;
    return interpreter.str();
  }
}

template <typename T> inline string to_string(const vector<T>& t)
{
  vector<string> ts;
  for (auto it = t.begin(); it != t.end(); it++) {
    ts.push_back(to_string(*it));
  }

  return join(ts, ",");
}

template <typename T> inline T from_string(const string& s)
{
  assert(!s.empty());
  T t;
  istringstream iss(s);
  iss >> boolalpha >> t;
  return t;
}

//! Returns true if exp matches anywhere in input
//! true regex expressions ie !^[\s](\w)+ do not work.
inline bool regex_m(const string& exp, const string& input)
{
  return input.find(exp) != string::npos;
}

inline string to_upper(const string& input)
{
  string str = input;
  transform(str.begin(), str.end(), str.begin(), ::toupper);
  return str;
}

inline string to_lower(const string& input)
{
  string str = input;
  transform(str.begin(), str.end(), str.begin(), ::tolower);
  return str;
}

inline string substitute(const string& in_s, const string& replace_this, const string& with_this)
{
  string s = in_s;
  size_t pos = s.find(replace_this);
  while (pos != string::npos) {
    s.replace(pos, replace_this.size(), with_this);
    pos += with_this.size();
    pos = s.find(replace_this, pos);
  }
  return s;
}

// Return the path of the file without the trailing forward-slash
inline string dirname(const string& file_name)
{
  size_t found = file_name.rfind('/');
  return ((found != string::npos) ? file_name.substr(0, found) : "");
}

//! Create a directory and any missing parents. Exits the program on failure.
//! Permissions follow the process umask (the old implementation cleared the
//! umask first, which made every output directory world-writable).
inline string create_path(const string& path)
{
  std::error_code ec;
  std::filesystem::create_directories(path, ec);
  if (ec && !std::filesystem::is_directory(path)) {
    cerr << "Could not create path: '" << path << "' (" << ec.message() << ")\n";
    exit(-1);
  }
  return path;
}
