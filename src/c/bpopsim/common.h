#pragma once

// Shared includes and small string/number helpers used throughout bpopsim.
// These helpers are at global scope and fully qualified; nothing here brings
// std names into scope for the files that include it.

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
#include "tree.hh"
#include "tree_util.hh"
#include "random.h"

//!< Split a string on a delimiter into a vector
inline std::vector<std::string> split(const std::string& theString, const std::string& theDelimiter)
{
  assert(!theDelimiter.empty());  // My own ASSERT macro.

  std::size_t start = 0, end = 0;
  std::vector<std::string> theStringVector;

  if (theString.empty()) return theStringVector;
  // return empty list if the string is empty

  while (end != std::string::npos) {
    end = theString.find(theDelimiter, start);

    // If at end, use length=maxLength.  Else use length=end-start.
    theStringVector.push_back(theString.substr(start, (end == std::string::npos) ? std::string::npos : end - start));

    // If at end, use start=maxSize.  Else use start=end+delimiter.
    start = (end > (std::string::npos - theDelimiter.size())) ? std::string::npos : end + theDelimiter.size();
  }
  return theStringVector;
}

inline std::string join(const std::vector<std::string>& values, const std::string& separator)
{
  if (values.empty()) return "";

  std::string::size_type size = separator.length() * values.size();
  for (const auto& value : values) size += value.size();

  std::string retval;
  retval.reserve(size);
  retval = values[0];
  for (uint32_t i = 1; i < values.size(); i++) retval += separator + values[i];

  return retval;
}

inline std::string chomp(const std::string& str)
{
  std::size_t last = str.find_last_not_of("\n \t");
  return (last == std::string::npos) ? "" : str.substr(0, last + 1);
}

inline std::ostream& operator<<(std::ostream& stream, const std::vector<std::string>& lhs)
{
  stream << join(lhs, ",");
  return stream;
}

template <typename T> inline std::istream& operator>>(std::istream& stream, std::vector<T>& rhs)
{
  rhs.clear();
  // the different values are separated by newlines "\n" in the string
  // so they are read out correctly this way! @JEB
  while (!stream.eof()) {
    T t;
    stream >> std::boolalpha >> t;
    rhs.push_back(t);
  }
  return stream;
}

template <typename T> inline std::string to_string(const T& t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}

inline std::string to_string(const std::pair<int, int>& t)
{
  return to_string(t.first) + '/' + to_string(t.second);
}
inline std::string to_string(const double& t, const uint32_t precision)
{
  if (std::isnan(t)) {
    return "NA";
  } else {
    std::ostringstream interpreter;
    interpreter << std::fixed << std::setprecision(precision) << t;
    return interpreter.str();
  }
}

template <typename T> inline std::string to_string(const std::vector<T>& t)
{
  std::vector<std::string> ts;
  for (auto it = t.begin(); it != t.end(); it++) {
    ts.push_back(to_string(*it));
  }

  return join(ts, ",");
}

template <typename T> inline T from_string(const std::string& s)
{
  assert(!s.empty());
  T t;
  std::istringstream iss(s);
  iss >> std::boolalpha >> t;
  return t;
}

//! Returns true if exp matches anywhere in input
//! true regex expressions ie !^[\s](\w)+ do not work.
inline bool regex_m(const std::string& exp, const std::string& input)
{
  return input.find(exp) != std::string::npos;
}

inline std::string to_upper(const std::string& input)
{
  std::string str = input;
  std::transform(str.begin(), str.end(), str.begin(), ::toupper);
  return str;
}

inline std::string to_lower(const std::string& input)
{
  std::string str = input;
  std::transform(str.begin(), str.end(), str.begin(), ::tolower);
  return str;
}

inline std::string substitute(const std::string& in_s, const std::string& replace_this, const std::string& with_this)
{
  std::string s = in_s;
  std::size_t pos = s.find(replace_this);
  while (pos != std::string::npos) {
    s.replace(pos, replace_this.size(), with_this);
    pos += with_this.size();
    pos = s.find(replace_this, pos);
  }
  return s;
}

// Return the path of the file without the trailing forward-slash
inline std::string dirname(const std::string& file_name)
{
  std::size_t found = file_name.rfind('/');
  return ((found != std::string::npos) ? file_name.substr(0, found) : "");
}

//! Create a directory and any missing parents. Exits the program on failure.
//! Permissions follow the process umask.
inline std::string create_path(const std::string& path)
{
  std::error_code ec;
  std::filesystem::create_directories(path, ec);
  if (ec && !std::filesystem::is_directory(path)) {
    std::cerr << "Could not create path: '" << path << "' (" << ec.message() << ")\n";
    std::exit(-1);
  }
  return path;
}
