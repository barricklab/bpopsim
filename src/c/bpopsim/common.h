#ifndef common_h
#define common_h

// System headers

// C
#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

// C++
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip> 
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <limits>

//Other Includes
#include "tree.hh"
#include "tree_util.hh"
#include "icsilog.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>

using namespace std;

// Global variable for keeping track of verbosity
extern bool g_verbose;


// Utility functions
inline bool file_exists(const char *filename)
{
  ifstream ifile(filename);
  return ifile;
}

inline bool file_empty(const char *filename)
{
  ifstream ifile(filename);
  string test;
  getline(ifile, test);
  return !ifile.eof();
}

inline uint32_t fix_flags(uint32_t flags)
{
  flags = ((flags >> 9) << 9) + flags % 128;
  return flags;
}

template <typename T, typename U> struct make_map : public map<T,U>
{
public:
  make_map(const T& key, const U& val) { (*this)(key, val); }
  make_map<T, U>& operator()(const T& key, const U& val)
  {
    (*this)[key] = val;
    return *this;
  }
};

template <typename T> struct make_list : public vector<T>
{
public:
  make_list(const T& t) { (*this)(t); }
  make_list& operator()(const T& t) {
    this->push_back(t);
    return *this;
  }
};


template <class T, class U> inline vector<T> get_keys(const map<T,U>& input)
{
  vector<T> retval;
  for (class map<T,U>::const_iterator it = input.begin(); it != input.end(); it++)
    retval.push_back(it->first);
  return retval;
}

//!< Split a string on a delimiter into a vector
inline vector<string> split(
                            const  string  & theString,
                            const  string  & theDelimiter
                            ) {
  assert(theDelimiter.size() > 0); // My own ASSERT macro.
  assert(theString.size() != 1); //@GRC check that programmer hasn't switched variables in declaration.
  
  size_t start = 0, end = 0;
  vector<string> theStringVector;
  
  if (theString.size() == 0) return theStringVector;
  // return empty list if the string is empty
  
  while (end != string::npos)
  {
    end = theString.find( theDelimiter, start );
    
    // If at end, use length=maxLength.  Else use length=end-start.
    theStringVector.push_back(
                              theString.substr(
                                               start,
                                               (end == string::npos) ? string::npos : end - start
                                               )
                              );
    
    // If at end, use start=maxSize.  Else use start=end+delimiter.
    start =
    (end > (string::npos - theDelimiter.size()))
    ? string::npos
    : end + theDelimiter.size();
  }
  return theStringVector;
}

//!< Split a string on any char in string of delimiters into a vector
inline vector<string> split_on_any(
                                   const  string  & theString,
                                   const  string  & theDelimiters
                                   ) {
  assert(theDelimiters.size() > 0); // My own ASSERT macro.
  
  size_t start = 0, end = 0;
  vector<string> theStringVector;
  
  while (end != string::npos)
  {
    end = theString.find_first_of( theDelimiters, start );
    
    // If at end, use length=maxLength.  Else use length=end-start.
    theStringVector.push_back(
                              theString.substr(
                                               start,
                                               (end == string::npos) ? string::npos : end - start
                                               )
                              );
    
    // If at end, use start=maxSize.  Else use start=end+delimiter.
    start =
    (end > (string::npos - 1))
    ? string::npos
    : end + 1;
  }
  return theStringVector;
}

inline string join(const vector<string>& values, const string& separator)
{
  if(values.size() == 0)
    return "";

  string::size_type size = separator.length() * values.size();
  for(uint32_t i=0; i < values.size(); i++)
    size += values[i].size();

  string retval;
  retval.reserve(size);
  retval = values[0];
  for(uint32_t i = 1; i < values.size(); i++)
    retval += separator + values[i];

  return retval;
}

inline string chomp(const string& str)
{
  return str.substr(0, str.find_last_not_of("\n \t")-1);
}


inline ostream &operator << (ostream &stream, vector<string> lhs)
{
  stream << join(lhs, ",");
  return stream;
}

template <typename T> inline istream &operator >> (istream &stream, vector<T>& rhs)
{
  rhs.clear();
  string value;
  // the different values are separated by newlines "\n" in the string
  // so they are read out correctly this way! @JEB
  while (!stream.eof())
  {
    T t;
    stream >> boolalpha >> t;
    rhs.push_back(t);
  }
  return stream;
}


template <typename T> inline string to_string (const T& t)
{
  stringstream ss;
  ss << t;
  return ss.str();
}

inline string to_string (const pair<int,int>& t)
{
  return to_string(t.first) + '/' + to_string(t.second);
}
inline string to_string (const double& t, const uint32_t precision)
{
  if(isnan(t)) {
    return "NA";
  } else {
    ostringstream interpreter;
    interpreter << fixed << setprecision(precision) << t;
    return interpreter.str();
  }
}


template <typename T> inline string to_string (const vector<T>& t)
{
  vector<string> ts;
  for (typename vector<T>::const_iterator it = t.begin(); it!=t.end(); it++) {
    ts.push_back(to_string(*it));
  }
  
  return join(ts, ",");
}


// handle bool as either TRUE/FALSE or zero/non-zero number
// Does not handle single-character T/F correctly 
// @JEB this is never called
inline bool from_string(const string& s)
{
  bool t = false;
  istringstream iss1(s);
  iss1 >> boolalpha >> t;
  
  int32_t t2;
  istringstream iss2(s);
  iss2 >> noboolalpha >> t2;
  t = t || (t2 != 0);
  
  return t;
}

template <typename T> inline T from_string(const string &s)
{
  assert(!s.empty());
  T t;
  istringstream iss(s);
  iss >> boolalpha >> t;
  return t;
}

//! Returns true if exp matches anywhere in input
//! true regex expressions ie !^[\s](\w)+ do not work.
inline	bool regex_m(string exp, string input)
{
  if(input.find(exp) != string::npos)
    return true;
  else
    return false;
}

inline string to_upper(const string& input)
{
  string str = input;
  transform(str.begin(), str.end(),str.begin(), ::toupper);
  return str;
}

inline string to_lower(const string& input)
{
  string str = input;
  transform(str.begin(), str.end(),str.begin(), ::tolower);
  return str;
}

inline string reverse_string(const string &in_string)
{
  string rev_string("");
  for(string::const_reverse_iterator rit=in_string.rbegin(); rit<in_string.rend(); rit++)
  {
    rev_string+=*rit;
  }
  return rev_string;
}

inline string repeat_char(const char char_to_append, const uint32_t num_times)
{
  string s;
  // reserve memory to be a little more efficient
  s.reserve(num_times);
  for (uint32_t i=0; i< num_times; i++)
  {
    s += char_to_append;
  }
  return s;
}

inline string substitute(const string& in_s, const string& replace_this, const string& with_this)
{
  string s = in_s;
  size_t pos = s.find(replace_this);
  while (pos != string::npos)
  {
    s.replace(pos, replace_this.size(), with_this);
    pos += with_this.size();
    pos = s.find(replace_this, pos);
  }
  return s;
}

//! Applies regex_m to a vector of strings, returns a vector of strings
//! that do/don't (dependent on bool match) contain an exp.
//! true regex expressions ie !^[\s](\w)+ do not work.
inline vector<string> grep(bool match, string exp, vector<string> lines)
{
  typedef vector<string> Lines;
  Lines matching_lines;
  
  for(Lines::iterator line = lines.begin();
      line != lines.end(); line ++)
    if(match == regex_m(exp, (*line)))
      matching_lines.push_back((*line));
  return matching_lines;    
}
///! Returns first element, then removes it from container.
template <typename T> inline T shift(vector<T>& input)
{
  assert(!input.empty());
  class vector<T>::iterator first = input.begin();
  T retval = (*first);
  input.erase(first);
  
  return retval;
}

// Return the path of the file without the trailing forward-slash
inline string dirname(string file_name)
{
  size_t found = file_name.rfind("/");
  return ((found != string::npos) ? file_name.substr(0, found) : "");
}

inline string create_path(string path)
{
  int status = umask(0);
  
  if (path.find("/") != string::npos) {
    string accumulate_path;
    vector<string> directories = split(path, "/");
    
    for (vector<string>::iterator itr = directories.begin();
         itr != directories.end(); itr ++) {
      accumulate_path.append(*itr);
      status = mkdir(accumulate_path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      accumulate_path.append("/");
    }
  } else {
    status = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  }
  
  if (status && (errno != EEXIST))
  {
    cerr << "Could not create path: '" << path << "'" << endl;
    exit(-1);
  }
  
  return path;
}

// Must be included last because it relies on these functions
#include "anyoption.h"

#endif
