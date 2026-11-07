#pragma once

#include <iosfwd>
#include <sstream>
#include <string>
#include <vector>

namespace bpopsim {

  //! Tag for registering an option that takes no value (a flag).
  struct NoArgumentTag
  {
  };
  constexpr NoArgumentTag NoArgument{};

  //! Minimal GNU/POSIX style command-line parser.
  //!
  //! Options are registered as "long,s" (long name plus one-character short
  //! name) or just "long". Accepted spellings on the command line:
  //!   --long value   --long=value   -s value   -svalue   -abc (combined flags)
  //! An option given several times accumulates its values, joined with "\n",
  //! which is how from_string<vector<T>>() reads them. Flags report "true".
  //! Unknown options, a missing value, or a stray positional argument throw
  //! std::invalid_argument; the caller decides how to report it.
  class OptionParser
  {
  public:
    explicit OptionParser(std::string usage_header = "");

    //! Register a flag.
    OptionParser& operator()(const std::string& names, const std::string& description, NoArgumentTag);
    //! Register an option that takes a value and has no default.
    OptionParser& operator()(const std::string& names, const std::string& description);
    //! Register an option that takes a value, with a default shown in the usage text.
    template <class T>
    OptionParser& operator()(const std::string& names, const std::string& description, const T& default_value)
    {
      std::ostringstream ss;
      ss << default_value;
      return Add(names, description, true, true, ss.str());
    }

    //! Append a line of free text to the usage output, in registration order
    //! (section headings between options, or error messages after parsing).
    void AddUsage(const std::string& line);

    //! Parse argv (argv[0] is skipped). Throws std::invalid_argument on error.
    void Parse(int argc, char* const* argv);

    //! True if the option or flag appeared on the command line.
    bool count(const std::string& long_name) const;

    //! The option's value(s) joined with "\n"; its default if not given; "true"
    //! for a flag that was given; "" otherwise.
    std::string operator[](const std::string& long_name) const;

    void PrintUsage(std::ostream& out) const;
    void PrintUsage() const;

  private:
    struct Option
    {
      std::string long_name;
      char short_name = 0;
      std::string description;
      bool takes_value = false;
      bool has_default = false;
      std::string default_value;
      bool given = false;
      std::vector<std::string> values;
    };
    // Usage output is a sequence of option entries and free-text lines.
    struct UsageEntry
    {
      int option_index;  // -1 for a text line
      std::string text;
    };

    OptionParser& Add(const std::string& names, const std::string& description, bool takes_value, bool has_default,
                      const std::string& default_value);
    const Option& Find(const std::string& long_name) const;
    Option* FindLong(const std::string& long_name);
    Option* FindShort(char short_name);

    std::string usage_header_;
    std::vector<Option> options_;
    std::vector<UsageEntry> usage_;
  };

}  // namespace bpopsim
