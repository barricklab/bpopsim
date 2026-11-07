#include "option_parser.h"

#include <sys/ioctl.h>
#include <unistd.h>

#include <algorithm>
#include <iostream>
#include <stdexcept>

namespace bpopsim {

  namespace {

    constexpr int kLeftColumnWidth = 35;
    constexpr int kDefaultTerminalWidth = 80;
    constexpr int kMinTerminalWidth = kLeftColumnWidth + 20;

    int TerminalWidth()
    {
      if (!isatty(STDOUT_FILENO)) return kDefaultTerminalWidth;
      struct winsize ws = {};
      if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &ws) != 0 || ws.ws_col == 0) return kDefaultTerminalWidth;
      return std::max<int>(ws.ws_col, kMinTerminalWidth);
    }

    // Greedy word wrap; a single word longer than width stays on its own line.
    std::vector<std::string> Wrap(const std::string& text, int width)
    {
      std::vector<std::string> lines;
      std::string line;
      std::istringstream words(text);
      std::string word;
      while (words >> word) {
        if (line.empty()) {
          line = word;
        } else if (static_cast<int>(line.size() + 1 + word.size()) <= width) {
          line += ' ';
          line += word;
        } else {
          lines.push_back(line);
          line = word;
        }
      }
      if (!line.empty() || lines.empty()) lines.push_back(line);
      return lines;
    }

  }  // namespace

  OptionParser::OptionParser(std::string usage_header)
  : usage_header_(std::move(usage_header))
  {}

  OptionParser& OptionParser::operator()(const std::string& names, const std::string& description, NoArgumentTag)
  {
    return Add(names, description, false, false, "");
  }

  OptionParser& OptionParser::operator()(const std::string& names, const std::string& description)
  {
    return Add(names, description, true, false, "");
  }

  OptionParser& OptionParser::Add(const std::string& names, const std::string& description, bool takes_value,
                                  bool has_default, const std::string& default_value)
  {
    Option option;
    size_t comma = names.find(',');
    option.long_name = names.substr(0, comma);
    if (comma != std::string::npos) {
      std::string short_part = names.substr(comma + 1);
      if (short_part.size() != 1) {
        throw std::logic_error("OptionParser: short name must be one character in \"" + names + "\"");
      }
      option.short_name = short_part[0];
    }
    if (option.long_name.empty()) throw std::logic_error("OptionParser: empty option name");
    if (FindLong(option.long_name) != nullptr) {
      throw std::logic_error("OptionParser: duplicate option --" + option.long_name);
    }
    if (option.short_name != 0 && FindShort(option.short_name) != nullptr) {
      throw std::logic_error(std::string("OptionParser: duplicate option -") + option.short_name);
    }
    option.description = description;
    option.takes_value = takes_value;
    option.has_default = has_default;
    option.default_value = default_value;
    options_.push_back(option);
    usage_.push_back({static_cast<int>(options_.size()) - 1, ""});
    return *this;
  }

  void OptionParser::AddUsage(const std::string& line)
  {
    usage_.push_back({-1, line});
  }

  OptionParser::Option* OptionParser::FindLong(const std::string& long_name)
  {
    for (Option& o : options_) {
      if (o.long_name == long_name) return &o;
    }
    return nullptr;
  }

  OptionParser::Option* OptionParser::FindShort(char short_name)
  {
    for (Option& o : options_) {
      if (o.short_name == short_name) return &o;
    }
    return nullptr;
  }

  const OptionParser::Option& OptionParser::Find(const std::string& long_name) const
  {
    for (const Option& o : options_) {
      if (o.long_name == long_name) return o;
    }
    throw std::logic_error("OptionParser: no option named --" + long_name);
  }

  void OptionParser::Parse(int argc, char* const* argv)
  {
    for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];

      if (arg.size() > 2 && arg.compare(0, 2, "--") == 0) {
        // --long or --long=value
        size_t equals = arg.find('=');
        std::string name = arg.substr(2, equals == std::string::npos ? std::string::npos : equals - 2);
        Option* option = FindLong(name);
        if (option == nullptr) throw std::invalid_argument("Unknown option: --" + name);
        option->given = true;
        if (!option->takes_value) {
          if (equals != std::string::npos) throw std::invalid_argument("Option --" + name + " does not take a value");
          continue;
        }
        if (equals != std::string::npos) {
          option->values.push_back(arg.substr(equals + 1));
        } else if (i + 1 < argc) {
          option->values.emplace_back(argv[++i]);
        } else {
          throw std::invalid_argument("Option --" + name + " requires a value");
        }
        continue;
      }

      if (arg.size() > 1 && arg[0] == '-') {
        // -s value, -svalue, or combined flags -abc (an option char ends the run
        // and takes the rest of the argument, or the next argument, as its value)
        for (size_t j = 1; j < arg.size(); ++j) {
          Option* option = FindShort(arg[j]);
          if (option == nullptr) throw std::invalid_argument(std::string("Unknown option: -") + arg[j]);
          option->given = true;
          if (!option->takes_value) continue;
          if (j + 1 < arg.size()) {
            option->values.push_back(arg.substr(j + 1));
          } else if (i + 1 < argc) {
            option->values.emplace_back(argv[++i]);
          } else {
            throw std::invalid_argument(std::string("Option -") + arg[j] + " requires a value");
          }
          break;
        }
        continue;
      }

      throw std::invalid_argument("Unexpected argument: " + arg);
    }
  }

  bool OptionParser::count(const std::string& long_name) const
  {
    return Find(long_name).given;
  }

  std::string OptionParser::operator[](const std::string& long_name) const
  {
    const Option& option = Find(long_name);
    if (!option.takes_value) return option.given ? "true" : "";
    if (!option.values.empty()) {
      std::string joined = option.values[0];
      for (size_t i = 1; i < option.values.size(); ++i) {
        joined += '\n';
        joined += option.values[i];
      }
      return joined;
    }
    return option.has_default ? option.default_value : "";
  }

  void OptionParser::PrintUsage() const
  {
    PrintUsage(std::cout);
  }

  void OptionParser::PrintUsage(std::ostream& out) const
  {
    const int width = TerminalWidth();
    const int right_width = width - kLeftColumnWidth;
    const std::string indent(kLeftColumnWidth, ' ');

    out << '\n' << usage_header_ << '\n';
    out << "Allowed Options\n";

    for (const UsageEntry& entry : usage_) {
      if (entry.option_index < 0) {
        out << entry.text << '\n';
        continue;
      }
      const Option& o = options_[entry.option_index];

      std::string left = "  ";
      if (o.short_name != 0) left += std::string("-") + o.short_name + ",";
      left += "--" + o.long_name;
      if (o.takes_value) left += " <arg>";

      std::string description = o.description;
      if (o.has_default) description += " (DEFAULT=" + o.default_value + ")";
      std::vector<std::string> lines = Wrap(description, right_width);

      if (static_cast<int>(left.size()) < kLeftColumnWidth) {
        left += std::string(kLeftColumnWidth - left.size(), ' ');
        out << left << lines[0] << '\n';
      } else {
        // Name too long for the left column: put the description on the next line.
        out << left << '\n' << indent << lines[0] << '\n';
      }
      for (size_t i = 1; i < lines.size(); ++i) out << indent << lines[i] << '\n';
    }
    out << '\n';
  }

}  // namespace bpopsim
