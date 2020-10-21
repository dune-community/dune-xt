// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2018)
//   René Fritze     (2013 - 2019)
//   Sven Kaulmann   (2013)
//   Tobias Leibner  (2014, 2018 - 2020)

#include "config.h"

#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include <dune/common/parametertreeparser.hh>

#include <dune/xt/common/filesystem.hh>

#include "configuration.hh"

namespace Dune {
namespace XT::Common {


ConfigurationDefaults::ConfigurationDefaults(bool warn_on_default_access_in,
                                             bool log_on_exit_in,
                                             std::string logfile_in)
  : warn_on_default_access(warn_on_default_access_in)
  , log_on_exit(log_on_exit_in)
  , logfile(logfile_in)
{}

Configuration::Configuration()
  : BaseType()
  , warn_on_default_access_(ConfigurationDefaults().warn_on_default_access)
  , log_on_exit_(ConfigurationDefaults().log_on_exit)
  , logfile_(ConfigurationDefaults().logfile)
{
  setup_();
}

Configuration::Configuration(const Dune::ParameterTree& tree_in, ConfigurationDefaults defaults)
  : BaseType(tree_in)
  , warn_on_default_access_(defaults.warn_on_default_access)
  , log_on_exit_(defaults.log_on_exit)
  , logfile_(defaults.logfile)
{
  setup_();
}

Configuration::Configuration(const ParameterTree& tree_in, const std::string sub_id)
  : BaseType()
  , warn_on_default_access_(ConfigurationDefaults().warn_on_default_access)
  , log_on_exit_(ConfigurationDefaults().log_on_exit)
  , logfile_(ConfigurationDefaults().logfile)
{
  setup_();
  add(tree_in, sub_id);
}

Configuration::Configuration(const Configuration& other)
  : BaseType(other)
  , warn_on_default_access_(other.warn_on_default_access_)
  , log_on_exit_(other.log_on_exit_)
  , logfile_(other.logfile_)
{}

Configuration::Configuration(const std::initializer_list<std::pair<std::string, std::string>>& key_value_pairs)
  : BaseType()
  , warn_on_default_access_(ConfigurationDefaults().warn_on_default_access)
  , log_on_exit_(ConfigurationDefaults().log_on_exit)
  , logfile_(ConfigurationDefaults().logfile)
{
  for (const auto& key_value_pair : key_value_pairs) {
    const auto& key = key_value_pair.first;
    const auto& value = key_value_pair.second;
    set(key, value);
  }
}

Configuration::Configuration(std::istream& in, ConfigurationDefaults defaults)
  : Configuration(initialize(in), defaults)
{}

Configuration::Configuration(const std::string& in, ConfigurationDefaults defaults)
  : Configuration(initialize(in), defaults)
{}

Configuration::Configuration(int argc, char** argv, ConfigurationDefaults defaults)
  : Configuration::Configuration(initialize(argc, argv), defaults)
{}

Configuration::~Configuration()
{
  if (log_on_exit_ && !empty()) {
    test_create_directory(directory_only(logfile_));
    report(*make_ofstream(logfile_));
  }
}

void Configuration::set_warn_on_default_access(const bool value)
{
  warn_on_default_access_ = value;
}

void Configuration::set_log_on_exit(const bool value)
{
  if (!log_on_exit_ && value)
    test_create_directory(directory_only(logfile_));
  log_on_exit_ = value;
}

void Configuration::set_logfile(const std::string logfile)
{
  if (logfile.empty())
    DUNE_THROW(Exceptions::wrong_input_given, "logfile must not be empty!");
  if (log_on_exit_)
    test_create_directory(directory_only(logfile_));
}

// method definitions for Configuration
bool Configuration::has_key(const std::string& key) const
{
  return BaseType::hasKey(key);
}

Configuration Configuration::sub(const std::string& sub_id, bool fail_if_missing, Configuration default_value) const
{
  if ((empty() || !has_sub(sub_id)) && !fail_if_missing)
    return default_value;
  if (empty()) {
    DUNE_THROW(Exceptions::configuration_error,
               "You can not get anything from an empty Configuration, use has_sub(\"" << sub_id
                                                                                      << "\") to check first!");
  }
  if (sub_id.empty())
    DUNE_THROW(Exceptions::configuration_error, "Given sub_id must not be empty!");
  if (!has_sub(sub_id))
    DUNE_THROW(Exceptions::configuration_error,
               "Subtree '" << sub_id << "' does not exist in this Configuration (see below), use has_sub(\"" << sub_id
                           << "\") to check first!"
                           << "\n======================\n"
                           << report_string());
  return Configuration(BaseType::sub(sub_id));
} // ... sub(...)

bool Configuration::has_sub(const std::string& subTreeName) const
{
  return BaseType::hasSub(subTreeName);
}

void Configuration::set(const std::string& key, const char* value, const bool overwrite)
{
  set(key, std::string(value), overwrite);
}

Configuration& Configuration::add(const Configuration& other, const std::string sub_id, const bool overwrite)
{
  add_tree_(other, sub_id, overwrite);
  return *this;
} // ... add(...)

Configuration& Configuration::add(const ParameterTree& other, const std::string sub_id, const bool overwrite)
{
  add_tree_(Configuration(other), sub_id, overwrite);
  return *this;
} // ... add(...)

Configuration& Configuration::operator+=(const Configuration& other)
{
  add(other);
  return *this;
} // ... operator+=

Configuration Configuration::operator+(const Configuration& other) const
{
  Configuration ret(*this);
  ret += other;
  return ret;
} // ... operator+

Configuration& Configuration::operator=(const Configuration& other)
{
  if (this != &other) {
    BaseType::operator=(other);
    warn_on_default_access_ = other.warn_on_default_access_;
    log_on_exit_ = other.log_on_exit_;
    logfile_ = other.logfile_;
  }
  return *this;
} // ... operator=(...)

bool Configuration::empty() const
{
  return this->getValueKeys().empty() && this->getSubKeys().empty();
}

void Configuration::report(std::ostream& out, const std::string& prefix) const
{
  if (empty())
    return;

  if (getSubKeys().size() == 0) {
    report_as_sub(out, prefix, "");
  } else if (getValueKeys().size() == 0) {
    const std::string common_prefix = find_common_prefix(*this, "");
    if (!common_prefix.empty()) {
      out << prefix << "[" << common_prefix << "]" << std::endl;
      const Configuration& commonSub = sub(common_prefix);
      report_flatly(commonSub, prefix, out);
    } else
      report_as_sub(out, prefix, "");
  } else {
    report_as_sub(out, prefix, "");
  }
} // ... report(...)

std::string Configuration::report_string(const std::string& prefix) const
{
  std::stringstream stream;
  report(stream, prefix);
  return stream.str();
} // ... report_string(...)

void Configuration::read_command_line(int argc, char* argv[])
{
  if (argc < 2) {
    boost::format usage("usage: %s parameter.file *[-section.key override-value]");
    DUNE_THROW(Dune::Exception, (usage % argv[0]).str());
  }
  if (boost::filesystem::exists(argv[1]))
    Dune::ParameterTreeParser::readINITree(argv[1], *this);
  Dune::ParameterTreeParser::readOptions(argc, argv, *this);
  // datadir and logdir may be given from the command line...
  setup_();
} // readCommandLine

void Configuration::read_options(int argc, char* argv[])
{
  Dune::ParameterTreeParser::readOptions(argc, argv, *this);
}

void Configuration::setup_()
{
  if (logfile_.empty())
    logfile_ = boost::filesystem::path(ConfigurationDefaults().logfile).string();
  if (has_key("global.datadir") && has_key("logging.dir"))
    logfile_ = (boost::filesystem::path(get<std::string>("global.datadir")) / get<std::string>("logging.dir")
                / "dxtc_parameter.log")
                   .string();
  logfile_ = boost::filesystem::path(logfile_).string();
} // ... setup_(...)

void Configuration::add_tree_(const Configuration& other, const std::string sub_id, const bool overwrite)
{
  for (const auto& element : other.flatten()) {
    auto key = element.first;
    if (!sub_id.empty())
      key = sub_id + "." + key;
    try {
      set(key, element.second, overwrite);
    } catch (Exceptions::configuration_error& ee) {
      DUNE_THROW(Exceptions::configuration_error,
                 "There was an error adding other (see below) to this:\n\n"
                     << ee.what() << "\n\n"
                     << "======================\n"
                     << other.report_string() << "\n");
    }
  }
} // ... add_tree_(...)

ParameterTree Configuration::initialize(const std::string filename)
{
  ParameterTree param_tree;
  Dune::ParameterTreeParser::readINITree(filename, param_tree);
  return param_tree;
} // ... initialize(...)

ParameterTree Configuration::initialize(std::istream& in)
{
  ParameterTree param_tree;
  Dune::ParameterTreeParser::readINITree(in, param_tree);
  return param_tree;
} // ... initialize(...)

ParameterTree Configuration::initialize(int argc, char** argv)
{
  ParameterTree param_tree;
  if (argc == 2) {
    Dune::ParameterTreeParser::readINITree(argv[1], param_tree);
  } else if (argc > 2) {
    Dune::ParameterTreeParser::readOptions(argc, argv, param_tree);
  }
  if (param_tree.hasKey("paramfile")) {
    Dune::ParameterTreeParser::readINITree(param_tree.get<std::string>("paramfile"), param_tree, false);
  }
  return param_tree;
} // ... initialize(...)

ParameterTree Configuration::initialize(int argc, char** argv, std::string filename)
{
  ParameterTree param_tree;
  if (argc == 1) {
    Dune::ParameterTreeParser::readINITree(filename, param_tree);
  } else if (argc == 2) {
    Dune::ParameterTreeParser::readINITree(argv[1], param_tree);
  } else {
    Dune::ParameterTreeParser::readOptions(argc, argv, param_tree);
  }
  if (param_tree.hasKey("paramfile")) {
    Dune::ParameterTreeParser::readINITree(param_tree.get<std::string>("paramfile"), param_tree, false);
  }
  return param_tree;
} // ... initialize(...)

void Configuration::report_as_sub(std::ostream& out, const std::string& prefix, const std::string& sub_path) const
{
  for (const auto& key : getValueKeys()) {
    out << prefix << key << " = " << ParameterTree::get<std::string>(key) << std::endl;
  }
  for (const auto& subkey : getSubKeys()) {
    Configuration sub_tree(sub(subkey));
    if (sub_tree.getValueKeys().size())
      out << prefix << "[" << sub_path << subkey << "]" << std::endl;
    sub_tree.report_as_sub(out, prefix, sub_path + subkey + ".");
  }
} // ... report_as_sub(...)

std::string Configuration::find_common_prefix(const BaseType& subtree, const std::string previous_prefix) const
{
  const auto& valuekeys = subtree.getValueKeys();
  const auto& subkeys = subtree.getSubKeys();
  if (valuekeys.size() == 0 && subkeys.size() == 1) {
    // we append the subs name
    if (previous_prefix.empty())
      return find_common_prefix(subtree.sub(subkeys[0]), subkeys[0]);
    return find_common_prefix(subtree.sub(subkeys[0]), previous_prefix + "." + subkeys[0]);
  } else {
    // end of the recursion, return the previous prefix
    return previous_prefix;
  }
} // ... find_common_prefix(...)

void Configuration::report_flatly(const BaseType& subtree, const std::string& prefix, std::ostream& out) const
{
  // report all the keys
  for (auto key : subtree.getValueKeys())
    out << prefix << key << " = " << subtree[key] << std::endl;
  // report all the subs
  for (auto subkey : subtree.getSubKeys()) {
    if (prefix.empty())
      report_flatly(subtree.sub(subkey), subkey + ".", out);
    else
      report_flatly(subtree.sub(subkey), prefix + subkey + ".", out);
  }
} // ... report_flatly(...)

std::map<std::string, std::string> Configuration::flatten() const
{
  std::map<std::string, std::string> ret;
  for (const auto& kk : getValueKeys())
    ret[kk] = get<std::string>(kk);
  for (const auto& ss : getSubKeys()) {
    const auto flat_sub = sub(ss).flatten();
    for (const auto& element : flat_sub)
      ret[ss + "." + element.first] = element.second;
  }
  return ret;
} // ... flatten(...)

std::ostream& operator<<(std::ostream& out, const Configuration& config)
{
  config.report(out);
  return out;
}

bool operator==(const Configuration& left, const Configuration& right)
{
  return left.flatten() == right.flatten();
}

bool operator!=(const Configuration& left, const Configuration& right)
{
  return !(left == right);
}

} // namespace XT::Common

bool operator==(const ParameterTree& left, const ParameterTree& right)
{
  return Dune::XT::Common::Configuration(left).flatten() == Dune::XT::Common::Configuration(right).flatten();
}

bool operator!=(const ParameterTree& left, const ParameterTree& right)
{
  return !(left == right);
}


} // namespace Dune
namespace std {

bool less<Dune::ParameterTree>::operator()(const Dune::ParameterTree& lhs, const Dune::ParameterTree& rhs) const
{
  return Dune::XT::Common::Configuration(lhs).flatten() < Dune::XT::Common::Configuration(rhs).flatten();
}

bool less<Dune::XT::Common::Configuration>::operator()(const Dune::XT::Common::Configuration& lhs,
                                                       const Dune::XT::Common::Configuration& rhs) const
{
  return lhs.flatten() < rhs.flatten();
}

} // namespace std
