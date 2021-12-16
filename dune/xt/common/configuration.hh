// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012 - 2017)
//   René Fritze     (2012 - 2020)
//   Sven Kaulmann   (2013)
//   Tim Keil        (2018)
//   Tobias Leibner  (2014, 2016, 2018 - 2020)

#ifndef DUNE_XT_COMMON_CONFIGURATION_HH
#define DUNE_XT_COMMON_CONFIGURATION_HH

#include <iosfwd>
#include <type_traits>

#include <boost/lexical_cast/bad_lexical_cast.hpp>

#include <dune/common/visibility.hh>
#include <dune/common/parametertree.hh>

#include <dune/xt/common/algorithm.hh>
#include <dune/xt/common/validation.hh>

namespace Dune {
namespace XT::Common {


struct ConfigurationDefaults
{
  ConfigurationDefaults(bool warn_on_default_access_in = false,
                        bool log_on_exit_in = false,
                        std::string logfile_in = std::string("data/log/dxtc_parameter.log"));
  const bool warn_on_default_access;
  const bool log_on_exit;
  const std::string logfile;
};

namespace internal {
template <class T>
struct Typer
{
  using type = typename std::conditional<std::is_same<T, const char*>::value, std::string, T>::type;
};

} // namespace internal

class Configuration : public Dune::ParameterTree
{
  using BaseType = Dune::ParameterTree;

public:
  Configuration();

  Configuration(const ParameterTree& tree, const ConfigurationDefaults& defaults = ConfigurationDefaults());

  Configuration(const ParameterTree& tree_in, const std::string& sub_id);

  Configuration(const Configuration& other);

  /**
   * \brief Allows to construct from a list of key/value pairs.
   *
   *        This ctor is not marked explicit on purpose to allow for code such as
\code
Configuration default_config()
{
  return {{"type", "custom"},
          {"tolerance", "1e-10"}};
}
\endcode
   *        or:
\code
some_function_which_expects_a_config({{"type", "custom"}, {"tolerance", "1e-10"}});
\endcode
   */
  Configuration(const std::initializer_list<std::pair<std::string, std::string>>& key_value_pairs);

  explicit Configuration(std::istream& in, ConfigurationDefaults defaults = ConfigurationDefaults());

  Configuration(const std::string& in, ConfigurationDefaults defaults = ConfigurationDefaults());

  //! read ParameterTree from given arguments and call Configuration(const ParameterTree& tree)
  Configuration(int argc, char** argv, ConfigurationDefaults defaults = ConfigurationDefaults());

  template <class T>
  Configuration(const std::vector<std::string>& keys,
                const std::initializer_list<T> values_in,
                const ConfigurationDefaults& defaults = ConfigurationDefaults())
    : warn_on_default_access_(defaults.warn_on_default_access)
    , log_on_exit_(defaults.log_on_exit)
    , logfile_(defaults.logfile)
  {
    const auto values = make_string_sequence(values_in.begin(), values_in.end());
    if (keys.size() != values.size()) {

      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The size of 'keys' (" << keys.size() << ") does not match the size of 'values' (" << values_in.size()
                                        << ")!");
    }
    size_t ii = 0;
    for (auto value : values)
      set(keys[ii++], value);
    setup_();
  }

  ~Configuration();

  Configuration& operator=(const Configuration& other);

  /**
   * \defgroup base ´´These methods replace or override those from Dune::ParameterTree.``
   * \{
   */

  //! check if key is existing in tree_
  bool has_key(const std::string& key) const;

  //! check if sub is existing in tree_
  bool has_sub(const std::string& subTreeName) const;

  /** \brief print the ParameterTree
   *  \param out output stream
   *  \param prefix to be prepended to each line
   */
  void report(std::ostream& out = std::cout, const std::string& prefix = "") const;

  /**
   * @attention Please note the difference to Dune::ParameterTree::sub (return: value vs. reference)!
   */
  Configuration sub(const std::string& sub_id,
                    bool fail_if_missing = true,
                    const Configuration& default_value = Configuration()) const;

  /**
   * \}
   */

  //! const get without default value, with validation
  template <class T, class Validator = ValidateAny<typename internal::Typer<T>::type>>
  typename internal::Typer<T>::type
  get(std::string key,
      const size_t size,
      const size_t cols = 0,
      const ValidatorInterface<T, Validator>& validator = ValidateAny<typename internal::Typer<T>::type>()) const
  {
    if (!has_key(key))
      DUNE_THROW(Exceptions::configuration_error,
                 "This Configuration (see below) does not contain the key '"
                     << key << "' and there was no default value provided!\n======================\n"
                     << report_string());
    return get_valid_value(key, T(), validator, size, cols);
  } // ... get(...)

  //! const get without default value, with validation
  template <class T, class Validator = ValidateAny<typename internal::Typer<T>::type>>
  typename internal::Typer<T>::type
  get(std::string key,
      const ValidatorInterface<typename internal::Typer<T>::type, Validator>& validator =
          ValidateAny<typename internal::Typer<T>::type>()) const
  {
    if (!has_key(key))
      DUNE_THROW(Exceptions::configuration_error,
                 "This Configuration (see below) does not contain the key '"
                     << key << "' and there was no default value provided!\n======================\n"
                     << report_string());
    return get_valid_value(key, T(), validator, 0, 0);
  } // ... get(...)

  //! get variation with default value, validation
  template <typename T, class Validator = ValidateAny<typename internal::Typer<T>::type>>
  typename internal::Typer<T>::type
  get(std::string key,
      T def,
      const size_t size,
      const size_t cols = 0,
      const ValidatorInterface<typename internal::Typer<T>::type, Validator>& validator =
          ValidateAny<typename internal::Typer<T>::type>()) const
  {
    return get_(key, def, validator, size, cols);
  } // ... get(...)

  //! get variation with default value, validation
  template <typename T, class Validator = ValidateAny<typename internal::Typer<T>::type>>
  typename internal::Typer<T>::type
  get(std::string key,
      T def,
      const ValidatorInterface<typename internal::Typer<T>::type, Validator>& validator =
          ValidateAny<typename internal::Typer<T>::type>()) const
  {
    return get_(key, def, validator, 0, 0);
  } // ... get(...)

  /**
   * \defgroup set ´´These methods allow to set key: value pairs.``
   * \{
   */

  //! set value to key in Configuration
  template <class T>
  void set(const std::string& key, const T& value, const bool overwrite = false)
  {
    if (has_key(key) && !overwrite)
      DUNE_THROW(Exceptions::configuration_error,
                 "While setting '" << key << "' in this configuration (see below), it already exists and you requested "
                                   << "no overwrite!\n======================\n"
                                   << report_string());
    BaseType::operator[](key) = to_string(value);
  } // ... set(..., T, ...)

  void set(const std::string& key, const char* value, const bool overwrite = false);

  /**
   * \}
   */

  /**
   * \defgroup add ´´These methods allow to merge trees.``
   * \{
   */

  /** \brief add another Configuration to this (merge tree_s and requests_map_s)
   *  \param other Configuration to add
   *  \param sub_id if not empty, other.tree_ is merged in as a sub "sub_id" of tree_
   *  \param overwrite if true, existing values are overwritten by other's values to the same key
   */
  Configuration& add(const Configuration& other, const std::string& sub_id = "", const bool overwrite = false);

  /** \brief add a Dune::ParameterTree paramtree to this (merge tree_ and paramtree)
   *  \param paramtree ParameterTree to add
   *  \param sub_id if not empty, paramtree is merged in as a sub "sub_id" of tree_
   *  \param overwrite if true, existing values are overwritten by paramtree's values to the same key
   */
  Configuration& add(const ParameterTree& other, const std::string& sub_id = "", const bool overwrite = false);

  //! add another Configuration to this (merge tree_s and requests_map_s)
  Configuration& operator+=(const Configuration& other);

  //! add this and another Configuration (merge tree_s and requests_map_s)
  Configuration operator+(const Configuration& other) const;

  /**
   * \}
   */
  void set_warn_on_default_access(const bool value);
  void set_log_on_exit(const bool value);
  void set_logfile(const std::string& logfile);

  //! check if tree_ is empty
  bool empty() const;

  //! store output of report(..., prefix) in std::string
  std::string report_string(const std::string& prefix = "") const;

  std::map<std::string, std::string> flatten() const;

  /** get parameters from parameter file or key-value pairs given on the command line and store in Configuration (and
  load into fem parameter, if available) */
  void read_command_line(int argc, char** argv);

  //! search command line options for key-value pairs and add them to Configuration
  void read_options(int argc, char** argv);

private:
  void setup_();

  void add_tree_(const Configuration& other, const std::string& sub_id, const bool overwrite);

  //! get value from tree and validate with validator
  template <typename T, class Validator>
  T get_valid_value(const std::string& key,
                    T def,
                    const ValidatorInterface<typename internal::Typer<T>::type, Validator>& validator,
                    const size_t size,
                    const size_t cols) const
  {
    std::string valstring = BaseType::get(key, to_string(def));
    T val;
    try {
      val = from_string<T>(valstring, size, cols);
    } catch (boost::bad_lexical_cast& e) {
      DUNE_THROW(Exceptions::external_error,
                 "Error in boost while converting the string '"
                     << valstring << "' to type '" << Typename<T>::value() << "':\n"
                     << e.what() << "\non accessing key " << key << " with default " << to_string(def));
    } catch (std::exception& e) {
      DUNE_THROW(Exceptions::external_error,
                 "Error in the stl while converting the string '"
                     << valstring << "' to type '" << Typename<T>::value() << "':\n"
                     << e.what() << "\non accessing key " << key << " with default " << to_string(def));
    }
    if (validator(val))
      return val;
    DUNE_THROW(Exceptions::configuration_error, validator.msg(val));
  } // ... get_valid_value(...)

  /** \brief all public get signatures call this one
   *  \param key requested key
   *  \param def default value
   *  \param validator validator that is used to validate the value before it is returned
   *  \param request Request that is stored in requests_map_
   *  \param size Determines the size of the returning container (size if T is a vector type,
   *  rows if T is a matrix type, 0 means automatic).
   *  \param cols Determines the number of columns of the returning matrix if T is a matrix type
   *  (0 means automatic, ignored, if T is a vector or scalar type).
   *  \param def_provided bool to indicate if the calling get method provides a default value
   *  \return value associated to key in Configuration (interpreted as type T),
   *  def if key does not exist in Configuration
   */
  template <typename T, class Validator>
  typename internal::Typer<T>::type get_(const std::string& key,
                                         typename internal::Typer<T>::type def,
                                         const ValidatorInterface<T, Validator>& validator,
                                         const size_t size,
                                         const size_t cols) const
  {
#ifndef NDEBUG
    if (warn_on_default_access_ && !has_key(key)) {
      std::cerr << color_string("WARNING:", Colors::brown) << " using default value for parameter \"" << key << "\""
                << std::endl;
    }
#endif // ifndef NDEBUG
    return get_valid_value(key, def, validator, size, cols);
  } // ... get_(...)

  //! read Dune::ParameterTree from file
  static ParameterTree initialize(const std::string& filename);

  //! read Dune::ParameterTree from istream
  static ParameterTree initialize(std::istream& in);

  //! read Dune::ParameterTree from arguments
  static ParameterTree initialize(int argc, char** argv);

  //! read Dune::ParameterTree from arguments and file
  static ParameterTree initialize(int argc, char** argv, std::string filename);

  void report_as_sub(std::ostream& out, const std::string& prefix, const std::string& sub_path) const;

  std::string find_common_prefix(const BaseType& subtree, const std::string& previous_prefix) const;

  void report_flatly(const BaseType& subtree, const std::string& prefix, std::ostream& out) const;

  bool warn_on_default_access_;
  bool log_on_exit_;
  std::string logfile_;
}; // class Configuration

std::ostream& operator<<(std::ostream& out, const Configuration& config);

bool operator==(const Configuration& left, const Configuration& right);

bool operator!=(const Configuration& left, const Configuration& right);

//! global Configuration instance
DUNE_EXPORT inline Configuration& Config()
{
  static Configuration parameters;
  return parameters;
}

} // namespace XT::Common

bool operator==(const ParameterTree& left, const ParameterTree& right);

bool operator!=(const ParameterTree& left, const ParameterTree& right);

} // namespace Dune
namespace std {

template <>
struct less<Dune::ParameterTree>
{
  using result_type = bool;
  using first_argument_type = Dune::ParameterTree;
  using second_argument_type = Dune::ParameterTree;

  bool operator()(const Dune::ParameterTree& lhs, const Dune::ParameterTree& rhs) const;
}; // struct less< ParameterTree >

template <>
struct less<Dune::XT::Common::Configuration>
{
  using result_type = bool;
  using first_argument_type = Dune::XT::Common::Configuration;
  using second_argument_type = Dune::XT::Common::Configuration;

  bool operator()(const Dune::XT::Common::Configuration& lhs, const Dune::XT::Common::Configuration& rhs) const;
}; // struct less< ParameterTree >

} // namespace std

#define DXTC_CONFIG Dune::XT::Common::Config()

template <typename T>
static auto DXTC_CONFIG_GET(std::string key, T def) -> decltype(DXTC_CONFIG.get(key, def))
{
  return DXTC_CONFIG.get(key, def);
}

template <typename T, class V>
static auto
DXTC_CONFIG_GETV(std::string key,
                 T def,
                 const Dune::XT::Common::ValidatorInterface<typename Dune::XT::Common::internal::Typer<T>::type, V>& v)
    -> decltype(DXTC_CONFIG.get(key, def, v))
{
  return DXTC_CONFIG.get(key, def, v);
}

#endif // DUNE_XT_COMMON_CONFIGURATION_HH
