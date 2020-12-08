// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2018)
//   René Fritze     (2014 - 2016, 2018 - 2019)
//   Tim Keil        (2018)
//   Tobias Leibner  (2014, 2018 - 2020)

#ifndef DUNE_XT_FUNCTIONS_SPE10_MODEL1_HH
#define DUNE_XT_FUNCTIONS_SPE10_MODEL1_HH

#include <dune/xt/common/color.hh>
#include <dune/xt/common/configuration.hh>
#if HAVE_DUNE_XT_DATA
#  include <dune/xt/data/paths.hh>
#endif

#include "../checkerboard.hh"

namespace Dune::XT::Functions::Spe10 {
namespace internal {


static constexpr size_t model1_x_elements = 100;
static constexpr size_t model1_y_elements = 1;
static constexpr size_t model1_z_elements = 20;
static constexpr double model_1_length_x = 762.0;
static constexpr double model_1_length_y = 7.62;
static constexpr double model_1_length_z = 152.4;
static constexpr double model1_min_value = 0.001;
static constexpr double model1_max_value = 998.915;


template <class E, size_t r, size_t rC, class R>
class Model1Base : public CheckerboardFunction<E, r, rC, R>
{
  using BaseType = CheckerboardFunction<E, r, rC, R>;

public:
  using typename BaseType::ElementType;
  using typename BaseType::LocalFunctionType;

  using RangeType = typename LocalFunctionType::RangeType;
  using DomainType = typename LocalFunctionType::DomainType;
  using RangeFieldType = typename LocalFunctionType::RangeFieldType;
  using DomainFieldType = typename LocalFunctionType::DomainFieldType;

  static constexpr bool available = true;

  static std::string static_id()
  {
    return GridFunctionInterface<E, 1, 1, R>::static_id() + ".spe10.model1";
  } // ... static_id(...)

private:
  static std::shared_ptr<std::vector<RangeType>> read_values_from_file(const std::string& filename,
                                                                       const RangeFieldType& min,
                                                                       const RangeFieldType& max,
                                                                       const RangeType& unit_range)

  {
    if (!(max > min))
      DUNE_THROW(Dune::RangeError, "max (is " << max << ") has to be larger than min (is " << min << ")!");
    const RangeFieldType scale = (max - min) / (internal::model1_max_value - internal::model1_min_value);
    const RangeFieldType shift = min - scale * internal::model1_min_value;
    // read all the data from the file
    std::ifstream datafile(filename);
    if (!datafile.is_open())
      DUNE_THROW(Exceptions::spe10_data_file_missing, "could not open '" << filename << "'!");
    static constexpr size_t entriesPerDim = model1_x_elements * model1_y_elements * model1_z_elements;
    // create storage (there should be exactly 6000 values in the file, but we only read the first 2000)
    auto data = std::make_shared<std::vector<RangeType>>(entriesPerDim, unit_range);
    double tmp = 0;
    size_t counter = 0;
    while (datafile >> tmp && counter < entriesPerDim)
      (*data)[counter++] *= (tmp * scale) + shift;
    datafile.close();
    if (counter != entriesPerDim)
      DUNE_THROW(Dune::IOError,
                 "wrong number of entries in '" << filename << "' (are " << counter << ", should be " << entriesPerDim
                                                << ")!");
    return data;
  } // ... read_values_from_file(...)

public:
  static Common::Configuration defaults()
  {
    Common::Configuration config;
    config["filename"] =
#if HAVE_DUNE_XT_DATA
        XT::Data::spe10_model1_filename();
#else
        "perm_case1.dat";
#endif
    config["lower_left"] = "[0.0 0.0]";
    config["upper_right"] =
        "[" + Common::to_string(internal::model_1_length_x) + " " + Common::to_string(internal::model_1_length_z) + "]";
    config.set("min_value", internal::model1_min_value);
    config.set("max_value", internal::model1_max_value);
    config["name"] = static_id();
    return config;
  } // ... defaults(...)

  Model1Base(const std::string& filename,
             const DomainType& lowerLeft,
             const DomainType& upperRight,
             const RangeFieldType min,
             const RangeFieldType max,
             const std::string& nm,
             const RangeType& unit_range)
    : BaseType(lowerLeft,
               upperRight,
               {model1_x_elements, model1_z_elements},
               read_values_from_file(filename, min, max, unit_range),
               nm)
  {}
}; // class Model1Base


} // namespace internal


// default, to allow for specialization
template <class E, size_t r, size_t rC = 1, class R = double>
class Model1Function : public GridFunctionInterface<E, r, rC, R>
{
  Model1Function()
  {
    static_assert(AlwaysFalse<E>::value, "Not available for these dimensions!");
  }
};


/**
 * We read only the Kx values from file and scale the unit matrix atm.
 */
template <class E, size_t r, class R>
class Model1Function<E, r, r, R> : public internal::Model1Base<E, r, r, R>
{
  using BaseType = internal::Model1Base<E, r, r, R>;
  using ThisType = Model1Function;

public:
  using DomainFieldType = typename BaseType::DomainFieldType;
  static constexpr size_t domain_dim = E::dimension;
  using RangeFieldType = typename BaseType::RangeFieldType;
  using RangeType = typename BaseType::RangeType;
  using BaseType::defaults;

  Model1Function(const std::string& filename,
                 const Common::FieldVector<DomainFieldType, domain_dim>& lower_left,
                 const Common::FieldVector<DomainFieldType, domain_dim>& upper_right,
                 const RangeFieldType min = internal::model1_min_value,
                 const RangeFieldType max = internal::model1_max_value,
                 const std::string& nm = BaseType::static_id())
    : BaseType(filename, lower_left, upper_right, min, max, nm, unit_matrix())
  {}

private:
  template <size_t d, bool anything = true>
  struct Call
  {
    static RangeType unit_matrix()
    {
      RangeType ret(0.0);
      for (size_t rr = 0; rr < r; ++rr)
        ret[rr][rr] = 1.0;
      return ret;
    }
  };

  template <bool anything>
  struct Call<1, anything>
  {
    static RangeType unit_matrix()
    {
      return RangeType(1);
    }
  };

  static RangeType unit_matrix()
  {
    return Call<r>::unit_matrix();
  }
}; // class Model1< ..., 2, ..., r, r >


} // namespace Dune::XT::Functions::Spe10


#endif // DUNE_XT_FUNCTIONS_SPE10_MODEL1_HH
