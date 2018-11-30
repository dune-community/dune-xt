// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2018)
//   René Fritze     (2015 - 2016, 2018)
//   TiKeil          (2018)
//   Tobias Leibner  (2017 - 2018)

#ifndef DUNE_XT_FUNCTIONS_SPE10_MODEL2_HH
#define DUNE_XT_FUNCTIONS_SPE10_MODEL2_HH

#include <dune/xt/common/color.hh>
#include <dune/xt/common/configuration.hh>
#include <dune/xt/data/paths.hh>

#include "../checkerboard.hh"

namespace Dune {
namespace XT {
namespace Functions {
namespace Spe10 {
namespace internal {


static const size_t model2_x_elements = 60;
static const size_t model2_y_elements = 220;
static const size_t model2_z_elements = 85;
static const double model_2_length_x = 1;
static const double model_2_length_y = 3.667;
static const double model_2_length_z = 1.417;

} // namespace internal


// default, to allow for specialization
template <class E, size_t r, size_t rC = 1, class R = double>
class Model2Function : public GridFunctionInterface<E, r, rC, R>
{
  Model2Function()
  {
    static_assert(AlwaysFalse<E>::value, "Not available for these dimensions!");
  }
};


template <class E, class R>
class Model2Function<E, 3, 3, R> : public CheckerboardFunction<E, 3, 3, R>
{
  using BaseType = CheckerboardFunction<E, 3, 3, R>;
  using ThisType = Model2Function<E, 3, 3, R>;

public:
  using typename BaseType::ElementType;
  using typename BaseType::LocalFunctionType;

  using RangeType = typename LocalFunctionType::RangeType;
  using DomainType = typename LocalFunctionType::DomainType;
  using RangeFieldType = typename LocalFunctionType::RangeFieldType;
  using DomainFieldType = typename LocalFunctionType::DomainFieldType;
  static const constexpr size_t domain_dim = E::dimension;

  static const bool available = true;

  static std::string static_id()
  {
    return GridFunctionInterface<E, 3, 3, R>::static_id() + ".spe10.model1";
  } // ... static_id(...)

private:
  static std::vector<RangeType> read_values_from_file(const std::string& filename,
                                                      const Common::FieldVector<size_t, domain_dim>& number_of_elements)

  {
    // read all the data from the file
    std::ifstream file(filename);
    if (!file.is_open()) {
      DXTC_LOG_ERROR_0 << "The SPE10-permeability data file could not be opened. This file does\n"
                       << "not come with the dune-multiscale repository due to file size. To download it\n"
                       << "execute\n"
                       << "wget http://www.spe.org/web/csp/datasets/por_perm_case2a.zip\n"
                       << "unzip the file and move the file 'spe_perm.dat' to\n"
                       << "dune-multiscale/dune/multiscale/problems/spe10_permeability.dat!\n";
      DUNE_THROW(IOError, "Data file for Groundwaterflow permeability could not be opened!");
    }
    const size_t entries_per_coordinate =
        number_of_elements[0] /*x*/ * number_of_elements[1] /*y*/ * number_of_elements[2] /*z*/;

    std::vector<double> values_in_file(3 * entries_per_coordinate);
    double val;
    file >> val;
    size_t filecounter = 0;
    while (!file.eof() && filecounter < 3 * entries_per_coordinate) {
      values_in_file[filecounter++] = val;
      file >> val; // sets EOF flag if no value found
    }
    file.close();
    /* todo */
    // if (filecounter != 3366000)
    //#warning you are not using the entire data file. Use default number_of_elements instead.

    std::vector<RangeType> data(entries_per_coordinate, RangeType(0));

    for (size_t ii = 0; ii < entries_per_coordinate; ++ii) {
      for (size_t dim = 0; dim < domain_dim; ++dim) {
        const auto idx = ii + dim * entries_per_coordinate;
        data[ii][dim][dim] = values_in_file[idx];
      }
    }
    return data;
  }


public:
  static Common::Configuration defaults()
  {
    Common::Configuration config;
    config["type"] = static_id();
    config["filename"] = Data::spe10_model2_filename();
    config["lower_left"] = "[0.0 0.0 0.0]";
    config["upper_right"] = "[" + Common::to_string(internal::model_2_length_x) + " "
                            + Common::to_string(internal::model_2_length_y) + " "
                            + Common::to_string(internal::model_2_length_z) + "]";
    config["number_of_elements"] = "[" + Common::to_string(internal::model2_x_elements) + " "
                                   + Common::to_string(internal::model2_y_elements) + " "
                                   + Common::to_string(internal::model2_z_elements) + "]";
    config["name"] = static_id();
    return config;
  } // ... defaults(...)

  Model2Function(const std::string& filename,
                 const Common::FieldVector<DomainFieldType, domain_dim>& lower_left,
                 const Common::FieldVector<DomainFieldType, domain_dim>& upper_right,
                 const Common::FieldVector<size_t, domain_dim>& number_of_elements = {internal::model2_x_elements,
                                                                                      internal::model2_y_elements,
                                                                                      internal::model2_z_elements},
                 const std::string nm = BaseType::static_id())
    : BaseType(lower_left, upper_right, number_of_elements, read_values_from_file(filename, number_of_elements), nm)

  {
  }
}; // class Model2Function


} // namespace Spe10
} // namespace Functions
} // namespace XT
} // namespace Dune


#endif // DUNE_XT_FUNCTIONS_SPE10_MODEL2_HH
