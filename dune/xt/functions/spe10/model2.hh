// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2015 - 2016)
//   Tobias Leibner  (2017)

#ifndef DUNE_XT_FUNCTIONS_SPE10_MODEL2_HH
#define DUNE_XT_FUNCTIONS_SPE10_MODEL2_HH

#include <iostream>
#include <memory>

#include <dune/xt/common/color.hh>
#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/common/type_traits.hh>

#include <dune/xt/functions/interfaces/global-function.hh>

namespace Dune {
namespace XT {
namespace Functions {
namespace Spe10 {


/**
 * Grid originally had LL (0,0,0) to UR (365.76, 670.56, 51.816) corners
 *
 */
template <class EntityImp, class DomainFieldImp, size_t dim_domain, class RangeFieldImp, size_t r, size_t rC>
class Model2Function : public GlobalFunctionInterface<EntityImp, DomainFieldImp, dim_domain, RangeFieldImp, r, rC>
{
  static_assert(r == rC, "");
  static_assert(dim_domain == rC, "");
  static_assert(dim_domain == 3, "");
  typedef GlobalFunctionInterface<EntityImp, DomainFieldImp, dim_domain, RangeFieldImp, r, rC> BaseType;

public:
  Model2Function(std::string data_filename = "perm_case2a.dat",
                 Common::FieldVector<double, dim_domain> upper_right = default_upper_right)
    : deltas_{{upper_right[0] / num_elements[0], upper_right[1] / num_elements[1], upper_right[2] / num_elements[2]}}
    , permeability_(nullptr)
    , permMatrix_(0.0)
    , filename_(data_filename)
  {
    readPermeability();
  }

  static const Common::FieldVector<double, dim_domain> default_upper_right;
  // unsigned int mandated by CubeGrid provider
  static const Common::FieldVector<unsigned int, dim_domain> num_elements;

  virtual ~Model2Function()
  {
    delete permeability_;
    permeability_ = nullptr;
  }

  //! currently used in gdt assembler
  virtual void evaluate(const typename BaseType::DomainType& x,
                        typename BaseType::RangeType& diffusion,
                        const Common::Parameter& /*mu*/ = {}) const final override
  {

    if (!permeability_) {
      DXTC_LOG_ERROR_0 << "The SPE10-permeability data file could not be opened. This file does\n"
                       << "not come with the dune-multiscale repository due to file size. To download it\n"
                       << "execute\n"
                       << "wget http://www.spe.org/web/csp/datasets/por_perm_case2a.zip\n"
                       << "unzip the file and move the file 'spe_perm.dat' to\n"
                       << "dune-multiscale/dune/multiscale/problems/spe10_permeability.dat!\n";
      DUNE_THROW(IOError, "Data file for Groundwaterflow permeability could not be opened!");
    }

    // 3 is the maximum space dimension
    for (size_t dim = 0; dim < dim_domain; ++dim)
      permIntervalls_[dim] = std::min(unsigned(std::floor(x[dim] / deltas_[dim])), num_elements[dim] - 1);

    const int offset = permIntervalls_[0] + permIntervalls_[1] * num_elements[0]
                       + permIntervalls_[2] * num_elements[1] * num_elements[0];
    for (size_t dim = 0; dim < dim_domain; ++dim) {
      const auto idx = offset + dim * 1122000;
      diffusion[dim][dim] = permeability_[idx];
    }
  }

  virtual size_t order(const XT::Common::Parameter& /*mu*/ = {}) const override
  {
    return 0u;
  }

private:
  void readPermeability()
  {
    std::ifstream file(filename_);
    double val;
    if (!file) { // file couldn't be opened
      return;
    }
    file >> val;
    int counter = 0;
    permeability_ = new double[3366000];
    while (!file.eof()) {
      // keep reading until end-of-file
      permeability_[counter++] = val;
      file >> val; // sets EOF flag if no value found
    }
    file.close();
  }

  std::array<double, dim_domain> deltas_;
  double* permeability_; //! TODO automatic memory
  mutable typename BaseType::DomainType permIntervalls_;
  mutable Dune::FieldMatrix<double, BaseType::DomainType::dimension, BaseType::DomainType::dimension> permMatrix_;
  const std::string filename_;
};

template <class EntityImp, class DomainFieldImp, size_t dim_domain, class RangeFieldImp, size_t r, size_t rC>
const Common::FieldVector<unsigned int, dim_domain>
    Model2Function<EntityImp, DomainFieldImp, dim_domain, RangeFieldImp, r, rC>::num_elements{{60, 220, 85}};
template <class EntityImp, class DomainFieldImp, size_t dim_domain, class RangeFieldImp, size_t r, size_t rC>
const Common::FieldVector<double, dim_domain>
    Model2Function<EntityImp, DomainFieldImp, dim_domain, RangeFieldImp, r, rC>::default_upper_right{{1, 3.667, 1.417}};


} // namespace Spe10
} // namespace Functions
} // namespace XT
} // namespace Dune


#include "model2.lib.hh"


#endif // DUNE_XT_FUNCTIONS_SPE10_MODEL2_HH
