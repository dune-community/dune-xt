// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_FUNCTIONS_SPE10MODEL2_HH
#define DUNE_STUFF_FUNCTIONS_SPE10MODEL2_HH

#include <iostream>
#include <memory>

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/color.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/common/fvector.hh>
#include <dune/stuff/common/type_utils.hh>

namespace Dune {
namespace Stuff {
namespace Functions {
namespace Spe10 {

/**
 * Grid is currently mandated to have LL (0,0,0) to UR (365.76, 670.56, 51.816) corners
 */
template <class EntityImp, class DomainFieldImp, size_t dim_domain, class RangeFieldImp, size_t r, size_t rC>
class Model2 : public Stuff::GlobalFunctionInterface<EntityImp, DomainFieldImp, dim_domain, RangeFieldImp, r, rC>
{
  static_assert(r == rC, "");
  static_assert(dim_domain == rC, "");
  static_assert(dim_domain == 3, "");
  typedef Stuff::GlobalFunctionInterface<EntityImp, DomainFieldImp, dim_domain, RangeFieldImp, r, rC> BaseType;

public:
  Model2(std::string data_filename = "perm_case2a.dat")
    : deltas_{{6.096, 3.048, 0.6096}}
    , permeability_(nullptr)
    , permMatrix_(0.0)
    , filename_(data_filename)
  {
    readPermeability();
  }

  virtual ~Model2()
  {
    delete permeability_;
    permeability_ = nullptr;
  }

  //! currently used in gdt assembler
  virtual void evaluate(const typename BaseType::DomainType& x,
                        typename BaseType::RangeType& diffusion) const final override
  {

    if (!permeability_) {
      DSC_LOG_ERROR_0 << "The SPE10-permeability data file could not be opened. This file does\n"
                      << "not come with the dune-multiscale repository due to file size. To download it\n"
                      << "execute\n"
                      << "wget http://www.spe.org/web/csp/datasets/por_perm_case2a.zip\n"
                      << "unzip the file and move the file 'spe_perm.dat' to\n"
                      << "dune-multiscale/dune/multiscale/problems/spe10_permeability.dat!\n";
      DUNE_THROW(IOError, "Data file for Groundwaterflow permeability could not be opened!");
    }
    // 3 is the maximum space dimension
    for (size_t dim = 0; dim < dim_domain; ++dim)
      permIntervalls_[dim] = std::floor(x[dim] / deltas_[dim]);

    const int offset = permIntervalls_[0] + permIntervalls_[1] * 60 + permIntervalls_[2] * 220 * 60;
    for (size_t dim = 0; dim < dim_domain; ++dim)
      diffusion[dim][dim] = permeability_[offset + dim * 1122000];
  }

  virtual size_t order() const
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
    int counter   = 0;
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

} // namespace Spe10
} // namespace Functions
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTIONS_SPE10MODEL2_HH
