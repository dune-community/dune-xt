// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_FUNCTIONS_SPE10_HH
#define DUNE_STUFF_FUNCTIONS_SPE10_HH

#include <iostream>
#include <memory>

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/color.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/common/fvector.hh>
#include <dune/stuff/common/type_utils.hh>

#include "checkerboard.hh"


namespace Dune {
namespace Stuff {
namespace Exceptions {


class spe10_data_file_missing : public Dune::IOError
{
};


} // namespace Exceptions
namespace Functions {
namespace Spe10 {
namespace internal {


static const std::string model1_filename = "perm_case1.dat";
static const size_t model1_x_elements    = 100;
static const size_t model1_y_elements    = 1;
static const size_t model1_z_elements    = 20;
static const double model_1_length_x     = 762.0;
static const double model_1_length_y     = 7.62;
static const double model_1_length_z     = 15.24;
static const double model1_min_value     = 0.001;
static const double model1_max_value     = 998.915;


template <class EntityImp, class DomainFieldImp, class RangeFieldImp, size_t r, size_t rC>
class Model1Base : public Checkerboard<EntityImp, DomainFieldImp, 2, RangeFieldImp, r, rC>
{
  typedef Checkerboard<EntityImp, DomainFieldImp, 2, RangeFieldImp, r, rC> BaseType;

public:
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::LocalfunctionType LocalfunctionType;

  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeFieldType RangeFieldType;
  typedef typename BaseType::RangeType RangeType;

  static const bool available = true;

  static std::string static_id()
  {
    return LocalizableFunctionInterface<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1, 1>::static_id()
           + ".spe10.model1";
  } // ... static_id(...)

private:
  static std::vector<RangeType> read_values_from_file(const std::string& filename, const RangeFieldType& min,
                                                      const RangeFieldType& max, const RangeType& unit_range)

  {
    if (!(max > min))
      DUNE_THROW(Dune::RangeError, "max (is " << max << ") has to be larger than min (is " << min << ")!");
    const RangeFieldType scale = (max - min) / (internal::model1_max_value - internal::model1_min_value);
    const RangeFieldType shift = min - scale * internal::model1_min_value;
    // read all the data from the file
    std::ifstream datafile(filename);
    if (datafile.is_open()) {
      static const size_t entriesPerDim = model1_x_elements * model1_y_elements * model1_z_elements;
      // create storage (there should be exactly 6000 values in the file, but we onyl read the first 2000)
      std::vector<RangeType> data(entriesPerDim, unit_range);
      double tmp     = 0;
      size_t counter = 0;
      while (datafile >> tmp && counter < entriesPerDim)
        data[counter++] *= (tmp * scale) + shift;
      datafile.close();
      if (counter != entriesPerDim)
        DUNE_THROW(Dune::IOError,
                   "wrong number of entries in '" << filename << "' (are " << counter << ", should be " << entriesPerDim
                                                  << ")!");
      return data;
    } else
      DUNE_THROW(Exceptions::spe10_data_file_missing, "could not open '" << filename << "'!");
  } // ... read_values_from_file(...)

public:
  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["filename"]   = internal::model1_filename;
    config["lower_left"] = "[0.0 0.0]";
    config["upper_right"] =
        "[" + Common::toString(internal::model_1_length_x) + " " + Common::toString(internal::model_1_length_z) + "]";
    config.set("min_value", internal::model1_min_value);
    config.set("max_value", internal::model1_max_value);
    config["name"] = static_id();
    if (sub_name.empty())
      return config;
    else {
      Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  template <class DerivedType>
  static std::unique_ptr<DerivedType> create_derived(const Common::Configuration config, const std::string sub_name)
  {
    // get correct config
    const Common::Configuration cfg         = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration default_cfg = default_config();
    // create
    return Common::make_unique<DerivedType>(cfg.get("filename", default_cfg.get<std::string>("filename")),
                                            cfg.get("lower_left", default_cfg.get<DomainType>("lower_left")),
                                            cfg.get("upper_right", default_cfg.get<DomainType>("upper_right")),
                                            cfg.get("min_value", default_cfg.get<RangeFieldType>("min_value")),
                                            cfg.get("max_value", default_cfg.get<RangeFieldType>("max_value")),
                                            cfg.get("name", default_cfg.get<std::string>("name")));
  } // ... create(...)

  Model1Base(const std::string& filename, const DomainType& lowerLeft, const DomainType& upperRight,
             const RangeFieldType min, const RangeFieldType max, const std::string nm, const RangeType& unit_range)
    : BaseType(lowerLeft, upperRight, {model1_x_elements, model1_z_elements},
               read_values_from_file(filename, min, max, unit_range), nm)
  {
  }

  virtual std::string type() const override
  {
    return BaseType::static_id() + ".spe10.model1";
  }
}; // class Model1Base


} // namespace internal


// default, to allow for specialization
template <class E, class D, size_t d, class R, size_t r, size_t rC = 1>
class Model1 : public LocalizableFunctionInterface<E, D, d, R, r, rC>
{
  Model1()
  {
    static_assert(AlwaysFalse<E>::value, "Not available for these dimensions!");
  }
};


/**
 * We read only the Kx values from file and scale the unit matrix atm.
 */
template <class EntityImp, class DomainFieldImp, class RangeFieldImp, size_t r>
class Model1<EntityImp, DomainFieldImp, 2, RangeFieldImp, r, r>
    : public internal::Model1Base<EntityImp, DomainFieldImp, RangeFieldImp, r, r>
{
  typedef internal::Model1Base<EntityImp, DomainFieldImp, RangeFieldImp, r, r> BaseType;
  typedef Model1<EntityImp, DomainFieldImp, 2, RangeFieldImp, r, r> ThisType;

public:
  using typename BaseType::DomainFieldType;
  using BaseType::dimDomain;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;

  static std::unique_ptr<ThisType> create(const Common::Configuration config = BaseType::default_config(),
                                          const std::string sub_name = BaseType::static_id())
  {
    return BaseType::template create_derived<ThisType>(config, sub_name);
  } // ... create(...)

  Model1(const std::string& filename, const Common::FieldVector<DomainFieldType, dimDomain>& lower_left,
         const Common::FieldVector<DomainFieldType, dimDomain>& upper_right,
         const RangeFieldType min = internal::model1_min_value, const RangeFieldType max = internal::model1_max_value,
         const std::string nm = BaseType::static_id())
    : BaseType(filename, lower_left, upper_right, min, max, nm, unit_matrix())
  {
  }

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

#endif // DUNE_STUFF_FUNCTIONS_SPE10_HH
