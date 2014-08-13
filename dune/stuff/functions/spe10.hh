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

#include "checkerboard.hh"


namespace Dune {
namespace Stuff {
namespace Functions {


// default, to allow for specialization
template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class Spe10Model1
{
public:
  Spe10Model1() = delete;
}; // class FunctionSpe10Model1


/**
 *  \note For dimRange 1 we only read the Kx values from the file.
 */
template <class EntityImp, class DomainFieldImp, class RangeFieldImp>
class Spe10Model1<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1, 1>
    : public Checkerboard<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1, 1>
{
  typedef Checkerboard<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1, 1> BaseType;
  typedef Spe10Model1<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1, 1> ThisType;

public:
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::LocalfunctionType LocalfunctionType;

  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const int dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType DomainType;

  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const int dimRange = BaseType::dimRange;
  typedef typename BaseType::RangeType RangeType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".spe10.model1";
  }

private:
  static const size_t numXelements;
  static const size_t numYelements;
  static const size_t numZelements;
  static const RangeFieldType minValue;
  static const RangeFieldType maxValue;

  static std::vector<RangeType> read_values_from_file(const std::string& filename, const RangeFieldType& min,
                                                      const RangeFieldType& max)

  {
    if (!(max > min))
      DUNE_THROW(Dune::RangeError, "max (is " << max << ") has to be larger than min (is " << min << ")!");
    const RangeFieldType scale = (max - min) / (maxValue - minValue);
    const RangeType shift      = min - scale * minValue;
    // read all the data from the file
    std::ifstream datafile(filename);
    if (datafile.is_open()) {
      static const size_t entriesPerDim = numXelements * numYelements * numZelements;
      // create storage (there should be exactly 6000 values in the file, but we onyl read the first 2000)
      std::vector<RangeType> data(entriesPerDim, RangeFieldType(0));
      double tmp     = 0;
      size_t counter = 0;
      while (datafile >> tmp && counter < entriesPerDim) {
        data[counter] = (tmp * scale) + shift;
        ++counter;
      }
      datafile.close();
      if (counter != entriesPerDim)
        DUNE_THROW(Dune::IOError,
                   "wrong number of entries in '" << filename << "' (are " << counter << ", should be " << entriesPerDim
                                                  << ")!");
      return data;
    } else
      DUNE_THROW(Dune::IOError, "could not open '" << filename << "'!");
  } // ... read_values_from_file(...)

public:
  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["filename"]    = "perm_case1.dat";
    config["lower_left"]  = "[0.0 0.0]";
    config["upper_right"] = "[762.0 15.24]";
    config["min_value"]   = "0.001";
    config["max_value"]   = "998.915";
    config["name"] = static_id();
    if (sub_name.empty())
      return config;
    else {
      Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr<ThisType> create(const Common::Configuration config = default_config(),
                                          const std::string sub_name = static_id())
  {
    // get correct config
    const Common::Configuration cfg         = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration default_cfg = default_config();
    // create
    return Common::make_unique<ThisType>(
        cfg.get("filename", default_cfg.get<std::string>("filename")),
        cfg.get("lower_left", default_cfg.get<std::vector<DomainFieldType>>("lower_left"), dimDomain),
        cfg.get("upper_right", default_cfg.get<std::vector<DomainFieldType>>("upper_right"), dimDomain),
        cfg.get("min_val", minValue),
        cfg.get("max_val", maxValue),
        cfg.get("name", default_cfg.get<std::string>("name")));
  } // ... create(...)

  Spe10Model1(const std::string& filename, std::vector<DomainFieldType>&& lowerLeft,
              std::vector<DomainFieldType>&& upperRight, const RangeFieldType min = minValue,
              const RangeFieldType max = maxValue, const std::string nm = static_id())
    : BaseType(std::move(lowerLeft), std::move(upperRight), {numXelements, numZelements},
               read_values_from_file(filename, min, max), nm)
  {
  }

  virtual std::string type() const DS_OVERRIDE
  {
    return BaseType::static_id() + ".spe10.model1";
  }

  virtual ThisType* copy() const DS_OVERRIDE
  {
    return new ThisType(*this);
  }
}; // class Spe10Model1< ..., 2, ..., 1, 1 >


/**
 * We read only the Kx values from file and scale the unit matrix. Thus we should be the same as the scalar variant.
 */
template <class EntityImp, class DomainFieldImp, class RangeFieldImp>
class Spe10Model1<EntityImp, DomainFieldImp, 2, RangeFieldImp, 2, 2>
    : public Checkerboard<EntityImp, DomainFieldImp, 2, RangeFieldImp, 2, 2>
{
  typedef Checkerboard<EntityImp, DomainFieldImp, 2, RangeFieldImp, 2, 2> BaseType;
  typedef Spe10Model1<EntityImp, DomainFieldImp, 2, RangeFieldImp, 2, 2> ThisType;

public:
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::LocalfunctionType LocalfunctionType;

  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const int dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType DomainType;

  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const int dimRange = BaseType::dimRange;
  typedef typename BaseType::RangeType RangeType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".spe10.model1";
  }

private:
  static const size_t numXelements;
  static const size_t numYelements;
  static const size_t numZelements;
  static const RangeFieldType minValue;
  static const RangeFieldType maxValue;

  static std::vector<RangeType> read_values_from_file(const std::string& filename, const RangeFieldType& min,
                                                      const RangeFieldType& max)

  {
    if (!(max > min))
      DUNE_THROW(Dune::RangeError, "max (is " << max << ") has to be larger than min (is " << min << ")!");
    const RangeFieldType scale = (max - min) / (maxValue - minValue);
    const RangeFieldType shift = min - scale * minValue;
    // read all the data from the file
    std::ifstream datafile(filename);
    if (datafile.is_open()) {
      static const size_t entriesPerDim = numXelements * numYelements * numZelements;
      // create storage (there should be exactly 6000 values in the file, but we onyl read the first 2000)
      RangeType unit_matrix(0.0);
      for (size_t dd = 0; dd < dimDomain; ++dd)
        unit_matrix[dd][dd] = 1.0;
      std::vector<RangeType> data(entriesPerDim, unit_matrix);
      double tmp     = 0;
      size_t counter = 0;
      while (datafile >> tmp && counter < entriesPerDim) {
        data[counter] *= (tmp * scale) + shift;
        ++counter;
      }
      datafile.close();
      if (counter != entriesPerDim)
        DUNE_THROW(Dune::IOError,
                   "wrong number of entries in '" << filename << "' (are " << counter << ", should be " << entriesPerDim
                                                  << ")!");
      return data;
    } else
      DUNE_THROW(Dune::IOError, "could not open '" << filename << "'!");
  } // ... read_values_from_file(...)

public:
  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["filename"]    = "perm_case1.dat";
    config["lower_left"]  = "[0.0 0.0]";
    config["upper_right"] = "[762.0 15.24]";
    config["min_value"]   = "0.001";
    config["max_value"]   = "998.915";
    config["name"] = static_id();
    if (sub_name.empty())
      return config;
    else {
      Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr<ThisType> create(const Common::Configuration config = default_config(),
                                          const std::string sub_name = static_id())
  {
    // get correct config
    const Common::Configuration cfg         = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration default_cfg = default_config();
    // create
    return Common::make_unique<ThisType>(
        cfg.get("filename", default_cfg.get<std::string>("filename")),
        cfg.get("lower_left", default_cfg.get<std::vector<DomainFieldType>>("lower_left"), dimDomain),
        cfg.get("upper_right", default_cfg.get<std::vector<DomainFieldType>>("upper_right"), dimDomain),
        cfg.get("min_val", minValue),
        cfg.get("max_val", maxValue),
        cfg.get("name", default_cfg.get<std::string>("name")));
  } // ... create(...)

  Spe10Model1(const std::string& filename, std::vector<DomainFieldType>&& lowerLeft,
              std::vector<DomainFieldType>&& upperRight, const RangeFieldType min = minValue,
              const RangeFieldType max = maxValue, const std::string nm = static_id())
    : BaseType(std::move(lowerLeft), std::move(upperRight), {numXelements, numZelements},
               read_values_from_file(filename, min, max), nm)
  {
  }

  virtual std::string type() const DS_OVERRIDE
  {
    return BaseType::static_id() + ".spe10.model1";
  }

  virtual ThisType* copy() const DS_OVERRIDE
  {
    return new ThisType(*this);
  }
}; // class Spe10Model1< ..., 2, ..., 1, 1 >


template <class E, class D, class R>
const size_t Spe10Model1<E, D, 2, R, 1, 1>::numXelements = 100;

template <class E, class D, class R>
const size_t Spe10Model1<E, D, 2, R, 1, 1>::numYelements = 1;

template <class E, class D, class R>
const size_t Spe10Model1<E, D, 2, R, 1, 1>::numZelements = 20;

template <class E, class D, class R>
const typename Spe10Model1<E, D, 2, R, 1, 1>::RangeFieldType Spe10Model1<E, D, 2, R, 1, 1>::minValue = 0.001;

template <class E, class D, class R>
const typename Spe10Model1<E, D, 2, R, 1, 1>::RangeFieldType Spe10Model1<E, D, 2, R, 1, 1>::maxValue = 998.915;


template <class E, class D, class R>
const size_t Spe10Model1<E, D, 2, R, 2, 2>::numXelements = 100;

template <class E, class D, class R>
const size_t Spe10Model1<E, D, 2, R, 2, 2>::numYelements = 1;

template <class E, class D, class R>
const size_t Spe10Model1<E, D, 2, R, 2, 2>::numZelements = 20;

template <class E, class D, class R>
const typename Spe10Model1<E, D, 2, R, 2, 2>::RangeFieldType Spe10Model1<E, D, 2, R, 2, 2>::minValue = 0.001;

template <class E, class D, class R>
const typename Spe10Model1<E, D, 2, R, 2, 2>::RangeFieldType Spe10Model1<E, D, 2, R, 2, 2>::maxValue = 998.915;


} // namespace Functions
} // namespace Stuff
} // namespace Dune

#if 0 // def DUNE_STUFF_FUNCTIONS_TO_LIB
the whole  DUNE_STUFF_FUNCTIONS_TO_LIB stuff threw loads of errors, so I disabled it everywhere (Rene)

#define DUNE_STUFF_FUNCTIONS_SPE10_LIST_DOMAINFIELDTYPES(etype, ddim, rdim, rcdim)                                     \
  DUNE_STUFF_FUNCTIONS_SPE10_LIST_RANGEFIELDTYPES(etype, double, ddim, rdim, rcdim)

#define DUNE_STUFF_FUNCTIONS_SPE10_LIST_RANGEFIELDTYPES(etype, dftype, ddim, rdim, rcdim)                              \
  DUNE_STUFF_FUNCTIONS_SPE10_LAST_EXPANSION(etype, dftype, ddim, double, rdim, rcdim)                                  \
  DUNE_STUFF_FUNCTIONS_SPE10_LAST_EXPANSION(etype, dftype, ddim, long double, rdim, rcdim)

#define DUNE_STUFF_FUNCTIONS_SPE10_LAST_EXPANSION(etype, dftype, ddim, rftype, rdim, rcdim)                            \
  extern template class Dune::Stuff::Functions::Spe10Model1<etype, dftype, ddim, rftype, rdim, rcdim>;

#if HAVE_DUNE_GRID

DUNE_STUFF_FUNCTIONS_SPE10_LIST_DOMAINFIELDTYPES(DuneStuffFunctionsInterfacesSGrid2dEntityType, 2, 1, 1)

DUNE_STUFF_FUNCTIONS_SPE10_LIST_DOMAINFIELDTYPES(DuneStuffFunctionsInterfacesYaspGrid2dEntityType, 2, 1, 1)

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H

DUNE_STUFF_FUNCTIONS_SPE10_LIST_DOMAINFIELDTYPES(DuneStuffFunctionsInterfacesAluSimplexGrid2dEntityType, 2, 1, 1)

#endif // HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
#endif // HAVE_DUNE_GRID

#undef DUNE_STUFF_FUNCTIONS_SPE10_LAST_EXPANSION
#undef DUNE_STUFF_FUNCTIONS_SPE10_LIST_RANGEFIELDTYPES
#undef DUNE_STUFF_FUNCTIONS_SPE10_LIST_DOMAINFIELDTYPES
#endif // DUNE_STUFF_FUNCTIONS_TO_LIB

#endif // DUNE_STUFF_FUNCTIONS_SPE10_HH
