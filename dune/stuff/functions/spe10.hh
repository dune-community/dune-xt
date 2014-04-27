// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_FUNCTIONS_SPE10_HH
#define DUNE_STUFF_FUNCTIONS_SPE10_HH

#include <iostream>
#include <memory>

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/parameter/tree.hh>
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
      DUNE_THROW_COLORFULLY(Dune::RangeError, "max (is " << max << ") has to be larger than min (is " << min << ")!");
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
        DUNE_THROW_COLORFULLY(Dune::IOError,
                              "wrong number of entries in '" << filename << "' (are " << counter << ", should be "
                                                             << entriesPerDim
                                                             << ")!");
      return data;
    } else
      DUNE_THROW_COLORFULLY(Dune::IOError, "could not open '" << filename << "'!");
  } // Spe10Model1()

public:
  Spe10Model1(const std::string& filename, std::vector<DomainFieldType>&& lowerLeft,
              std::vector<DomainFieldType>&& upperRight, const RangeFieldType min = minValue,
              const RangeFieldType max = maxValue, const std::string nm = static_id())
    : BaseType(std::move(lowerLeft), std::move(upperRight), {numXelements, numZelements},
               read_values_from_file(filename, min, max), nm)
  {
  }

  virtual ThisType* copy() const DS_OVERRIDE
  {
    return new ThisType(*this);
  }

  static Dune::ParameterTree defaultSettings(const std::string subName = "")
  {
    Dune::ParameterTree description;
    description["filename"]   = "perm_case1.dat";
    description["lowerLeft"]  = "[0.0; 0.0]";
    description["upperRight"] = "[762.0; 15.24]";
    description["minValue"]   = "0.001";
    description["maxValue"]   = "998.915";
    description["name"] = static_id();
    if (subName.empty())
      return description;
    else {
      Dune::Stuff::Common::ExtendedParameterTree extendedDescription;
      extendedDescription.add(description, subName);
      return extendedDescription;
    }
  } // ... defaultSettings(...)

  static ThisType* create(const Dune::Stuff::Common::ExtendedParameterTree settings = defaultSettings())
  {
    // get data
    const std::string filename              = settings.get<std::string>("filename");
    std::vector<DomainFieldType> lowerLeft  = settings.getVector<DomainFieldType>("lowerLeft", dimDomain);
    std::vector<DomainFieldType> upperRight = settings.getVector<DomainFieldType>("upperRight", dimDomain);
    const RangeFieldType minVal             = settings.get<RangeFieldType>("minValue", minValue);
    const RangeFieldType maxVal             = settings.get<RangeFieldType>("maxValue", maxValue);
    const std::string nm                    = settings.get<std::string>("name", static_id());
    // create and return, leave the checks to the constructor
    return new ThisType(filename, std::move(lowerLeft), std::move(upperRight), minVal, maxVal, nm);
  } // ... create(...)
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


} // namespace Functions
} // namespace Stuff
} // namespace Dune

#ifdef DUNE_STUFF_FUNCTIONS_TO_LIB
#define DUNE_STUFF_FUNCTIONS_SPE10_LIST_DOMAINFIELDTYPES(etype, ddim, rdim, rcdim)                                     \
  DUNE_STUFF_FUNCTIONS_SPE10_LIST_RANGEFIELDTYPES(etype, double, ddim, rdim, rcdim)

#define DUNE_STUFF_FUNCTIONS_SPE10_LIST_RANGEFIELDTYPES(etype, dftype, ddim, rdim, rcdim)                              \
  DUNE_STUFF_FUNCTIONS_SPE10_LAST_EXPANSION(etype, dftype, ddim, double, rdim, rcdim)                                  \
  DUNE_STUFF_FUNCTIONS_SPE10_LAST_EXPANSION(etype, dftype, ddim, long double, rdim, rcdim)

#define DUNE_STUFF_FUNCTIONS_SPE10_LAST_EXPANSION(etype, dftype, ddim, rftype, rdim, rcdim)                            \
  extern template class Dune::Stuff::Function::Spe10Model1<etype, dftype, ddim, rftype, rdim, rcdim>;

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
