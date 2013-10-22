// This file is part of the dune-stuff project:
//   http://users.dune-project.org/projects/dune-stuff/
// Copyright Holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_FUNCTIONS_SPE10_HH
#define DUNE_STUFF_FUNCTIONS_SPE10_HH

#include <iostream>
#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/color.hh>
#include <dune/stuff/common/string.hh>

#include "checkerboard.hh"


namespace Dune {
namespace Stuff {
namespace Function {


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

  static std::string static_id();

private:
  static const size_t numXelements;
  static const size_t numYelements;
  static const size_t numZelements;
  static const RangeFieldType minValue;
  static const RangeFieldType maxValue;

  static std::vector<RangeType> read_values_from_file(const std::string& filename, const RangeFieldType& min,
                                                      const RangeFieldType& max);

public:
  Spe10Model1(const std::string& filename, std::vector<DomainFieldType>&& lowerLeft,
              std::vector<DomainFieldType>&& upperRight, const RangeFieldType min = minValue,
              const RangeFieldType max = maxValue, const std::string nm = static_id());

  virtual ThisType* copy() const DS_OVERRIDE;

  static Dune::ParameterTree defaultSettings(const std::string subName = "");

  static ThisType* create(const Dune::Stuff::Common::ExtendedParameterTree settings = defaultSettings());
}; // class Spe10Model1< ..., 2, ..., 1, 1 >


} // namespace Function
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTIONS_SPE10_HH
