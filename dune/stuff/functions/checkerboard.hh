// This file is part of the dune-stuff project:
//   http://users.dune-project.org/projects/dune-stuff/
// Copyright Holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_FUNCTION_CHECKERBOARD_HH
#define DUNE_STUFF_FUNCTION_CHECKERBOARD_HH

#include <vector>
#include <cmath>

#include <dune/common/exceptions.hh>
#include <dune/common/static_assert.hh>

#include <dune/stuff/common/parameter/tree.hh>

#include "interfaces.hh"

namespace Dune {
namespace Stuff {
namespace Function {


template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class Checkerboard
    : public LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
{
  typedef LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
      BaseType;
  typedef Checkerboard<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> ThisType;

  class Localfunction
      : public LocalfunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
  {
    typedef LocalfunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
        BaseType;

  public:
    typedef typename BaseType::EntityType EntityType;

    typedef typename BaseType::DomainFieldType DomainFieldType;
    static const unsigned int dimDomain = BaseType::dimDomain;
    typedef typename BaseType::DomainType DomainType;

    typedef typename BaseType::RangeFieldType RangeFieldType;
    static const unsigned int dimRange     = BaseType::dimRange;
    static const unsigned int dimRangeCols = BaseType::dimRangeCols;
    typedef typename BaseType::RangeType RangeType;

    typedef typename BaseType::JacobianRangeType JacobianRangeType;

    Localfunction(const EntityType& ent, const RangeType value)
      : BaseType(ent)
      , value_(value)
    {
    }

    Localfunction(const Localfunction& /*other*/) = delete;

    Localfunction& operator=(const Localfunction& /*other*/) = delete;

    virtual size_t order() const DS_OVERRIDE
    {
      return 0;
    }

    virtual void evaluate(const DomainType& xx, RangeType& ret) const DS_OVERRIDE
    {
      assert(this->is_a_valid_point(xx));
      ret = value_;
    }

    virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const DS_OVERRIDE
    {
      assert(this->is_a_valid_point(xx));
      ret *= RangeFieldType(0);
    }

  private:
    const RangeType value_;
  }; // class Localfunction

public:
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::LocalfunctionType LocalfunctionType;

  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const int dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType DomainType;

  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const int dimRange     = BaseType::dimRange;
  static const int dimRangeCols = BaseType::dimRangeCols;
  typedef typename BaseType::RangeType RangeType;

  typedef typename BaseType::JacobianRangeType JacobianRangeType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".checkerboard";
  }

  static Dune::ParameterTree defaultSettings(const std::string subName = "")
  {
    Dune::ParameterTree description;
    description["lowerLeft"]   = "[0.0; 0.0; 0.0]";
    description["upperRight"]  = "[1.0; 1.0; 1.0]";
    description["numElements"] = "[2; 2; 2]";
    description["values"]      = "[1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0]";
    description["name"] = static_id();
    if (subName.empty())
      return description;
    else {
      Dune::Stuff::Common::ExtendedParameterTree extendedDescription;
      extendedDescription.add(description, subName);
      return extendedDescription;
    }
  } // ... defaultSettings(...)

  static ThisType* create(const DSC::ExtendedParameterTree settings = defaultSettings())
  {
    // get data
    std::vector<DomainFieldType> lowerLeft  = settings.getVector("lowerLeft", DomainFieldType(0), dimDomain);
    std::vector<DomainFieldType> upperRight = settings.getVector("upperRight", DomainFieldType(1), dimDomain);
    std::vector<size_t> numElements         = settings.getVector("numElements", size_t(1), dimDomain);
    std::vector<RangeFieldType> values_rf = settings.getVector("values", RangeFieldType(1), 1);
    std::vector<RangeType> values(values_rf.size());
    for (size_t ii = 0; ii < values_rf.size(); ++ii)
      values[ii] = RangeType(values_rf[ii]);
    // create and return, leave the checks to the base constructor
    return new ThisType(std::move(lowerLeft), std::move(upperRight), std::move(numElements), std::move(values));
  } // ... create(...)

  Checkerboard(std::vector<DomainFieldType>&& lowerLeft, std::vector<DomainFieldType>&& upperRight,
               std::vector<size_t>&& numElements, std::vector<RangeType>&& values, std::string nm = static_id())
    : lowerLeft_(new std::vector<DomainFieldType>(std::move(lowerLeft)))
    , upperRight_(new std::vector<DomainFieldType>(std::move(upperRight)))
    , numElements_(new std::vector<size_t>(std::move(numElements)))
    , values_(new std::vector<RangeType>(std::move(values)))
    , name_(nm)
  {
    // checks
    if (lowerLeft_->size() < dimDomain)
      DUNE_THROW(Dune::RangeError,
                 "lowerLeft too small (is " << lowerLeft_->size() << ", should be " << dimDomain << ")");
    if (upperRight_->size() < dimDomain)
      DUNE_THROW(Dune::RangeError,
                 "upperRight too small (is " << upperRight_->size() << ", should be " << dimDomain << ")");
    if (numElements_->size() < dimDomain)
      DUNE_THROW(Dune::RangeError,
                 "numElements too small (is " << numElements_->size() << ", should be " << dimDomain << ")");
    size_t totalSubdomains = 1;
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      const auto& ll = (*lowerLeft_)[dd];
      const auto& ur = (*upperRight_)[dd];
      const auto& ne = (*numElements_)[dd];
      if (!(ll < ur))
        DUNE_THROW(Dune::RangeError, "lowerLeft has to be elementwise smaller than upperRight!");
      totalSubdomains *= ne;
    }
    if (values_->size() < totalSubdomains)
      DUNE_THROW(Dune::RangeError,
                 "values too small (is " << values_->size() << ", should be " << totalSubdomains << ")");
  } // Checkerboard(...)

  Checkerboard(const ThisType& other)
    : lowerLeft_(other.lowerLeft_)
    , upperRight_(other.upperRight_)
    , numElements_(other.numElements_)
    , values_(other.values_)
    , name_(other.name_)
  {
  }

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      lowerLeft_   = other.lowerLeft_;
      upperRight_  = other.upperRight_;
      numElements_ = other.numElements_;
      values_      = other.values_;
      name_        = other.name_;
    }
    return *this;
  }

  virtual ThisType* copy() const DS_OVERRIDE
  {
    return new ThisType(*this);
  }

  virtual std::string name() const DS_OVERRIDE
  {
    return name_;
  }

  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const DS_OVERRIDE
  {
    // decide on the subdomain the center of the entity belongs to
    const auto center = entity.geometry().center();
    std::vector<size_t> whichPartition(dimDomain, 0);
    const auto& ll = *lowerLeft_;
    const auto& ur = *upperRight_;
    const auto& ne = *numElements_;
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      // for points that are on upperRight_[d], this selects one partition too much
      // so we need to cap this
      whichPartition[dd] =
          std::min(size_t(std::floor(ne[dd] * ((center[dd] - ll[dd]) / (ur[dd] - ll[dd])))), ne[dd] - 1);
    }
    size_t subdomain = 0;
    if (dimDomain == 1)
      subdomain = whichPartition[0];
    else if (dimDomain == 2)
      subdomain = whichPartition[0] + whichPartition[1] * ne[0];
    else
      subdomain = whichPartition[0] + whichPartition[1] * ne[0] + whichPartition[2] * ne[1] * ne[0];
    // return the component that belongs to the subdomain
    return std::unique_ptr<Localfunction>(new Localfunction(entity, (*values_)[subdomain]));
  } // ... local_function(...)

private:
  std::shared_ptr<const std::vector<DomainFieldType>> lowerLeft_;
  std::shared_ptr<const std::vector<DomainFieldType>> upperRight_;
  std::shared_ptr<const std::vector<size_t>> numElements_;
  std::shared_ptr<const std::vector<RangeType>> values_;
  std::string name_;
}; // class Checkerboard


} // namespace Function
} // namespace Stuff
} // namespace Dune

#define DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGE(etype, ddim)                                                   \
  DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGECOLS(etype, ddim, 1)                                                  \
  DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGECOLS(etype, ddim, 2)                                                  \
  DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGECOLS(etype, ddim, 3)

#define DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGECOLS(etype, ddim, rdim)                                         \
  DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DOMAINFIELDTYPES(etype, ddim, rdim, 1)                                        \
  DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DOMAINFIELDTYPES(etype, ddim, rdim, 2)                                        \
  DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DOMAINFIELDTYPES(etype, ddim, rdim, 3)

#define DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DOMAINFIELDTYPES(etype, ddim, rdim, rcdim)                              \
  DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_RANGEFIELDTYPES(etype, double, ddim, rdim, rcdim)

#define DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_RANGEFIELDTYPES(etype, dftype, ddim, rdim, rcdim)                       \
  DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LAST_EXPANSION(etype, dftype, ddim, double, rdim, rcdim)                           \
  DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LAST_EXPANSION(etype, dftype, ddim, long double, rdim, rcdim)

#define DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LAST_EXPANSION(etype, dftype, ddim, rftype, rdim, rcdim)                     \
  extern template class Dune::Stuff::Function::Checkerboard<etype, dftype, ddim, rftype, rdim, rcdim>;

#ifdef HAVE_DUNE_GRID

#include <dune/grid/sgrid.hh>

typedef typename Dune::SGrid<1, 1>::template Codim<0>::Entity DuneFunctionsCheckerboardSGrid1dEntityType;
typedef typename Dune::SGrid<2, 2>::template Codim<0>::Entity DuneFunctionsCheckerboardSGrid2dEntityType;
typedef typename Dune::SGrid<3, 3>::template Codim<0>::Entity DuneFunctionsCheckerboardSGrid3dEntityType;

DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGE(DuneFunctionsCheckerboardSGrid1dEntityType, 1)
DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGE(DuneFunctionsCheckerboardSGrid2dEntityType, 2)
DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGE(DuneFunctionsCheckerboardSGrid3dEntityType, 3)

#include <dune/grid/yaspgrid.hh>

typedef typename Dune::YaspGrid<1>::template Codim<0>::Entity DuneFunctionsCheckerboardYaspGrid1dEntityType;
typedef typename Dune::YaspGrid<2>::template Codim<0>::Entity DuneFunctionsCheckerboardYaspGrid2dEntityType;
typedef typename Dune::YaspGrid<3>::template Codim<0>::Entity DuneFunctionsCheckerboardYaspGrid3dEntityType;

DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGE(DuneFunctionsCheckerboardYaspGrid1dEntityType, 1)
DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGE(DuneFunctionsCheckerboardYaspGrid2dEntityType, 2)
DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGE(DuneFunctionsCheckerboardYaspGrid3dEntityType, 3)

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
#ifdef ALUGRID_CONFORM
#define DUNE_STUFF_FUNCTION_INTERFACE_ALUGRID_CONFORM_WAS_DEFINED_BEFORE
#else
#define ALUGRID_CONFORM 1
#endif
#ifdef ENABLE_ALUGRID
#define DUNE_STUFF_FUNCTION_INTERFACE_ENABLE_ALUGRID_WAS_DEFINED_BEFORE
#else
#define ENABLE_ALUGRID 1
#endif

#include <dune/grid/alugrid.hh>

typedef
    typename Dune::ALUSimplexGrid<2, 2>::template Codim<0>::Entity DuneFunctionsCheckerboardAluSimplexGrid2dEntityType;
typedef
    typename Dune::ALUSimplexGrid<3, 3>::template Codim<0>::Entity DuneFunctionsCheckerboardAluSimplexGrid3dEntityType;
typedef typename Dune::ALUCubeGrid<3, 3>::template Codim<0>::Entity DuneFunctionsCheckerboardAluCubeGrid3dEntityType;

DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGE(DuneFunctionsCheckerboardAluSimplexGrid2dEntityType, 2)
DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGE(DuneFunctionsCheckerboardAluSimplexGrid3dEntityType, 3)
DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGE(DuneFunctionsCheckerboardAluCubeGrid3dEntityType, 3)

#ifdef DUNE_STUFF_FUNCTION_INTERFACE_ALUGRID_CONFORM_WAS_DEFINED_BEFORE
#undef DUNE_STUFF_FUNCTION_INTERFACE_ALUGRID_CONFORM_WAS_DEFINED_BEFORE
#else
#undef ALUGRID_CONFORM
#endif
#ifdef DUNE_STUFF_FUNCTION_INTERFACE_ENABLE_ALUGRID_WAS_DEFINED_BEFORE
#undef DUNE_STUFF_FUNCTION_INTERFACE_ENABLE_ALUGRID_WAS_DEFINED_BEFORE
#else
#undef ENABLE_ALUGRID
#endif
#endif // HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
#endif // HAVE_DUNE_GRID

#undef DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LAST_EXPANSION
#undef DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_RANGEFIELDTYPES
#undef DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DOMAINFIELDTYPES
#undef DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGECOLS
#undef DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGE

#endif // DUNE_STUFF_FUNCTION_CHECKERBOARD_HH
