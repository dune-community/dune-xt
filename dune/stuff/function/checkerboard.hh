#ifndef DUNE_STUFF_FUNCTION_CHECKERBOARD_HH
#define DUNE_STUFF_FUNCTION_CHECKERBOARD_HH

#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/static_assert.hh>
#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/color.hh>

#include "interface.hh"

namespace Dune {
namespace Stuff {


template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols>
class FunctionCheckerboardBase
    : public FunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols>
{
  typedef FunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols> BaseType;

public:
  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const int dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const int dimRangeRows = BaseType::dimRangeRows;
  static const int dimRangeCols = BaseType::dimRangeCols;
  typedef typename BaseType::RangeType RangeType;

  static const std::string id()
  {
    return BaseType::id() + ".checkerboard";
  }

  FunctionCheckerboardBase(const DomainType _lowerLeft, const DomainType _upperRight,
                           const std::vector<size_t> _numElements, const std::vector<RangeType> _values,
                           const std::string _name = id())
    : lowerLeft_(_lowerLeft)
    , upperRight_(_upperRight)
    , numElements_(_numElements)
    , values_(_values)
    , name_(_name)
  {
    // checks
    dune_static_assert((dimDomain > 0), "Really?");
    dune_static_assert((dimDomain <= 3), "Not implemented!");
    dune_static_assert((dimRangeRows > 0), "Really?");
    // get total number of subdomains
    size_t totalSubdomains = 1;
    for (int dd = 0; dd < dimDomain; ++dd) {
      totalSubdomains *= numElements_[dd];
    }
    if (values_.size() < totalSubdomains)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                      << " please provide at least as many '_values' as subdomains given by '_numElements'!");
  } // FunctionCheckerboard(...)

  const DomainType& lowerLeft() const
  {
    return lowerLeft_;
  }

  const DomainType& upperRight() const
  {
    return upperRight_;
  }

  const std::vector<size_t>& numElements() const
  {
    return numElements_;
  }

  const std::vector<RangeType>& values() const
  {
    return values_;
  }

  virtual std::string name() const
  {
    return name_;
  }

  virtual int order() const
  {
    return 0;
  }

  /**
   *  \todo put this stuff in an expression, that is expanded at compile time (dimension dependent)
   */
  virtual void evaluate(const DomainType& x, RangeType& ret) const
  {
    // decide on the subdomain the point x belongs to
    std::vector<size_t> whichPartition;
    for (int d = 0; d < dimDomain; ++d) {
      // for points that are on upperRight_[d], this selects one partition too much
      // so we need to cap this
      whichPartition.push_back(
          std::min(size_t(std::floor(numElements_[d] * ((x[d] - lowerLeft_[d]) / (upperRight_[d] - lowerLeft_[d])))),
                   numElements_[d] - 1));
    }
    size_t subdomain = 0;
    if (dimDomain == 1)
      subdomain = whichPartition[0];
    else if (dimDomain == 2)
      subdomain = whichPartition[0] + whichPartition[1] * numElements_[0];
    else // has to be 3, see checks in constructor
      subdomain = whichPartition[0] + whichPartition[1] * numElements_[0]
                  + whichPartition[2] * numElements_[1] * numElements_[0];
    // return the component that belongs to the subdomain of x
    ret = values_[subdomain];
  } // virtual void evaluate(const DomainType& x, RangeType& ret) const

private:
  DomainType lowerLeft_;
  DomainType upperRight_;
  std::vector<size_t> numElements_;
  std::vector<RangeType> values_;
  std::string name_;
}; // class FunctionCheckerboardBase


// forward, to allow for specialization
template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols = 1>
class FunctionCheckerboard
{
public:
  FunctionCheckerboard() = delete;
};


template <class DomainFieldImp, int domainDim, class RangeFieldImp>
class FunctionCheckerboard<DomainFieldImp, domainDim, RangeFieldImp, 1, 1>
    : public FunctionCheckerboardBase<DomainFieldImp, domainDim, RangeFieldImp, 1, 1>
{
  typedef FunctionCheckerboardBase<DomainFieldImp, domainDim, RangeFieldImp, 1, 1> BaseType;

public:
  typedef FunctionCheckerboard<DomainFieldImp, domainDim, RangeFieldImp, 1, 1> ThisType;

  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeFieldType RangeFieldType;
  typedef typename BaseType::RangeType RangeType;

  using BaseType::id;

  FunctionCheckerboard(const DomainType _lowerLeft, const DomainType _upperRight,
                       const std::vector<size_t> _numElements, const std::vector<RangeType> _values,
                       const std::string _name = id())
    : BaseType(_lowerLeft, _upperRight, _numElements, _values, _name)
  {
  }

  static Dune::ParameterTree defaultSettings(const std::string subName = "")
  {
    Dune::ParameterTree description;
    description["lowerLeft"]   = "[0.0; 0.0; 0.0]";
    description["upperRight"]  = "[1.0; 1.0; 1.0]";
    description["numElements"] = "[2; 2; 2]";
    description["values"]      = "[1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0]";
    description["name"] = id();
    if (subName.empty())
      return description;
    else {
      Dune::Stuff::Common::ExtendedParameterTree extendedDescription;
      extendedDescription.add(description, subName);
      return extendedDescription;
    }
  } // ... defaultSettings(...)

  static ThisType* create(const DSC::ExtendedParameterTree settings)
  {
    // get data
    const std::vector<DomainFieldType> lowerLefts  = settings.getVector("lowerLeft", DomainFieldType(0), dimDomain);
    const std::vector<DomainFieldType> upperRights = settings.getVector("upperRight", DomainFieldType(1), dimDomain);
    const std::vector<size_t> numElements          = settings.getVector("numElements", size_t(1), dimDomain);
    size_t subdomains = 1;
    for (size_t ii = 0; ii < numElements.size(); ++ii)
      subdomains *= numElements[ii];
    const std::vector<RangeFieldType> rangeFieldTypeValues =
        settings.getVector("values", RangeFieldType(1), subdomains);
    std::vector<RangeType> rangeTypeValues(rangeFieldTypeValues.size(), RangeType(0));
    for (size_t ii = 0; ii < rangeFieldTypeValues.size(); ++ii)
      rangeTypeValues[ii] = rangeFieldTypeValues[ii];
    // convert and leave the checks to the base constructor
    DomainType lowerLeft;
    DomainType upperRight;
    for (int dd = 0; dd < dimDomain; ++dd) {
      lowerLeft[dd]  = lowerLefts[dd];
      upperRight[dd] = upperRights[dd];
    }
    // create and return
    return new ThisType(lowerLeft, upperRight, numElements, rangeTypeValues);
  } // ... create(...)
};


} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTION_CHECKERBOARD_HH
