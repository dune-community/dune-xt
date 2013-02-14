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
namespace Function {


template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class Checkerboard;


template <class DomainFieldImp, int domainDim, class RangeFieldImp>
class Checkerboard<DomainFieldImp, domainDim, RangeFieldImp, 1>
    : public Interface<DomainFieldImp, domainDim, RangeFieldImp, 1>
{
public:
  typedef Checkerboard<DomainFieldImp, domainDim, RangeFieldImp, 1> ThisType;
  typedef Interface<DomainFieldImp, domainDim, RangeFieldImp, 1> BaseType;

  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const int dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const int dimRange = BaseType::dimRange;
  typedef typename BaseType::RangeType RangeType;

  static const std::string id()
  {
    return "function.checkerboard";
  }

  Checkerboard(const DomainType _lowerLeft, const DomainType _upperRight, const std::vector<size_t> _numElements,
               const std::vector<RangeFieldType> _values, const std::string _name = id())
    : lowerLeft_(_lowerLeft)
    , upperRight_(_upperRight)
    , numElements_(_numElements)
    , values_(_values)
    , name_(_name)
  {
    // checks
    dune_static_assert((dimDomain > 0), "Really?");
    dune_static_assert((dimDomain <= 3), "Not implemented!");
    dune_static_assert((dimRange > 0), "Really?");
    // get total number of subdomains
    size_t totalSubdomains = 1;
    for (int dd = 0; dd < dimDomain; ++dd) {
      totalSubdomains *= numElements_[dd];
    }
    if (values_.size() < totalSubdomains)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                      << " please provide at least as many '_values' as subdomains given by '_numElements'!");
  }

  Checkerboard(const ThisType& _other)
    : lowerLeft_(_other.lowerLeft_)
    , upperRight_(_other.upperRight_)
    , numElements_(_other.numElements_)
    , values_(_other.values_)
    , name_(_other.name_)
  {
  }

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      lowerLeft_   = other.lowerLeft();
      upperRight_  = other.upperRight();
      numElements_ = other.numElements();
      values_      = other.values();
      name_        = other.name();
    }
    return this;
  }

  static Dune::ParameterTree createSampleDescription(const std::string subName = "")
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
  }

  static ThisType* createFromDescription(const DSC::ExtendedParameterTree description)
  {
    // get data
    const std::vector<DomainFieldType> lowerLefts  = description.getVector("lowerLeft", DomainFieldType(0), dimDomain);
    const std::vector<DomainFieldType> upperRights = description.getVector("upperRight", DomainFieldType(1), dimDomain);
    const std::vector<size_t> numElements          = description.getVector("numElements", size_t(1), dimDomain);
    size_t subdomains = 1;
    for (size_t ii = 0; ii < numElements.size(); ++ii)
      subdomains *= numElements[ii];
    const std::vector<RangeFieldType> values = description.getVector("values", RangeFieldType(0), subdomains);
    // convert and leave the checks to the constructor
    DomainType lowerLeft;
    DomainType upperRight;
    for (int dd = 0; dd < dimDomain; ++dd) {
      lowerLeft[dd]  = lowerLefts[dd];
      upperRight[dd] = upperRights[dd];
    }
    // create and return
    return new ThisType(lowerLeft, upperRight, numElements, values);
  } // static ThisType createFromParamTree(const Dune::ParameterTree paramTree)


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

  const std::vector<RangeFieldType>& values() const
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

  virtual void evaluate(const DomainType& x, RangeType& ret) const
  {
    // decide on the subdomain the point x belongs to
    std::vector<size_t> whichPartition;
    for (int d = 0; d < dimDomain; ++d) {
      whichPartition.push_back(
          std::floor(numElements_[d] * ((x[d] - lowerLeft_[d]) / (upperRight_[d] - lowerLeft_[d]))));
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
  std::vector<RangeFieldType> values_;
  std::string name_;
}; // class Checkerboard

} // namespace Function
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTION_CHECKERBOARD_HH
