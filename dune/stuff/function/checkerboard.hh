#ifndef DUNE_STUFF_FUNCTION_CHECKERBOARD_HH
#define DUNE_STUFF_FUNCTION_CHECKERBOARD_HH

#include <vector>

#include "interface.hh"

namespace Dune {
namespace Stuff {
namespace Function {

template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class Checkerboard : public Interface<DomainFieldImp, domainDim, RangeFieldImp, rangeDim>
{
public:
  typedef DomainFieldImp DomainFieldType;

  static const int dimDomain = domainDim;

  typedef RangeFieldImp RangeFieldType;

  static const int dimRange = rangeDim;

  typedef Interface<DomainFieldType, dimDomain, RangeFieldType, dimRange> BaseType;

  typedef Checkerboard<DomainFieldType, dimDomain, RangeFieldType, dimRange> ThisType;

  typedef typename BaseType::DomainType DomainType;

  typedef typename BaseType::RangeType RangeType;

  Checkerboard(const DomainType lowerLeft, const DomainType upperRight, const std::vector<unsigned int> numElements,
               const std::vector<RangeFieldType> components)
    : lowerLeft_(lowerLeft)
    , upperRight_(upperRight)
    , numElements_(numElements)
    , components_(components)
  {
    // get total number of subdomains
    unsigned int totalSubdomains = 1;
    for (int d = 0; d < dimDomain; ++d) {
      totalSubdomains *= numElements_[d];
    }
    assert(totalSubdomains <= components_.size() && "Please provide at least as many components as subdomains!");
  }

  Checkerboard(const ThisType& other)
    : lowerLeft_(other.lowerLeft_)
    , upperRight_(other.upperRight_)
    , numElements_(other.numElements_)
    , components_(other.components_)
  {
    // no cheks necessary, since they have been carried out in the constructor of other
  }

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      lowerLeft_   = other.lowerLeft();
      upperRight_  = other.upperRight();
      numElements_ = other.numElements();
      components_  = other.components();
    }
    return this;
  }

  const DomainType& lowerLeft() const
  {
    return lowerLeft_;
  }

  const DomainType& upperRight() const
  {
    return upperRight_;
  }

  const std::vector<unsigned int>& numElements() const
  {
    return numElements_();
  }

  const std::vector<RangeFieldType>& components() const
  {
    return components_();
  }

  virtual void evaluate(const DomainType& x, RangeType& ret) const
  {
    // decide on the subdomain the point x belongs to
    std::vector<unsigned int> whichPartition;
    for (int d = 0; d < dimDomain; ++d) {
      whichPartition.push_back(
          std::floor(numElements_[d] * ((x[d] - lowerLeft_[d]) / (upperRight_[d] - lowerLeft_[d]))));
    }
    unsigned int subdomain = 0;
    if (dimDomain == 1)
      subdomain = whichPartition[0];
    else if (dimDomain == 2)
      subdomain = whichPartition[0] + whichPartition[1] * numElements_[0];
    else if (dimDomain == 3)
      subdomain = whichPartition[0] + whichPartition[1] * numElements_[0]
                  + whichPartition[2] * numElements_[1] * numElements_[0];
    else {
      DUNE_THROW(Dune::NotImplemented, "\nError: not implemented for grid dimensions other than 1, 2 or 3!");
    } // decide on the subdomain the point x belongs to
    // return the component that belongs to the subdomain of x
    ret = components_[subdomain];
  } // virtual void evaluate(const DomainType& x, RangeType& ret) const

private:
  DomainType lowerLeft_;
  DomainType upperRight_;
  std::vector<unsigned int> numElements_;
  std::vector<RangeFieldType> components_;
}; // class Checkerboard

} // namespace Function
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTION_CHECKERBOARD_HH
