#ifndef DUNE_STUFF_FUNCTION_INTERFACE_hh
#define DUNE_STUFF_FUNCTION_INTERFACE_hh

// dune-common
#include <dune/common/fvector.hh>

namespace Dune {

namespace Stuff {

namespace Function {

template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class Interface
{
public:
  typedef DomainFieldImp DomainFieldType;

  static const int dimDomain = domainDim;

  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef RangeFieldImp RangeFieldType;

  static const int dimRange = rangeDim;

  typedef Dune::FieldVector<RangeFieldType, dimRange> RangeType;

  virtual void evaluate(const DomainType&, RangeType&) const = 0;
}; // class Interface

} // namespace Function

} // namespace Stuff

} // namespace Dune

#endif // DUNE_STUFF_FUNCTION_INTERFACE_hh
