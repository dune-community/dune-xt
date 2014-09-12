// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Sven Kaulmann

#ifndef DUNE_STUFF_FUNCTION_GLOBAL_HH
#define DUNE_STUFF_FUNCTION_GLOBAL_HH

#include <functional>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/common/memory.hh>

namespace Dune {
namespace Stuff {


/**
 * Global-valued function you can pass a lambda expression to that gets evaluated
 * \example LambdaType lambda([](DomainType x) { return x;}, 1 );
 */
template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class GlobalLambdaFunction
    : public GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
{
  typedef GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> BaseType;

public:
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

private:
  typedef std::function<RangeType(DomainType)> LambdaType;

public:
  GlobalLambdaFunction(LambdaType lambda, size_t order)
    : lambda_(lambda)
    , order_(order)
  {
  }

  virtual size_t order() const DS_OVERRIDE DS_FINAL
  {
    return order_;
  }

  virtual void evaluate(const DomainType& xx, RangeType& ret) const DS_OVERRIDE DS_FINAL
  {
    ret = lambda_(xx);
  }

  virtual RangeType evaluate(const DomainType& xx) const DS_OVERRIDE DS_FINAL
  {
    return lambda_(xx);
  }

private:
  const LambdaType lambda_;
  const size_t order_;
};


} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTION_GLOBAL_HH
