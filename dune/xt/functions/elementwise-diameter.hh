// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)
//   René Fritze     (2020)

#ifndef DUNE_XT_FUNCTIONS_ELEMENTWISE_DIAMETER_HH
#define DUNE_XT_FUNCTIONS_ELEMENTWISE_DIAMETER_HH

#include <dune/xt/functions/interfaces/grid-function.hh>

namespace Dune::XT::Functions {


template <class E>
class ElementwiseDiameterFunction : public GridFunctionInterface<E>
{
  using ThisType = ElementwiseDiameterFunction;
  using BaseType = GridFunctionInterface<E>;

  class LocalFunction : public XT::Functions::ElementFunctionInterface<E>
  {
    using BaseType = XT::Functions::ElementFunctionInterface<E>;

  public:
    using typename BaseType::DerivativeRangeReturnType;
    using typename BaseType::DomainType;
    using typename BaseType::ElementType;
    using typename BaseType::RangeReturnType;

    LocalFunction()
      : BaseType()
      , diameter_(0)
    {}

  protected:
    void post_bind(const ElementType& element) override final
    {
      diameter_ = Grid::diameter(element);
    }

  public:
    int order(const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      return 0;
    }

    RangeReturnType evaluate(const DomainType& /*xx*/, const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      return diameter_;
    }

    DerivativeRangeReturnType jacobian(const DomainType& /*xx*/,
                                       const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      return DerivativeRangeReturnType();
    }

  private:
    double diameter_;
  }; // class LocalFunction

public:
  using BaseType::d;
  using BaseType::r;
  using BaseType::rC;
  using typename BaseType::LocalFunctionType;

  ElementwiseDiameterFunction(const std::string& nm = "ElementwiseDiameterFunction")
    : BaseType()
    , name_(nm)
  {}

  ElementwiseDiameterFunction(const ThisType&) = default;

  ElementwiseDiameterFunction(ThisType&&) = default;


private:
  ThisType* copy_as_grid_function_impl() const override
  {
    return new ThisType(*this);
  }

public:
  std::unique_ptr<ThisType> copy_as_grid_function() const
  {
    return std::unique_ptr<ThisType>(this->copy_as_grid_function_impl());
  }
  std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    return std::make_unique<LocalFunction>();
  }

  std::string name() const override final
  {
    return name_;
  }

private:
  std::string name_;
}; // class ElementwiseDiameterFunction


} // namespace Dune::XT::Functions

#endif // DUNE_XT_FUNCTIONS_ELEMENTWISE_DIAMETER_HH
