// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Tobias Leibner  (2017)

#ifndef DUNE_XT_FUNCTIONS_BASE_SMOOTH_LOCALIZABLE_FUNCTION_HH
#define DUNE_XT_FUNCTIONS_BASE_SMOOTH_LOCALIZABLE_FUNCTION_HH

#include <dune/xt/common/memory.hh>

#include <dune/xt/functions/exceptions.hh>
#include <dune/xt/functions/interfaces/localizable-function.hh>
#include <dune/xt/functions/interfaces/smooth-function.hh>
#include <dune/xt/functions/type_traits.hh>

namespace Dune {
namespace XT {
namespace Functions {


template <class E_, size_t r_, size_t rC_, class R_>
class SmoothFunctionAsLocalizableWrapper : public LocalizableFunctionInterface<E_, r_, rC_, R_>
{
  using BaseType = LocalizableFunctionInterface<E_, r_, rC_, R_>;
  using ThisType = SmoothFunctionAsLocalizableWrapper<E_, r_, rC_, R_>;

public:
  using BaseType::d;
  using BaseType::r;
  using BaseType::rC;
  using typename BaseType::R;
  using typename BaseType::ElementType;
  using typename BaseType::LocalFunctionType;
  using SmoothFunctionType = SmoothFunctionInterface<d, r, rC, R>;

  SmoothFunctionAsLocalizableWrapper(const SmoothFunctionType& smooth_function)
    : smooth_function_storage_(smooth_function)
  {
  }

  SmoothFunctionAsLocalizableWrapper(const SmoothFunctionType*&& smooth_function_ptr)
    : smooth_function_storage_(std::move(smooth_function_ptr))
  {
  }

  SmoothFunctionAsLocalizableWrapper(std::unique_ptr<SmoothFunctionType>&& smooth_function_ptr)
    : smooth_function_storage_(std::move(smooth_function_ptr))
  {
  }

  /**
   * \name ´´These methods are required by LocalizableFunctionInterface.''
   * \{
   **/

  std::unique_ptr<LocalFunctionType> local_function(const ElementType& element) const override final
  {
    return std::make_unique<LocalFunction>(smooth_function_storage_.access(), element);
  }

  std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    return std::make_unique<LocalFunction>(smooth_function_storage_.access());
  }

  /**
   * \}
   * \name ´´These methods are optionally required by LocalizableFunctionInterface.''
   * \{
   **/

  std::string type() const override final
  {
    return smooth_function_storage_.access().type();
  }

  std::string name() const override final
  {
    return smooth_function_storage_.access().name();
  }

  /**
   * \}
   **/

private:
  class LocalFunction : public LocalFunctionType
  {
    using BaseType = LocalFunctionType;
    using GeometryType = typename ElementType::Geometry;

  public:
    using typename BaseType::DomainType;
    using typename BaseType::RangeType;
    using typename BaseType::DerivativeRangeType;

    LocalFunction(const SmoothFunctionType& smooth_function, const ElementType& el)
      : BaseType(el)
      , smooth_function_(smooth_function)
      , geometry_(nullptr) // <- happens in post_bind
    {
      post_bind(el);
    }

    LocalFunction(const SmoothFunctionType& smooth_function)
      : BaseType()
      , smooth_function_(smooth_function)
      , geometry_(nullptr) // <- happens in post_bind
    {
    }

    using BaseType::evaluate;
    using BaseType::jacobian;
    using BaseType::derivative;

    void post_bind(const ElementType& el) override final
    {
      geometry_ = std::make_unique<GeometryType>(el.geometry());
    }

    int order(const Common::Parameter& param = {}) const override final
    {
      if (!geometry_)
        DUNE_THROW(Exceptions::this_function_is_not_bound_to_an_element_yet, "");
      return smooth_function_.order(param);
    }

    RangeType evaluate(const DomainType& point_in_reference_element,
                       const Common::Parameter& param = {}) const override final
    {
      if (!geometry_)
        DUNE_THROW(Exceptions::this_function_is_not_bound_to_an_element_yet, "");
      this->assert_inside_reference_element(point_in_reference_element);
      return smooth_function_.evaluate(geometry_->global(point_in_reference_element), param);
    }

    DerivativeRangeType jacobian(const DomainType& point_in_reference_element,
                                 const Common::Parameter& param = {}) const override final
    {
      if (!geometry_)
        DUNE_THROW(Exceptions::this_function_is_not_bound_to_an_element_yet, "");
      this->assert_inside_reference_element(point_in_reference_element);
      return smooth_function_.jacobian(geometry_->global(point_in_reference_element), param);
    }

    DerivativeRangeType derivative(const std::array<size_t, d>& alpha,
                                   const DomainType& point_in_reference_element,
                                   const Common::Parameter& param = {}) const override final
    {
      if (!geometry_)
        DUNE_THROW(Exceptions::this_function_is_not_bound_to_an_element_yet, "");
      this->assert_inside_reference_element(point_in_reference_element);
      return smooth_function_.derivative(alpha, geometry_->global(point_in_reference_element), param);
    }

  private:
    const SmoothFunctionType& smooth_function_;
    std::unique_ptr<GeometryType> geometry_;
  }; // class LocalFunction

  const XT::Common::ConstStorageProvider<SmoothFunctionType> smooth_function_storage_;
}; // class SmoothFunctionAsLocalizableWrapper


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_BASE_SMOOTH_LOCALIZABLE_FUNCTION_HH
