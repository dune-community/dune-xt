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
  using typename BaseType::EntityType;
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

  std::unique_ptr<LocalFunctionType> local_function(const EntityType& entity) const override final
  {
    return std::make_unique<LocalFunction>(smooth_function_storage_.access(), entity);
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
    using GeometryType = typename EntityType::Geometry;

  public:
    using typename BaseType::DomainType;
    using typename BaseType::RangeType;
    using typename BaseType::DerivativeRangeType;

    LocalFunction(const SmoothFunctionType& smooth_function, const EntityType& en)
      : BaseType(en)
      , smooth_function_(smooth_function)
      , geometry_(nullptr) // <- happens in post_bind
    {
      post_bind(en);
    }

    LocalFunction(const SmoothFunctionType& smooth_function)
      : BaseType()
      , smooth_function_(smooth_function)
      , geometry_(nullptr)
    {
    }

    using BaseType::evaluate;
    using BaseType::jacobian;
    using BaseType::derivative;

    void post_bind(const EntityType& en) override final
    {
      geometry_ = std::make_unique<GeometryType>(en.geometry());
    }

    int order(const Common::Parameter& mu = {}) const override final
    {
      if (!geometry_)
        DUNE_THROW(Exceptions::this_function_is_not_bound_to_an_entity_yet, "");
      return smooth_function_.order(mu);
    }

    RangeType evaluate(const DomainType& xx, const Common::Parameter& mu = {}) const override final
    {
      if (!geometry_)
        DUNE_THROW(Exceptions::this_function_is_not_bound_to_an_entity_yet, "");
      this->ensure_this_is_a_valid_point(xx);
      return smooth_function_.evaluate(geometry_->global(xx), mu);
    }

    DerivativeRangeType jacobian(const DomainType& xx, const Common::Parameter& mu = {}) const override final
    {
      if (!geometry_)
        DUNE_THROW(Exceptions::this_function_is_not_bound_to_an_entity_yet, "");
      this->ensure_this_is_a_valid_point(xx);
      return smooth_function_.jacobian(geometry_->global(xx), mu);
    }

    DerivativeRangeType derivative(const std::array<size_t, d>& alpha,
                                   const DomainType& xx,
                                   const Common::Parameter& mu = {}) const override final
    {
      if (!geometry_)
        DUNE_THROW(Exceptions::this_function_is_not_bound_to_an_entity_yet, "");
      this->ensure_this_is_a_valid_point(xx);
      return smooth_function_.derivative(alpha, geometry_->global(xx), mu);
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
