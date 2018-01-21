// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)

#ifndef DUNE_XT_FUNCTIONS_SLICED_HH
#define DUNE_XT_FUNCTIONS_SLICED_HH

#include <dune/common/typetraits.hh>

#include <dune/xt/functions/type_traits.hh>
#include <dune/xt/functions/interfaces/localizable-function.hh>

namespace Dune {
namespace XT {
namespace Functions {


template <class LF, size_t r, size_t rC = 1>
class SlicedLocalizableFunction
{
  static_assert(AlwaysFalse<LF>::value, "Not available for matrix-valued functions (yet)!");
};


/**
 * \brief Can be used to extract components of a function as a new function.
 *
 *        Consider a given function u: R^d -> R^{d+2}, i.e. modeling [density (scalar), density*velocity
 *        (d-dimensional), energy (scalar)] in the context of the euler equations and one wishes to obtain the
 *        individual components for further postprocessing or visualization, the following code would achieve just that
 *        (in the case d = 2 with u already given):
\code
auto density = XT::Functions::make_sliced_function<1>(u, {0}, "density");
auto density_times_velocity = XT::Functions::make_sliced_function<d>(u, {1, 2}, "density_times_velocity");
auto energy = XT::Functions::make_sliced_function<1>(u, {3}, "energy");
\endcode
 */
template <class LF, size_t r>
class SlicedLocalizableFunction<LF, r, 1>
    : public XT::Functions::LocalizableFunctionInterface<typename LF::E, typename LF::D, LF::d, typename LF::R, r, 1>
{
  static_assert(is_localizable_function<LF>::value, "");
  static_assert(r <= LF::r, "Does not make sense!");
  using BaseType =
      XT::Functions::LocalizableFunctionInterface<typename LF::E, typename LF::D, LF::d, typename LF::R, r, 1>;

  class SlicedLocalFunction
      : public XT::Functions::LocalfunctionInterface<typename LF::E, typename LF::D, LF::d, typename LF::R, r, 1>
  {
    using BaseType = XT::Functions::LocalfunctionInterface<typename LF::E, typename LF::D, LF::d, typename LF::R, r, 1>;

  public:
    using typename BaseType::EntityType;
    using typename BaseType::DomainType;
    using typename BaseType::RangeType;
    using typename BaseType::JacobianRangeType;

    SlicedLocalFunction(const LF& function, const std::array<size_t, r>& dims, const EntityType& ent)
      : BaseType(ent)
      , local_function_(function.local_function(ent))
      , dims_(dims)
    {
    }

    size_t order(const XT::Common::Parameter& = {}) const override final
    {
      return local_function_->order();
    }

    void evaluate(const DomainType& xx, RangeType& ret, const XT::Common::Parameter& mu = {}) const override final
    {
      const auto value = local_function_->evaluate(xx, mu);
      for (size_t ii = 0; ii < r; ++ii)
        ret[ii] = value[dims_[ii]];
    }

    void jacobian(const DomainType& /*xx*/,
                  JacobianRangeType& /*ret*/,
                  const XT::Common::Parameter& /*mu*/ = {}) const override final
    {
      DUNE_THROW(NotImplemented, "Yet!");
    }

  private:
    const std::unique_ptr<typename LF::LocalfunctionType> local_function_;
    const std::array<size_t, r>& dims_;
  }; // class SlicedLocalFunction

public:
  using typename BaseType::EntityType;
  using typename BaseType::LocalfunctionType;

  SlicedLocalizableFunction(const LF& function, const std::array<size_t, r>& dims, const std::string& nm = "")
    : function_(function)
    , dims_(dims)
    , name_(nm)
  {
    for (size_t ii = 0; ii < r; ++ii)
      if (dims_[ii] >= LF::r)
        DUNE_THROW(InvalidStateException,
                   "LF::r = " << LF::r << "\n   "
                              << "r = "
                              << r
                              << "\n   "
                              << "dims["
                              << ii
                              << "] = "
                              << dims_[ii]);
  }

  std::string name() const override final
  {
    return name_.empty() ? "sliced " + function_.name() : name_;
  }

  std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const override final
  {
    return std::make_unique<SlicedLocalFunction>(function_, dims_, entity);
  }

private:
  const LF& function_;
  const std::array<size_t, r> dims_;
  const std::string name_;
}; // class SlicedLocalizableFunction


template <size_t sliced_r, class E, class D, size_t d, class R, size_t r>
SlicedLocalizableFunction<LocalizableFunctionInterface<E, D, d, R, r>, sliced_r>
make_sliced_function(const LocalizableFunctionInterface<E, D, d, R, r>& function,
                     const std::array<size_t, sliced_r>& dims,
                     const std::string& name = "")
{
  return SlicedLocalizableFunction<LocalizableFunctionInterface<E, D, d, R, r>, sliced_r>(function, dims, name);
}


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_SLICED_HH
