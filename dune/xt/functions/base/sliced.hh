// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   René Fritze     (2018 - 2019)
//   Tim Keil        (2018)
//   Tobias Leibner  (2020)

#ifndef DUNE_XT_FUNCTIONS_BASE_SLICED_HH
#define DUNE_XT_FUNCTIONS_BASE_SLICED_HH

#include <dune/common/typetraits.hh>

#include <dune/xt/functions/type_traits.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>


namespace Dune {
namespace XT {
namespace Functions {


template <class LF, size_t r, size_t rC = 1>
class SlicedGridFunction
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
auto density                = XT::Functions::make_sliced_function<1>(u, {0},    "density");
auto density_times_velocity = XT::Functions::make_sliced_function<d>(u, {1, 2}, "density_times_velocity");
auto energy                 = XT::Functions::make_sliced_function<1>(u, {3},    "energy");
\endcode
 */
template <class GF, size_t r>
class SlicedGridFunction<GF, r, 1> : public XT::Functions::GridFunctionInterface<typename GF::E, r, 1, typename GF::R>
{
  static_assert(is_grid_function<GF>::value, "");
  static_assert(r <= GF::r, "Does not make sense!");

  using ThisType = SlicedGridFunction;
  using BaseType = XT::Functions::GridFunctionInterface<typename GF::E, r, 1, typename GF::R>;

  class SlicedLocalFunction : public XT::Functions::ElementFunctionInterface<typename GF::E, r, 1, typename GF::R>
  {
    using BaseType = XT::Functions::ElementFunctionInterface<typename GF::E, r, 1, typename GF::R>;

  public:
    using typename BaseType::DerivativeRangeReturnType;
    using typename BaseType::DerivativeRangeType;
    using typename BaseType::DomainType;
    using typename BaseType::ElementType;
    using typename BaseType::RangeReturnType;
    using typename BaseType::RangeType;

    SlicedLocalFunction(const GF& function, const std::array<size_t, r>& dims)
      : BaseType(function.parameter_type())
      , function_(function.copy_as_grid_function())
      , local_function_(function_->local_function())
      , dims_(dims)
    {}

  protected:
    void post_bind(const ElementType& element) override final
    {
      local_function_->bind(element);
    }

  public:
    int order(const XT::Common::Parameter& = {}) const override final
    {
      return local_function_->order();
    }

    RangeReturnType evaluate(const DomainType& xx, const XT::Common::Parameter& param = {}) const override final
    {
      RangeReturnType ret;
      const auto value = local_function_->evaluate(xx, param);
      for (size_t ii = 0; ii < r; ++ii)
        ret[ii] = value[dims_[ii]];
      return ret;
    }

    DerivativeRangeReturnType jacobian(const DomainType& /*xx*/,
                                       const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      DUNE_THROW(NotImplemented, "Yet!");
    }

  private:
    std::unique_ptr<GridFunctionInterface<typename GF::E, GF::r, GF::rC, typename GF::R>> function_;
    std::unique_ptr<typename GF::LocalFunctionType> local_function_;
    const std::array<size_t, r> dims_;
  }; // class SlicedLocalFunction

public:
  using typename BaseType::E;
  using typename BaseType::LocalFunctionType;

  SlicedGridFunction(const GF& function, const std::array<size_t, r>& dims, const std::string& nm = "")
    : BaseType(function.parameter_type())
    , function_(function.copy_as_grid_function())
    , dims_(dims)
    , name_(nm)
  {
    for (size_t ii = 0; ii < r; ++ii)
      if (dims_[ii] >= GF::r)
        DUNE_THROW(InvalidStateException,
                   "GF::r = " << GF::r << "\n   "
                              << "r = " << r << "\n   "
                              << "dims[" << ii << "] = " << dims_[ii]);
  }

  SlicedGridFunction(const ThisType& other)
    : BaseType(other)
    , function_(other.function_->copy_as_grid_function())
    , dims_(other.dims_)
    , name_(other.name_)
  {}

  SlicedGridFunction(ThisType&&) = default;


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
  std::string name() const override final
  {
    return name_.empty() ? "sliced " + function_->name() : name_;
  }

  std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    return std::make_unique<SlicedLocalFunction>(*function_, dims_);
  }

private:
  std::unique_ptr<GridFunctionInterface<typename GF::E, GF::r, GF::rC, typename GF::R>> function_;
  const std::array<size_t, r> dims_;
  const std::string& name_;
}; // class SlicedGridFunction


template <size_t sliced_r, class E, size_t r>
SlicedGridFunction<GridFunctionInterface<E, r>, sliced_r> make_sliced_function(
    const GridFunctionInterface<E, r>& function, const std::array<size_t, sliced_r>& dims, const std::string& name = "")
{
  return SlicedGridFunction<GridFunctionInterface<E, r>, sliced_r>(function, dims, name);
}


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_SLICED_HH
