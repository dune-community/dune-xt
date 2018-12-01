// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_XT_GRID_FUNCTORS_GENERIC_HH
#define DUNE_XT_GRID_FUNCTORS_GENERIC_HH

#include <functional>

#include "interfaces.hh"

namespace Dune {
namespace XT {
namespace Grid {


template <class GL>
class GenericElementFunctor : public ElementFunctor<GL>
{
  using BaseType = ElementFunctor<GL>;

public:
  using typename BaseType::ElementType;
  using typename BaseType::GridViewType;

  using GenericPrepareFunctionType = std::function<void()>;
  using GenericApplyFunctionType = std::function<void(const ElementType&)>;
  using GenericFinalizeFunctionType = std::function<void()>;

  GenericElementFunctor(GenericPrepareFunctionType prepare_func,
                        GenericApplyFunctionType apply_func,
                        GenericFinalizeFunctionType finalize_func)
    : prepare_func_(prepare_func)
    , apply_func_(apply_func)
    , finalize_func_(finalize_func)
  {}

  BaseType* copy() override final
  {
    return new GenericElementFunctor<GL>(*this);
  }

  void prepare() override final
  {
    prepare_func_();
  }

  void apply_local(const ElementType& element) override final
  {
    apply_func_(element);
  }

  void finalize() override final
  {
    finalize_func_();
  }

private:
  const GenericPrepareFunctionType prepare_func_;
  const GenericApplyFunctionType apply_func_;
  const GenericFinalizeFunctionType finalize_func_;
}; // class GenericElementFunctor


template <class GL>
class GenericIntersectionFunctor : public IntersectionFunctor<GL>
{
  using BaseType = IntersectionFunctor<GL>;

public:
  using typename BaseType::ElementType;
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  using GenericPrepareFunctionType = std::function<void()>;
  using GenericApplyFunctionType =
      std::function<void(const IntersectionType& intersection, const ElementType&, const ElementType&)>;
  using GenericFinalizeFunctionType = std::function<void()>;

  GenericIntersectionFunctor(GenericPrepareFunctionType prepare_func,
                             GenericApplyFunctionType apply_func,
                             GenericFinalizeFunctionType finalize_func)
    : prepare_func_(prepare_func)
    , apply_func_(apply_func)
    , finalize_func_(finalize_func)
  {}

  BaseType* copy() override final
  {
    return new GenericIntersectionFunctor(*this);
  }

  void prepare() override final
  {
    prepare_func_();
  }

  void apply_local(const IntersectionType& intersection,
                   const ElementType& inside_element,
                   const ElementType& outside_element) override final
  {
    apply_func_(intersection, inside_element, outside_element);
  }

  void finalize() override final
  {
    finalize_func_();
  }

private:
  const GenericPrepareFunctionType prepare_func_;
  const GenericApplyFunctionType apply_func_;
  const GenericFinalizeFunctionType finalize_func_;
}; // class GenericIntersectionFunctor


template <class GL>
class GenericElementAndIntersectionFunctor : public ElementAndIntersectionFunctor<GL>
{
  using BaseType = ElementAndIntersectionFunctor<GL>;

public:
  using typename BaseType::ElementType;
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  using GenericPrepareFunctionType = std::function<void()>;
  using GenericElementApplyFunctionType = std::function<void(const ElementType&)>;
  using GenericIntersectionApplyFunctionType =
      std::function<void(const IntersectionType& intersection, const ElementType&, const ElementType&)>;
  using GenericFinalizeFunctionType = std::function<void()>;

  GenericElementAndIntersectionFunctor(GenericPrepareFunctionType prepare_func,
                                       GenericElementApplyFunctionType element_apply_func,
                                       GenericIntersectionApplyFunctionType intersection_apply_func,
                                       GenericFinalizeFunctionType finalize_func)
    : prepare_func_(prepare_func)
    , element_apply_func_(element_apply_func)
    , intersection_apply_func_(intersection_apply_func)
    , finalize_func_(finalize_func)
  {}

  void prepare() override final
  {
    prepare_func_();
  }

  void apply_local(const ElementType& element) override final
  {
    element_apply_func_(element);
  }

  void apply_local(const IntersectionType& intersection,
                   const ElementType& inside_element,
                   const ElementType& outside_element) override final
  {
    intersection_apply_func_(intersection, inside_element, outside_element);
  }

  void finalize() override final
  {
    finalize_func_();
  }

  BaseType* copy() override
  {
    return new GenericElementAndIntersectionFunctor<GL>(*this);
  }

private:
  const GenericPrepareFunctionType prepare_func_;
  const GenericElementApplyFunctionType element_apply_func_;
  const GenericIntersectionApplyFunctionType intersection_apply_func_;
  const GenericFinalizeFunctionType finalize_func_;
}; // class GenericElementAndIntersectionFunctor


///\todo Implement lambda variants of all other functor types from functors/interfaces.hh!


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_FUNCTORS_GENERIC_HH
