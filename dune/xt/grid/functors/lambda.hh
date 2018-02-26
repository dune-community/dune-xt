// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2017 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_XT_GRID_FUNCTORS_LAMBDA_HH
#define DUNE_XT_GRID_FUNCTORS_LAMBDA_HH

#include <functional>

#include "interfaces.hh"

namespace Dune {
namespace XT {
namespace Grid {


template <class GL>
class ElementLambdaFunctor : public ElementFunctor<GL>
{
  using BaseType = ElementFunctor<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::ElementType;
  using ApplyLambdaType = std::function<void(const ElementType&)>;
  using PrepareLambdaType = std::function<void()>;
  using FinalizeLambdaType = PrepareLambdaType;

  ElementLambdaFunctor(ApplyLambdaType apply_lambda,
                       PrepareLambdaType prepare_lambda = []() {},
                       FinalizeLambdaType finalize_lambda = []() {})
    : apply_lambda_(apply_lambda)
    , prepare_lambda_(prepare_lambda)
    , finalize_lambda_(finalize_lambda)
  {
  }

  void prepare() override final
  {
    prepare_lambda_();
  }

  void apply_local(const ElementType& element) override final
  {
    apply_lambda_(element);
  }

  void finalize() override final
  {
    finalize_lambda_();
  }

  BaseType* copy() override
  {
    return new ElementLambdaFunctor<GL>(*this);
  }

private:
  const ApplyLambdaType apply_lambda_;
  const PrepareLambdaType prepare_lambda_;
  const FinalizeLambdaType finalize_lambda_;
}; // class ElementLambdaFunctor


template <class GL>
class IntersectionLambdaFunctor : public IntersectionFunctor<GL>
{
  using BaseType = IntersectionFunctor<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::ElementType;
  using typename BaseType::IntersectionType;
  using ApplyLambdaType =
      std::function<void(const IntersectionType& intersection, const ElementType&, const ElementType&)>;
  using PrepareLambdaType = std::function<void()>;
  using FinalizeLambdaType = PrepareLambdaType;

  IntersectionLambdaFunctor(ApplyLambdaType apply_lambda,
                            PrepareLambdaType prepare_lambda = []() {},
                            FinalizeLambdaType finalize_lambda = []() {})
    : apply_lambda_(apply_lambda)
    , prepare_lambda_(prepare_lambda)
    , finalize_lambda_(finalize_lambda)
  {
  }

  void prepare() override final
  {
    prepare_lambda_();
  }

  void apply_local(const IntersectionType& intersection,
                   const ElementType& inside_element,
                   const ElementType& outside_element) override final
  {
    apply_lambda_(intersection, inside_element, outside_element);
  }

  void finalize() override final
  {
    finalize_lambda_();
  }

  BaseType* copy() override
  {
    return new IntersectionLambdaFunctor(*this);
  }

private:
  const ApplyLambdaType apply_lambda_;
  const PrepareLambdaType prepare_lambda_;
  const FinalizeLambdaType finalize_lambda_;
}; // class IntersectionLambdaFunctor


template <class GL>
class ElementAndIntersectionLambdaFunctor : public ElementAndIntersectionFunctor<GL>
{
  using BaseType = ElementAndIntersectionFunctor<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::ElementType;
  using typename BaseType::IntersectionType;
  using ElementApplyLambdaType = std::function<void(const ElementType&)>;
  using IntersectionApplyLambdaType =
      std::function<void(const IntersectionType& intersection, const ElementType&, const ElementType&)>;
  using PrepareLambdaType = std::function<void()>;
  using FinalizeLambdaType = PrepareLambdaType;

  ElementAndIntersectionLambdaFunctor(ElementApplyLambdaType element_apply_lambda,
                                      IntersectionApplyLambdaType intersection_apply_lambda,
                                      PrepareLambdaType prepare_lambda = []() {},
                                      FinalizeLambdaType finalize_lambda = []() {})
    : element_apply_lambda_(element_apply_lambda)
    , intersection_apply_lambda_(intersection_apply_lambda)
    , prepare_lambda_(prepare_lambda)
    , finalize_lambda_(finalize_lambda)
  {
  }

  void prepare() override final
  {
    prepare_lambda_();
  }

  void apply_local(const ElementType& element) override final
  {
    element_apply_lambda_(element);
  }

  void apply_local(const IntersectionType& intersection,
                   const ElementType& inside_element,
                   const ElementType& outside_element) override final
  {
    intersection_apply_lambda_(intersection, inside_element, outside_element);
  }

  void finalize() override final
  {
    finalize_lambda_();
  }

  BaseType* copy() override
  {
    return new ElementAndIntersectionLambdaFunctor<GL>(*this);
  }

private:
  const ElementApplyLambdaType element_apply_lambda_;
  const IntersectionApplyLambdaType intersection_apply_lambda_;
  const PrepareLambdaType prepare_lambda_;
  const FinalizeLambdaType finalize_lambda_;
}; // class ElementAndIntersectionLambdaFunctor


///\todo Implement lambda variants of all other functor types from functors/interfaces.hh!


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_FUNCTORS_LAMBDA_HH
