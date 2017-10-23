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
class EntityLambdaFunctor : public EntityFunctor<GL>
{
  using BaseType = EntityFunctor<GL>;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;
  using ApplyLambdaType = std::function<void(const EntityType&)>;
  using PrepareLambdaType = std::function<void()>;
  using FinalizeLambdaType = PrepareLambdaType;

  EntityLambdaFunctor(ApplyLambdaType apply_lambda,
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

  void apply_local(const EntityType& entity) override final
  {
    apply_lambda_(entity);
  }

  void finalize() override final
  {
    finalize_lambda_();
  }

private:
  const ApplyLambdaType apply_lambda_;
  const PrepareLambdaType prepare_lambda_;
  const FinalizeLambdaType finalize_lambda_;
}; // class EntityLambdaFunctor


template <class GL, class R>
class EntityReturnLambdaFunctor : public EntityReturnFunctor<GL, R>
{
  using BaseType = EntityReturnFunctor<GL, R>;

public:
  using typename BaseType::ResultType;
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;
  using ComputeLambdaType = std::function<ResultType(const EntityType&)>;
  using PrepareLambdaType = std::function<void()>;
  using FinalizeLambdaType = PrepareLambdaType;

  EntityReturnLambdaFunctor(ComputeLambdaType compute_lambda,
                            PrepareLambdaType prepare_lambda = []() {},
                            FinalizeLambdaType finalize_lambda = []() {})
    : compute_lambda_(compute_lambda)
    , prepare_lambda_(prepare_lambda)
    , finalize_lambda_(finalize_lambda)
    , result_(0)
  {
  }

  EntityReturnFunctor<GridLayerType, ResultType>* copy() const override final
  {
    return new EntityReturnLambdaFunctor<GridLayerType, ResultType>(compute_lambda_, prepare_lambda_, finalize_lambda_);
  }

  void prepare() override final
  {
    prepare_lambda_();
  }

  void apply_local(const EntityType& entity) override final
  {
    result_ += compute_locally(entity);
  }

  ResultType compute_locally(const EntityType& entity) override final
  {
    return compute_lambda_(entity);
  }

  void finalize() override final
  {
    finalize_lambda_();
  }

  ResultType result() const override final
  {
    return result_;
  }

private:
  const ComputeLambdaType compute_lambda_;
  const PrepareLambdaType prepare_lambda_;
  const FinalizeLambdaType finalize_lambda_;
  ResultType result_;
}; // class EntityReturnLambdaFunctor


template <class GL>
class IntersectionLambdaFunctor : public IntersectionFunctor<GL>
{
  using BaseType = IntersectionFunctor<GL>;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;
  using typename BaseType::IntersectionType;
  using ApplyLambdaType =
      std::function<void(const IntersectionType& intersection, const EntityType&, const EntityType&)>;
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
                   const EntityType& inside_entity,
                   const EntityType& outside_entity) override final
  {
    apply_lambda_(intersection, inside_entity, outside_entity);
  }

  void finalize() override final
  {
    finalize_lambda_();
  }

private:
  const ApplyLambdaType apply_lambda_;
  const PrepareLambdaType prepare_lambda_;
  const FinalizeLambdaType finalize_lambda_;
}; // class IntersectionLambdaFunctor


template <class GL>
class EntityAndIntersectionLambdaFunctor : public EntityAndIntersectionFunctor<GL>
{
  using BaseType = EntityAndIntersectionFunctor<GL>;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;
  using typename BaseType::IntersectionType;
  using EntityApplyLambdaType = std::function<void(const EntityType&)>;
  using IntersectionApplyLambdaType =
      std::function<void(const IntersectionType& intersection, const EntityType&, const EntityType&)>;
  using PrepareLambdaType = std::function<void()>;
  using FinalizeLambdaType = PrepareLambdaType;

  EntityAndIntersectionLambdaFunctor(EntityApplyLambdaType entity_apply_lambda,
                                     IntersectionApplyLambdaType intersection_apply_lambda,
                                     PrepareLambdaType prepare_lambda = []() {},
                                     FinalizeLambdaType finalize_lambda = []() {})
    : entity_apply_lambda_(entity_apply_lambda)
    , intersection_apply_lambda_(intersection_apply_lambda)
    , prepare_lambda_(prepare_lambda)
    , finalize_lambda_(finalize_lambda)
  {
  }

  void prepare() override final
  {
    prepare_lambda_();
  }

  void apply_local(const EntityType& entity) override final
  {
    entity_apply_lambda_(entity);
  }

  void apply_local(const IntersectionType& intersection,
                   const EntityType& inside_entity,
                   const EntityType& outside_entity) override final
  {
    intersection_apply_lambda_(intersection, inside_entity, outside_entity);
  }

  void finalize() override final
  {
    finalize_lambda_();
  }

private:
  const EntityApplyLambdaType entity_apply_lambda_;
  const IntersectionApplyLambdaType intersection_apply_lambda_;
  const PrepareLambdaType prepare_lambda_;
  const FinalizeLambdaType finalize_lambda_;
}; // class EntityAndIntersectionLambdaFunctor


///\todo Implement lambda variants of all other functor types from functors/interfaces.hh!


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_FUNCTORS_LAMBDA_HH
