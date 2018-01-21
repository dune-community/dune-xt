// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014, 2016 - 2017)
//   Rene Milk       (2014 - 2016, 2018)
//   Tobias Leibner  (2014)

#ifndef DUNE_XT_GRID_WALKER_WRAPPER_HH
#define DUNE_XT_GRID_WALKER_WRAPPER_HH

#include "apply-on.hh"
#include "functors.hh"

namespace Dune {
namespace XT {
namespace Grid {
namespace internal {

template <class GridLayerType>
class Codim0Object : public Functor::Codim0<GridLayerType>
{
  typedef Functor::Codim0<GridLayerType> BaseType;

public:
  typedef typename BaseType::EntityType EntityType;

  virtual ~Codim0Object()
  {
  }

  virtual bool apply_on(const GridLayerType& grid_layer, const EntityType& entity) const = 0;
};

template <class GridLayerImp, class ReturnType>
class Codim0ReturnObject : public Functor::Codim0Return<GridLayerImp, ReturnType>
{
  typedef Functor::Codim0Return<GridLayerImp, ReturnType> BaseType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;

  virtual ~Codim0ReturnObject() = default;

  virtual bool apply_on(const GridLayerType& grid_layer, const EntityType& entity) const = 0;
};

template <class GridLayerType, class Codim0FunctorType>
class Codim0FunctorWrapper : public Codim0Object<GridLayerType>
{
  typedef Codim0Object<GridLayerType> BaseType;

public:
  typedef typename BaseType::EntityType EntityType;

  Codim0FunctorWrapper(Codim0FunctorType& wrapped_functor, const ApplyOn::WhichEntity<GridLayerType>* where)
    : wrapped_functor_(wrapped_functor)
    , where_(where)
  {
  }

  virtual ~Codim0FunctorWrapper()
  {
  }

  virtual void prepare() override final
  {
    wrapped_functor_.prepare();
  }

  virtual bool apply_on(const GridLayerType& grid_layer, const EntityType& entity) const override final
  {
    return where_->apply_on(grid_layer, entity);
  }

  virtual void apply_local(const EntityType& entity) override final
  {
    wrapped_functor_.apply_local(entity);
  }

  virtual void finalize() override final
  {
    wrapped_functor_.finalize();
  }

private:
  Codim0FunctorType& wrapped_functor_;
  std::unique_ptr<const ApplyOn::WhichEntity<GridLayerType>> where_;
}; // class Codim0FunctorWrapper

template <class GridLayerType>
class Codim1Object : public Functor::Codim1<GridLayerType>
{
  typedef Functor::Codim1<GridLayerType> BaseType;

public:
  typedef typename BaseType::IntersectionType IntersectionType;

  virtual ~Codim1Object()
  {
  }

  virtual bool apply_on(const GridLayerType& grid_layer, const IntersectionType& intersection) const = 0;
};

template <class GridLayerType, class Codim1FunctorType>
class Codim1FunctorWrapper : public Codim1Object<GridLayerType>
{
  typedef Codim1Object<GridLayerType> BaseType;

public:
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::IntersectionType IntersectionType;

  Codim1FunctorWrapper(Codim1FunctorType& wrapped_functor, const ApplyOn::WhichIntersection<GridLayerType>* where)
    : wrapped_functor_(wrapped_functor)
    , where_(where)
  {
  }

  virtual void prepare() override final
  {
    wrapped_functor_.prepare();
  }

  virtual bool apply_on(const GridLayerType& grid_layer, const IntersectionType& intersection) const override final
  {
    return where_->apply_on(grid_layer, intersection);
  }

  virtual void apply_local(const IntersectionType& intersection,
                           const EntityType& inside_entity,
                           const EntityType& outside_entity) override final
  {
    wrapped_functor_.apply_local(intersection, inside_entity, outside_entity);
  }

  virtual void finalize() override final
  {
    wrapped_functor_.finalize();
  }

private:
  Codim1FunctorType& wrapped_functor_;
  std::unique_ptr<const ApplyOn::WhichIntersection<GridLayerType>> where_;
}; // class Codim1FunctorWrapper

template <class GridLayerType, class WalkerType>
class WalkerWrapper : public Codim0Object<GridLayerType>, public Codim1Object<GridLayerType>
{
public:
  typedef typename Codim1Object<GridLayerType>::EntityType EntityType;
  typedef typename Codim1Object<GridLayerType>::IntersectionType IntersectionType;

  WalkerWrapper(WalkerType& grid_walker, const ApplyOn::WhichEntity<GridLayerType>* which_entities)
    : grid_walker_(grid_walker)
    , which_entities_(which_entities)
    , which_intersections_(new ApplyOn::AllIntersections<GridLayerType>())
  {
  }

  WalkerWrapper(WalkerType& grid_walker, const ApplyOn::WhichIntersection<GridLayerType>* which_intersections)
    : grid_walker_(grid_walker)
    , which_entities_(new ApplyOn::AllEntities<GridLayerType>())
    , which_intersections_(which_intersections)
  {
  }

  virtual ~WalkerWrapper()
  {
  }

  virtual void prepare() override final
  {
    grid_walker_.prepare();
  }

  virtual bool apply_on(const GridLayerType& grid_layer, const EntityType& entity) const override final
  {
    return which_entities_->apply_on(grid_layer, entity) && grid_walker_.apply_on(entity);
  }

  virtual bool apply_on(const GridLayerType& grid_layer, const IntersectionType& intersection) const override final
  {
    return which_intersections_->apply_on(grid_layer, intersection) && grid_walker_.apply_on(intersection);
  }

  virtual void apply_local(const EntityType& entity) override final
  {
    grid_walker_.apply_local(entity);
  }

  virtual void apply_local(const IntersectionType& intersection,
                           const EntityType& inside_entity,
                           const EntityType& outside_entity) override final
  {
    grid_walker_.apply_local(intersection, inside_entity, outside_entity);
  }

  virtual void finalize() override final
  {
    grid_walker_.finalize();
  }

private:
  WalkerType& grid_walker_;
  std::unique_ptr<const ApplyOn::WhichEntity<GridLayerType>> which_entities_;
  std::unique_ptr<const ApplyOn::WhichIntersection<GridLayerType>> which_intersections_;
}; // class WalkerWrapper

template <class GridLayerType>
class Codim0LambdaWrapper : public Codim0Object<GridLayerType>
{
  typedef Codim0Object<GridLayerType> BaseType;

public:
  typedef typename BaseType::EntityType EntityType;
  typedef std::function<void(const EntityType&)> LambdaType;

  Codim0LambdaWrapper(LambdaType lambda, const ApplyOn::WhichEntity<GridLayerType>* where)
    : lambda_(lambda)
    , where_(where)
  {
  }

  virtual ~Codim0LambdaWrapper()
  {
  }

  virtual bool apply_on(const GridLayerType& grid_layer, const EntityType& entity) const override final
  {
    return where_->apply_on(grid_layer, entity);
  }

  virtual void apply_local(const EntityType& entity) override final
  {
    lambda_(entity);
  }

private:
  LambdaType lambda_;
  std::unique_ptr<const ApplyOn::WhichEntity<GridLayerType>> where_;
}; // class Codim0LambdaWrapper

template <class GridLayerType>
class Codim1LambdaWrapper : public Codim1Object<GridLayerType>
{
  typedef Codim1Object<GridLayerType> BaseType;

public:
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::IntersectionType IntersectionType;
  typedef std::function<void(const IntersectionType&, const EntityType&, const EntityType&)> LambdaType;

  Codim1LambdaWrapper(LambdaType lambda, const ApplyOn::WhichIntersection<GridLayerType>* where)
    : lambda_(lambda)
    , where_(where)
  {
  }

  virtual bool apply_on(const GridLayerType& grid_layer, const IntersectionType& intersection) const override final
  {
    return where_->apply_on(grid_layer, intersection);
  }

  virtual void apply_local(const IntersectionType& intersection,
                           const EntityType& inside_entity,
                           const EntityType& outside_entity) override final
  {
    lambda_(intersection, inside_entity, outside_entity);
  }

private:
  LambdaType lambda_;
  std::unique_ptr<const ApplyOn::WhichIntersection<GridLayerType>> where_;
}; // class Codim1FunctorWrapper

} // namespace internal
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_WALKER_WRAPPER_HH
