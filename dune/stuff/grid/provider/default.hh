// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_GRID_PROVIDER_DEFAULT_HH
#define DUNE_STUFF_GRID_PROVIDER_DEFAULT_HH

#include <memory>

#include "interface.hh"

namespace Dune {
namespace Stuff {
namespace Grid {
namespace Providers {

#if HAVE_DUNE_GRID


template <class GridImp>
class ConstDefault : public ConstProviderInterface<GridImp>
{
  typedef ConstProviderInterface<GridImp> BaseType;

public:
  using typename BaseType::GridType;

  static const std::string static_id()
  {
    return BaseType::static_id();
  }

  ConstDefault(const GridType* grid_ptr)
    : grid_(grid_ptr)
  {
  }

  ConstDefault(std::shared_ptr<const GridType> grid_ptr)
    : grid_(grid_ptr)
  {
  }

  ConstDefault(std::unique_ptr<const GridType>&& grid_ptr)
    : grid_(grid_ptr)
  {
  }

  virtual ~ConstDefault()
  {
  }

  virtual std::shared_ptr<const GridType> grid() const DS_OVERRIDE
  {
    return grid_;
  }

private:
  std::shared_ptr<const GridType> grid_;
}; // class ConstDefault


template <class GridImp>
class Default : public ProviderInterface<GridImp>
{
  typedef ProviderInterface<GridImp> BaseType;

public:
  using typename BaseType::GridType;

  static const std::string static_id()
  {
    return BaseType::static_id();
  }

  Default(GridType* grid_ptr)
    : grid_(grid_ptr)
  {
  }

  Default(std::shared_ptr<GridType> grid_ptr)
    : grid_(grid_ptr)
  {
  }

  Default(std::unique_ptr<GridType>&& grid_ptr)
    : grid_(grid_ptr)
  {
  }

  virtual ~Default()
  {
  }

  virtual std::shared_ptr<GridType> grid() DS_OVERRIDE
  {
    return grid_;
  }

  virtual std::shared_ptr<const GridType> grid() const DS_OVERRIDE
  {
    return grid_;
  }

private:
  std::shared_ptr<GridType> grid_;
}; // class Default


#else // HAVE_DUNE_GRID


template <class GridImp>
class ConstDefault
{
  static_assert(AlwaysFalse<GridImp>::value, "You are missing dune-grid!");
};


template <class GridImp>
class Default
{
  static_assert(AlwaysFalse<GridImp>::value, "You are missing dune-grid!");
};


#endif // HAVE_DUNE_GRID

} // namespace Providers
} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_GRID_PROVIDER_DEFAULT_HH
