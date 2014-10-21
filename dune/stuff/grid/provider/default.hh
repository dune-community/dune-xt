// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_GRID_PROVIDER_DEFAULT_HH
#define DUNE_STUFF_GRID_PROVIDER_DEFAULT_HH

#include <memory>

#include <dune/stuff/common/memory.hh>

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

  ConstDefault(const GridType& grid_in)
    : grid_(grid_in)
  {
  }

  virtual ~ConstDefault()
  {
  }

  virtual const GridType& grid() const override
  {
    return grid_;
  }

private:
  const GridType& grid_;
}; // class ConstDefault


template <class GridImp>
class Default : Common::StorageProvider<GridImp>, public ProviderInterface<GridImp>
{
  typedef Common::StorageProvider<GridImp> StorageProviderBaseType;
  typedef ProviderInterface<GridImp> GridProviderBaseType;

public:
  using typename GridProviderBaseType::GridType;

  static const std::string static_id()
  {
    return GridProviderBaseType::static_id();
  }

  Default(GridType& grid_in)
    : StorageProviderBaseType(grid_in)
  {
  }

  /**
   * \note Takes ownership of grid_ptr in the sense that you must not delete it manually!
   */
  Default(GridType* grid_ptr)
    : StorageProviderBaseType(grid_ptr)
  {
  }

  Default(std::shared_ptr<GridType> grid_ptr)
    : StorageProviderBaseType(grid_ptr)
  {
  }

  Default(std::unique_ptr<GridType>&& grid_ptr)
    : StorageProviderBaseType(grid_ptr)
  {
  }

  virtual ~Default()
  {
  }

  virtual GridType& grid() override
  {
    return this->storage_access();
  }

  virtual const GridType& grid() const override
  {
    return this->storage_access();
  }
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
