// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014)
//   Rene Milk       (2014 - 2015)
//   Sven Kaulmann   (2014)
//   Tobias Leibner  (2014)

#ifndef DUNE_XT_GRID_PROVIDER_DEFAULT_HH
#define DUNE_XT_GRID_PROVIDER_DEFAULT_HH

#include <memory>

#include <dune/xt/common/memory.hh>

#include "interface.hh"

namespace Dune {
namespace XT {
namespace Grid {
namespace Providers {


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

  explicit Default(GridType& grid_in)
    : StorageProviderBaseType(grid_in)
  {
  }

  /**
   * \note Takes ownership of grid_ptr in the sense that you must not delete it manually!
   */
  explicit Default(GridType* grid_ptr)
    : StorageProviderBaseType(grid_ptr)
  {
  }

  explicit Default(std::shared_ptr<GridType> grid_ptr)
    : StorageProviderBaseType(grid_ptr)
  {
  }

  explicit Default(std::unique_ptr<GridType>&& grid_ptr)
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


} // namespace Providers
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_PROVIDER_DEFAULT_HH
