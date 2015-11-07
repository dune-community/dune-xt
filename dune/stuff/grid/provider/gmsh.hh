// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Kirsten Weber

#ifndef DUNE_STUFF_GRID_PROVIDER_GMSH_HH
#define DUNE_STUFF_GRID_PROVIDER_GMSH_HH

#if HAVE_DUNE_GRID

#include <memory>
#include <type_traits>

#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/alugrid.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/grid/provider/interface.hh>

namespace Dune {
namespace Stuff {
namespace Grid {
namespace Providers {

/**
 * \brief   Gmsh grid provider
 */
template <class GridImp>
class Gmsh : public Grid::ProviderInterface<GridImp>
{
  static_assert(!(std::is_same<GridImp, Dune::YaspGrid<GridImp::dimension>>::value),
                "GmshReader does not work with YaspGrid!");
  static_assert(!(std::is_same<GridImp, Dune::SGrid<GridImp::dimension, GridImp::dimension>>::value),
                "GmshReader does not work with SGrid!");
  typedef Grid::ProviderInterface<GridImp> BaseType;
  typedef Gmsh<GridImp> ThisType;

public:
  using typename BaseType::GridType;

  static const std::string static_id()
  {
    return BaseType::static_id() + ".gmsh";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    std::string filename = "g.msh";
    if (std::is_same<ALUGrid<2, 2, simplex, conforming>, GridType>::value
        || std::is_same<ALUGrid<2, 2, simplex, nonconforming>, GridType>::value)
      filename = "gmsh_2d_simplices.msh";
    Common::Configuration config("filename", filename);
    if (sub_name.empty())
      return config;
    else {
      Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr<ThisType> create(const Common::Configuration config = default_config(),
                                          const std::string sub_name = static_id())
  {
    const Common::Configuration cfg         = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration default_cfg = default_config();
    return Common::make_unique<ThisType>(cfg.get("filename", default_cfg.get<std::string>("filename")));
  }

  Gmsh(const std::string filename)
  {
    grid_ = std::shared_ptr<GridType>(GmshReader<GridType>::read(filename));
  }

  Gmsh(ThisType&& source) = default;
  Gmsh(const ThisType& other) = default;

  virtual ~Gmsh() = default;

  ThisType& operator=(const ThisType& other) = default;
  ThisType& operator=(ThisType&& source) = default;

  virtual const GridType& grid() const override final
  {
    return *grid_;
  }

  virtual GridType& grid() override final
  {
    return *grid_;
  }

  const std::shared_ptr<const GridType> grid_ptr() const
  {
    return grid_;
  }

  std::shared_ptr<GridType> grid_ptr()
  {
    return grid_;
  }

private:
  std::shared_ptr<GridType> grid_;
}; // class Gmsh

} // namespace Providers
} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // HAVE_DUNE_GRID

#endif // DUNE_STUFF_GRID_PROVIDER_GMSH_HH
