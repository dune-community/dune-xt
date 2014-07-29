// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_GRIDS_PROVIDER_CUBE_HH
#define DUNE_STUFF_GRIDS_PROVIDER_CUBE_HH

#include <memory>
#include <sstream>
#include <type_traits>
#include <vector>
#include <array>

#include <boost/assign/list_of.hpp>

#include <dune/stuff/common/disable_warnings.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/sgrid.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/grid/yaspgrid.hh>
#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/configtree.hh>
#include <dune/stuff/common/memory.hh>

#include "default.hh"

namespace Dune {
namespace Stuff {
namespace Grid {
namespace Providers {
namespace internal {


template <typename GridType>
struct ElementVariant;


template <typename GridType>
struct ElementVariant
{
  static const int id = 2;
};


template <int dim>
struct ElementVariant<Dune::YaspGrid<dim>>
{
  static const int id = 1;
};


template <int dimGrid, int dimWorld>
struct ElementVariant<Dune::SGrid<dimGrid, dimWorld>>
{
  static const int id = 1;
};


#if HAVE_ALUGRID


template <int dimGrid, int dimWorld>
struct ElementVariant<Dune::ALUCubeGrid<dimGrid, dimWorld>>
{
  static const int id = 1;
};


template <int dimGrid, int dimWorld>
struct ElementVariant<Dune::ALUGrid<dimGrid, dimWorld, Dune::cube, Dune::conforming>>
{
  static const int id = 1;
};


template <int dimGrid, int dimWorld>
struct ElementVariant<Dune::ALUGrid<dimGrid, dimWorld, Dune::cube, Dune::nonconforming>>
{
  static const int id = 1;
};


#endif // HAVE_ALUGRID

} // namespace internal

#if HAVE_DUNE_GRID


/**
 *  \brief  Creates a grid of a cube in various dimensions.
 *
 *          Default implementation using the Dune::StructuredGridFactory to create a grid of a cube in 1, 2 or 3
 *          dimensions. Tested with
 *          <ul><li> \c YASPGRID, \c variant 1, dim = 1, 2, 3,
 *          <li> \c SGRID, \c variant 1, dim = 1, 2, 3,
 *          <li> \c ALUGRID_SIMPLEX, \c variant 2, dim = 2, 3,
 *          <li> \c ALUGRID_CONFORM, \c variant 2, dim = 2, 2 and
 *          <li> \c ALUGRID_CUBE, \c variant 1, dim = 2, 3.</ul>
 *  \tparam GridImp
 *          Type of the underlying grid.
 *  \tparam variant
 *          Type of the codim 0 elements:
 *          <ul><li>\c 1: cubes
 *          <li>2: simplices</ul>
 **/
template <typename GridImp, int variant = internal::ElementVariant<GridImp>::id>
class Cube : public ProviderInterface<GridImp>
{
  typedef ProviderInterface<GridImp> BaseType;
  typedef Cube<GridImp, variant> ThisType;

public:
  using typename BaseType::GridType;
  static const unsigned int dimDomain = BaseType::dimDomain;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;

  static const std::string static_id()
  {
    return BaseType::static_id() + ".cube";
  }

  static Common::ConfigTree default_config(const std::string sub_name = "")
  {
    Common::ConfigTree config;
    config["lower_left"]   = "[0.0 0.0 0.0]";
    config["upper_right"]  = "[1.0 1.0 1.0]";
    config["num_elements"] = "[8 8 8]";
    if (sub_name.empty())
      return config;
    else {
      Common::ConfigTree tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr<ThisType> create(const Common::ConfigTree config = default_config(),
                                          const std::string sub_name = static_id())
  {
    // get correct config
    const Common::ConfigTree cfg         = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::ConfigTree default_cfg = default_config();
    return Common::make_unique<ThisType>(
        cfg.get("lower_left", default_cfg.get<DomainType>("lower_left")),
        cfg.get("upper_right", default_cfg.get<DomainType>("upper_right")),
        cfg.get("num_elements", default_cfg.get<std::vector<unsigned int>>("num_elements"), dimDomain));
  } // ... create(...)

  /**
   *  \brief      Creates a cube.
   *  \param[in]  lower_left
   *              Used as a lower left corner (in each dimension, if scalar).
   *  \param[in]  upper_right
   *              Used as an upper right corner (in each dimension, if scalar).
   *  \param[in]  num_elements (optional)
   *              Number of elements.
   **/
  Cube(const DomainFieldType lower_left = default_config().get<DomainFieldType>("lower_left"),
       const DomainFieldType upper_right = default_config().get<DomainFieldType>("upper_right"),
       const unsigned int num_elements = default_config().get<std::vector<unsigned int>>("num_elements")[0])
    : grid_ptr_(create_grid(DomainType(lower_left), DomainType(upper_right), parse_array(num_elements)))
  {
  }

  Cube(const std::vector<DomainFieldType>& lower_left, const std::vector<DomainFieldType>& upper_right,
       const std::vector<unsigned int> num_elements = default_config().get<std::vector<unsigned int>>("num_elements"))
    : grid_ptr_(create_grid(parse_vector(lower_left), parse_vector(upper_right), parse_array(num_elements)))
  {
  }

  Cube(const DomainType& lower_left, const DomainType& upper_right,
       const unsigned int num_elements = default_config().get<std::vector<unsigned int>>("num_elements")[0])
    : grid_ptr_(create_grid(lower_left, upper_right, parse_array(num_elements)))
  {
  }

  Cube(const DomainType& lower_left, const DomainType& upper_right,
       const std::vector<unsigned int> num_elements = default_config().get<std::vector<unsigned int>>("num_elements"))
    : grid_ptr_(create_grid(lower_left, upper_right, parse_array(num_elements)))
  {
  }

  virtual GridType& grid() DS_OVERRIDE
  {
    return *grid_ptr_;
  }

  virtual const GridType& grid() const DS_OVERRIDE
  {
    return *grid_ptr_;
  }

private:
  static std::array<unsigned int, dimDomain> parse_array(const unsigned int in)
  {
    std::array<unsigned int, dimDomain> ret;
    std::fill(ret.begin(), ret.end(), in);
    return ret;
  } // ... parse_array(...)

  static std::array<unsigned int, dimDomain> parse_array(const std::vector<unsigned int>& in)
  {
    if (in.size() < dimDomain)
      DUNE_THROW_COLORFULLY(Exceptions::wrong_input_given,
                            "Given vector is too short: should be " << dimDomain << ", is " << in.size() << "!");
    std::array<unsigned int, dimDomain> ret;
    for (unsigned int ii = 0; ii < dimDomain; ++ii)
      ret[ii] = in[ii];
    return ret;
  } // ... parse_array(...)

  static DomainType parse_vector(const std::vector<DomainFieldType>& in)
  {
    if (in.size() < dimDomain)
      DUNE_THROW_COLORFULLY(Exceptions::wrong_input_given,
                            "Given vector is too short: should be " << dimDomain << ", is " << in.size() << "!");
    DomainType ret;
    for (unsigned int ii = 0; ii < dimDomain; ++ii)
      ret[ii] = in[ii];
    return ret;
  } // ... parse_vector(...)

  static std::shared_ptr<GridType> create_grid(DomainType lower_left, DomainType upper_right,
                                               const std::array<unsigned int, dimDomain>& num_elements)
  {
    static_assert(variant == 1 || variant == 2, "variant has to be 1 or 2!");
    for (unsigned int dd = 0; dd < dimDomain; ++dd) {
      if (!(lower_left[dd] < upper_right[dd]))
        DUNE_THROW_COLORFULLY(Exceptions::wrong_input_given,
                              "lower_left has to be elementwise smaller than upper_right!\n\n" << lower_left[dd]
                                                                                               << " vs. "
                                                                                               << upper_right[dd]);
    }
    switch (variant) {
      case 1:
        return Dune::StructuredGridFactory<GridType>::createCubeGrid(lower_left, upper_right, num_elements);
        break;
      case 2:
      default:
        return Dune::StructuredGridFactory<GridType>::createSimplexGrid(lower_left, upper_right, num_elements);
        break;
    }
  } // ... create_grid(...)

  std::shared_ptr<GridType> grid_ptr_;
}; // class Cube


#else // HAVE_DUNE_GRID


template <typename GridImp, int variant = 1>
class Cube
{
  static_assert(AlwaysFalse<GridImp>::value, "You are missing dune-grid!");
};


#endif // HAVE_DUNE_GRID

} // namespace Providers
} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_GRIDS_PROVIDER_CUBE_HH
