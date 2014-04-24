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

#include <dune/common/static_assert.hh>

#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/yaspgrid.hh>
#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/configtree.hh>
#include <dune/stuff/common/memory.hh>

#include "interface.hh"

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


template <int dim>
struct ElementVariant<Dune::SGrid<dim, dim>>
{
  static const int id = 1;
};


#if HAVE_ALUGRID


template <int dim>
struct ElementVariant<Dune::ALUCubeGrid<dim, dim>>
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

  virtual std::string id() /*DS_OVERRIDE*/
  {
    return BaseType::static_id() + ".cube";
  }

  /**
   *  \brief      Creates a cube.
   *  \param[in]  _lowerLeft
   *              A double that is used as a lower left corner in each dimension.
   *  \param[in]  _upperRight
   *              A double that is used as a upper right corner in each dimension.
   *  \param[in]  numElements (optional)
   *              number of elements.
   **/
  Cube(const double _lowerLeft = 0.0, const double _upperRight = 1.0, const unsigned int numElements = 1u)
    : lower_left_(_lowerLeft)
    , upper_right_(_upperRight)
  {
    std::array<unsigned int, dimDomain> tmpNumElements;
    std::fill(tmpNumElements.begin(), tmpNumElements.end(), numElements);
    buildGrid(tmpNumElements);
  }

  Cube(const std::vector<double>& lL, const std::vector<double>& uR, const std::vector<unsigned int>& nE)
  {
    if (lL.size() < dimDomain)
      DUNE_THROW_COLORFULLY(Exceptions::wrong_input_given,
                            "lL has to be at least of size " << dimDomain << " (is " << lL.size() << ")!");
    if (uR.size() < dimDomain)
      DUNE_THROW_COLORFULLY(Exceptions::wrong_input_given,
                            "uR has to be at least of size " << dimDomain << " (is " << uR.size() << ")!");
    if (nE.size() < dimDomain)
      DUNE_THROW_COLORFULLY(Exceptions::wrong_input_given,
                            "nE has to be at least of size " << dimDomain << " (is " << nE.size() << ")!");
    std::array<unsigned int, dimDomain> num_elements;
    for (unsigned int ii = 0; ii < dimDomain; ++ii) {
      lower_left_[ii]  = lL[ii];
      upper_right_[ii] = uR[ii];
      num_elements[ii] = nE[ii];
    }
    buildGrid(num_elements);
  }

  /**
   *  \brief      Creates a cube.
   *  \param[in]  _lowerLeft
   *              A vector that is used as a lower left corner.
   *  \param[in]  _upperRight
   *              A vector that is used as a upper right corner.
   *  \param[in]  numElements (optional)
   *              number of elements.
   **/
  Cube(const DomainType& _lowerLeft, const DomainType& _upperRight, const unsigned int numElements = 1u)
    : lower_left_(_lowerLeft)
    , upper_right_(_upperRight)
  {
    std::array<unsigned int, dimDomain> tmpNumElements;
    std::fill(tmpNumElements.begin(), tmpNumElements.end(), numElements);
    buildGrid(tmpNumElements);
  }

  /**
    \brief      Creates a cube. This signature allows to prescribe anisotopic refinement
    \param[in]  _lowerLeft
                A double that is used as a lower left corner in each dimension.
    \param[in]  _upperRight
                A double that is used as a upper right corner in each dimension.
    \param[in]  numElements
                number of elements in each dimension.
                can contain 0 to dim elements (missing dimension are initialized to 1)
    \tparam ContainerType some sequence type that functions with std::begin/end
    **/
  template <class ContainerType>
  Cube(const DomainType& _lowerLeft, const DomainType& _upperRight,
       const ContainerType numElements = boost::assign::list_of<typename ContainerType::value_type>().repeat(
           dimDomain, typename ContainerType::value_type(1u)))
    : lower_left_(_lowerLeft)
    , upper_right_(_upperRight)
  {
    std::array<unsigned int, dimDomain> tmpNumElements;
    static_assert(std::is_unsigned<typename ContainerType::value_type>::value
                      && std::is_integral<typename ContainerType::value_type>::value,
                  "only unsigned integral number of elements per dimension allowed");
    // base init in case input is shorter
    std::fill(tmpNumElements.begin(), tmpNumElements.end(), 1u);
    std::copy(numElements.begin(), numElements.end(), tmpNumElements.begin());
    buildGrid(tmpNumElements);
  }

  static ParameterTree default_config(const std::string sub_name = "")
  {
    Common::ConfigTree config;
    config["lower_left"]   = "[0.0; 0.0; 0.0]";
    config["upper_right"]  = "[1.0; 1.0; 1.0]";
    config["num_elements"] = "[8; 8; 8]";
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
    Common::ConfigTree cfg;
    if (config.has_sub(sub_name))
      cfg = config.sub(sub_name);
    else
      cfg = config;
    // extract needed data
    try {
      // try to get vectors
      auto lower_left   = cfg.get<DomainType>("lower_left", dimDomain);
      auto upper_right  = cfg.get<DomainType>("upper_right", dimDomain);
      auto num_elements = cfg.get<std::vector<unsigned int>>("num_elements", dimDomain);
      return Common::make_unique<ThisType>(lower_left, upper_right, num_elements);
    } catch (Exceptions::configuration_error&) {
      // get scalars
      auto lower_left   = cfg.get<DomainFieldType>("lower_left");
      auto upper_right  = cfg.get<DomainFieldType>("upper_right");
      auto num_elements = cfg.get<unsigned int>("num_elements");
      return Common::make_unique<ThisType>(lower_left, upper_right, num_elements);
    }
  } // ... create(...)

  virtual std::shared_ptr<GridType>& grid() /*DS_OVERRIDE*/
  {
    return grid_;
  }

  virtual const std::shared_ptr<const GridType>& grid() const /*DS_OVERRIDE*/
  {
    return grid_;
  }

  const DomainType& lower_left() const
  {
    return lower_left_;
  }

  const DomainType& upper_right() const
  {
    return upper_right_;
  }

private:
  void buildGrid(const std::array<unsigned int, dimDomain>& numElements)
  {
    static_assert(variant == 1 || variant == 2, "variant has to be 1 or 2!");
    for (unsigned int dd = 0; dd < dimDomain; ++dd) {
      if (!(lower_left_[dd] < upper_right_[dd]))
        DUNE_THROW_COLORFULLY(Exceptions::wrong_input_given,
                              "lower_left has to be elementwise smaller than upper_right!\n\n" << lower_left_[dd]
                                                                                               << " vs. "
                                                                                               << upper_right_[dd]);
    }
    switch (variant) {
      case 1:
        grid_ = Dune::StructuredGridFactory<GridType>::createCubeGrid(lower_left_, upper_right_, numElements);
        break;
      case 2:
      default:
        grid_ = Dune::StructuredGridFactory<GridType>::createSimplexGrid(lower_left_, upper_right_, numElements);
        break;
    }
  } // void buildGrid(const DomainType& lowerLeft, const DomainType& upperRight)

  DomainType lower_left_;
  DomainType upper_right_;
  std::shared_ptr<GridType> grid_;
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
