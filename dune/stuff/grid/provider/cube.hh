#ifndef DUNE_STUFF_GRID_PROVIDER_CUBE_HH
#define DUNE_STUFF_GRID_PROVIDER_CUBE_HH

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#elif defined(HAVE_CONFIG_H)
#include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

#if HAVE_DUNE_GRID

#include <sstream>
#include <type_traits>
#include <boost/assign/list_of.hpp>

#include <dune/common/parametertree.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/static_assert.hh>

#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>
#ifdef HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif
#include <dune/grid/sgrid.hh>

#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/color.hh>

#include "interface.hh"

namespace Dune {
namespace Stuff {
namespace Grid {
namespace Provider {

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
template <typename GridImp, int variant>
class GenericCube : public Interface<GridImp>
{
public:
  //! Type of the provided grid.
  typedef GridImp GridType;

  typedef Interface<GridType> BaseType;

  typedef GenericCube<GridType, variant> ThisType;

private:
  typedef typename GridType::LeafGridView GridViewType;

public:
  static const unsigned int dim = BaseType::dim;

  //! Type of the grids coordinates.
  typedef typename GridType::ctype ctype;
  typedef Dune::FieldVector<ctype, dim> CoordinateType;

  //! Unique identifier: \c stuff.grid.provider.cube
  static const std::string id()
  {
    return BaseType::id() + ".cube";
  }

  /**
   *  \brief      Creates a cube.
   *  \param[in]  lowerLeft
   *              A double that is used as a lower left corner in each dimension.
   *  \param[in]  upperRight
   *              A double that is used as a upper right corner in each dimension.
   *  \param[in]  numElements (optional)
   *              number of elements.
   **/
  GenericCube(const double _lowerLeft = 0.0, const double _upperRight = 1.0, const unsigned int numElements = 1u)
    : lowerLeft_(_lowerLeft)
    , upperRight_(_upperRight)
  {
    Dune::array<unsigned int, dim> tmpNumElements;
    std::fill(tmpNumElements.begin(), tmpNumElements.end(), numElements);
    buildGrid(tmpNumElements);
  }

  /**
   *  \brief      Creates a cube.
   *  \param[in]  lowerLeft
   *              A vector that is used as a lower left corner.
   *  \param[in]  upperRight
   *              A vector that is used as a upper right corner.
   *  \param[in]  numElements (optional)
   *              number of elements.
   **/
  GenericCube(const CoordinateType& _lowerLeft, const CoordinateType& _upperRight, const unsigned int numElements = 1u)
    : lowerLeft_(_lowerLeft)
    , upperRight_(_upperRight)
  {
    Dune::array<unsigned int, dim> tmpNumElements;
    std::fill(tmpNumElements.begin(), tmpNumElements.end(), numElements);
    buildGrid(tmpNumElements);
  }

  /**
    \brief      Creates a cube. This signature allows to prescribe anisotopic refinement
    \param[in]  lowerLeft
                A double that is used as a lower left corner in each dimension.
    \param[in]  upperRight
                A double that is used as a upper right corner in each dimension.
    \param[in]  numElements
                number of elements in each dimension.
                can contain 0 to dim elements (missing dimension are initialized to 1)
    \tparam Coord anything that CoordinateType allows to copy construct from
    \tparam ContainerType some sequence type that functions with std::begin/end
    \tparam T an unsigned integral Type
    **/
  template <class ContainerType>
  GenericCube(const CoordinateType& _lowerLeft, const CoordinateType& _upperRight,
              const ContainerType numElements = boost::assign::list_of<typename ContainerType::value_type>().repeat(
                  dim, typename ContainerType::value_type(1u)))
    : lowerLeft_(_lowerLeft)
    , upperRight_(_upperRight)
  {
    Dune::array<unsigned int, dim> tmpNumElements;
    static_assert(std::is_unsigned<typename ContainerType::value_type>::value
                      && std::is_integral<typename ContainerType::value_type>::value,
                  "only unsigned integral number of elements per dimension allowed");
    // base init in case input is shorter
    std::fill(tmpNumElements.begin(), tmpNumElements.end(), 1u);
    std::copy(numElements.begin(), numElements.end(), tmpNumElements.begin());
    buildGrid(tmpNumElements);
  }

  GenericCube(ThisType& other)
    : lowerLeft_(other.lowerLeft_)
    , upperRight_(other.upperRight_)
    , grid_(other.grid_)
  {
  }

  GenericCube(const ThisType& other)
    : lowerLeft_(other.lowerLeft_)
    , upperRight_(other.upperRight_)
    , grid_(other.grid_)
  {
  }

  /**
   *  \brief      Creates a cube.
   *  \param[in]  paramTree
   *              A Dune::ParameterTree containing
   *              <ul><li> the following keys directly or
   *              <li> a subtree named Cube::id, containing the following keys. If a subtree is present, it is always
   *selected. Also it is solely selceted, so that all keys in the supertree are ignored.</ul>
   *              The actual keys are:
   *              <ul><li> \c lowerLeft: \a double or a vector that is used as lower left corners.
   *              <li> \c upperRight: \a double or a vector that is used as upper right corners.
   *              <li> \c numElements: \a int or vector to denote the number of elements.
   *              </ul>
   **/
  static ThisType createFromParamTree(const Dune::ParameterTree& paramTree, const std::string subName = id())
  {
    // get correct paramTree
    Dune::Stuff::Common::ParameterTreeX extendedParamTree;
    if (paramTree.hasSub(subName))
      extendedParamTree = paramTree.sub(subName);
    else
      extendedParamTree = paramTree;
    // get lower left
    std::vector<ctype> lowerLefts;
    if (extendedParamTree.hasVector("lowerLeft")) {
      lowerLefts = extendedParamTree.getVector("lowerLeft", ctype(0), dim);
      assert(lowerLefts.size() >= dim && "Given vector too short!");
    } else if (extendedParamTree.hasKey("lowerLeft")) {
      const ctype lowerLeft = extendedParamTree.get("lowerLeft", ctype(0));
      lowerLefts            = std::vector<ctype>(dim, lowerLeft);
    } else {
      std::cout << "\n" << Dune::Stuff::Common::colorString("WARNING in " + id() + ":")
                << " neither vector nor key 'lowerLeft' given, defaulting to 0.0!" << std::flush;
      lowerLefts = std::vector<ctype>(dim, ctype(0));
    }
    // get upper right
    std::vector<ctype> upperRights;
    if (extendedParamTree.hasVector("upperRight")) {
      upperRights = extendedParamTree.getVector("upperRight", ctype(1), dim);
      assert(upperRights.size() >= dim && "Given vector too short!");
    } else if (extendedParamTree.hasKey("upperRight")) {
      const ctype upperRight = extendedParamTree.get("upperRight", ctype(1));
      upperRights            = std::vector<ctype>(dim, upperRight);
    } else {
      std::cout << "\n" << Dune::Stuff::Common::colorString("WARNING in " + id() + ":")
                << " neither vector nor key 'upperRight' given, defaulting to 1.0!" << std::flush;
      upperRights = std::vector<ctype>(dim, ctype(1));
    }
    // get number of elements
    std::vector<unsigned int> tmpNumElements;
    if (extendedParamTree.hasVector("numElements")) {
      tmpNumElements = extendedParamTree.getVector("numElements", 1u, dim);
      assert(tmpNumElements.size() >= dim && "Given vector too short!");
    } else if (extendedParamTree.hasKey("numElements")) {
      const unsigned int numElement = extendedParamTree.get("numElements", 1u);
      tmpNumElements                = std::vector<unsigned int>(dim, numElement);
    } else {
      std::cout << "\n" << Dune::Stuff::Common::colorString("WARNING in " + id() + ":")
                << " neither vector nor key 'numElements' given, defaulting to 1!" << std::flush;
      tmpNumElements = std::vector<unsigned int>(dim, 1u);
    }
    // check and save
    CoordinateType lowerLeft;
    CoordinateType upperRight;
    Dune::array<unsigned int, dim> numElements;
    for (unsigned int d = 0; d < dim; ++d) {
      assert(lowerLefts[d] < upperRights[d]
             && "Given 'upperRight' hast to be elementwise larger than given 'lowerLeft'!");
      lowerLeft[d]  = lowerLefts[d];
      upperRight[d] = upperRights[d];
      assert(tmpNumElements[d] > 0 && "Given 'numElements' has to be elementwise positive!");
      numElements[d] = tmpNumElements[d];
    }
    return GenericCube(lowerLeft, upperRight, numElements);
  } // static ThisType createFromParamTree(const Dune::ParameterTree& paramTree, const std::string subName = id())

  ThisType& operator=(ThisType& other)
  {
    if (this != &other) {
      lowerLeft_  = other.lowerLeft();
      upperRight_ = other.upperRight();
      grid_       = other.grid();
    }
    return this;
  } // ThisType& operator=(ThisType& other)

  //! access to shared ptr
  virtual Dune::shared_ptr<GridType> grid()
  {
    return grid_;
  }

  virtual const Dune::shared_ptr<const GridType> grid() const
  {
    return grid_;
  }

  const CoordinateType& lowerLeft() const
  {
    return lowerLeft_;
  }

  const CoordinateType& upperRight() const
  {
    return upperRight_;
  }

  /**
   *  \brief      Visualizes the grid using Dune::VTKWriter.
   *  \param[in]  filename
   **/

private:
  void buildGrid(const Dune::array<unsigned int, dim>& numElements)
  {
    dune_static_assert(variant >= 1 && variant <= 2, "only variant 1 and 2 are valid");
    switch (variant) {
      case 1:
        grid_ = Dune::StructuredGridFactory<GridType>::createCubeGrid(lowerLeft_, upperRight_, numElements);
        break;
      case 2:
      default:
        grid_ = Dune::StructuredGridFactory<GridType>::createSimplexGrid(lowerLeft_, upperRight_, numElements);
        break;
    }
  } // void buildGrid(const CoordinateType& lowerLeft, const CoordinateType& upperRight)

  CoordinateType lowerLeft_;
  CoordinateType upperRight_;
  Dune::shared_ptr<GridType> grid_;
}; // class GenericCube

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

// default implementation of a cube for any grid
// tested for
// dim = 2
//  ALUGRID_SIMPLEX, variant 2
//  ALUGRID_CONFORM, variant 2
// dim = 3
//  ALUGRID_SIMPLEX, variant 2
#if defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template <class GridType = Dune::GridSelector::GridType>
#else // defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template <class GridType = Dune::SGrid<2, 2>>
#endif // defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
class Cube : public GenericCube<GridType, ElementVariant<GridType>::id>
{
public:
  typedef GenericCube<GridType, ElementVariant<GridType>::id> BaseType;

  typedef Cube<GridType> ThisType;

  typedef typename BaseType::CoordinateType CoordinateType;

  Cube(const double _lowerLeft = 0.0, const double _upperRight = 1.0, const unsigned int _numElements = 1)
    : BaseType(_lowerLeft, _upperRight, _numElements)
  {
  }

  Cube(const CoordinateType& _lowerLeft, const CoordinateType& _upperRight, const unsigned int _numElements = 1)
    : BaseType(_lowerLeft, _upperRight, _numElements)
  {
  }

  template <class ContainerType>
  Cube(const CoordinateType& lowerLeft, const CoordinateType& upperRight,
       const ContainerType numElements = boost::assign::list_of<typename ContainerType::value_type>().repeat(
           GridType::dimensionworld, typename ContainerType::value_type(1u)))
    : BaseType(lowerLeft, upperRight, numElements)
  {
  }

  Cube(BaseType& other)
    : BaseType(other)
  {
  }

  static ThisType createFromParamTree(const Dune::ParameterTree& paramTree, const std::string subName = BaseType::id())
  {
    BaseType base = BaseType::createFromParamTree(paramTree, subName);
    return ThisType(base);
  }
}; // class Cube

} // namespace Provider
} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // HAVE_DUNE_GRID

#endif // DUNE_STUFF_GRID_PROVIDER_CUBE_HH
