#ifndef DUNE_STUFF_GRID_PROVIDER_CUBE_HH
#define DUNE_STUFF_GRID_PROVIDER_CUBE_HH

#if HAVE_DUNE_GRID

#include <memory>
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

template <typename GridType>
struct ElementVariant;

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
template <typename GridImp, int variant = ElementVariant<GridImp>::id>
class GridProviderCube : public GridProviderInterface<GridImp>
{
public:
  //! Type of the provided grid.
  typedef GridImp GridType;

  typedef GridProviderInterface<GridType> BaseType;

  typedef GridProviderCube<GridType, variant> ThisType;

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
   *  \param[in]  _lowerLeft
   *              A double that is used as a lower left corner in each dimension.
   *  \param[in]  _upperRight
   *              A double that is used as a upper right corner in each dimension.
   *  \param[in]  numElements (optional)
   *              number of elements.
   **/
  GridProviderCube(const double _lowerLeft = 0.0, const double _upperRight = 1.0, const unsigned int numElements = 1u)
    : lowerLeft_(_lowerLeft)
    , upperRight_(_upperRight)
  {
    Dune::array<unsigned int, dim> tmpNumElements;
    std::fill(tmpNumElements.begin(), tmpNumElements.end(), numElements);
    buildGrid(tmpNumElements);
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
  GridProviderCube(const CoordinateType& _lowerLeft, const CoordinateType& _upperRight,
                   const unsigned int numElements = 1u)
    : lowerLeft_(_lowerLeft)
    , upperRight_(_upperRight)
  {
    Dune::array<unsigned int, dim> tmpNumElements;
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
  GridProviderCube(const CoordinateType& _lowerLeft, const CoordinateType& _upperRight,
                   const ContainerType numElements = boost::assign::list_of<typename ContainerType::value_type>()
                                                         .repeat(dim, typename ContainerType::value_type(1u)))
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

  static Dune::ParameterTree defaultSettings(const std::string subName = "")
  {
    Dune::ParameterTree description;
    description["lowerLeft"]   = "[0.0; 0.0; 0.0]";
    description["upperRight"]  = "[1.0; 1.0; 1.0]";
    description["numElements"] = "[4; 4; 4]";
    if (subName.empty())
      return description;
    else {
      Dune::Stuff::Common::ExtendedParameterTree extendedDescription;
      extendedDescription.add(description, subName);
      return extendedDescription;
    }
  } // ... createDefaultSettings(...)

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
  static ThisType* create(const Dune::ParameterTree& _settings, const std::string subName = id())
  {
    // get correct _settings
    Dune::Stuff::Common::ExtendedParameterTree settings;
    if (_settings.hasSub(subName))
      settings = _settings.sub(subName);
    else
      settings = _settings;
    // get lower left
    std::vector<ctype> lowerLefts;
    if (settings.hasVector("lowerLeft")) {
      lowerLefts = settings.getVector("lowerLeft", ctype(0), dim);
      assert(lowerLefts.size() >= dim && "Given vector too short!");
    } else if (settings.hasKey("lowerLeft")) {
      const ctype lowerLeft = settings.get("lowerLeft", ctype(0));
      lowerLefts            = std::vector<ctype>(dim, lowerLeft);
    } else {
      std::cout << "\n" << Dune::Stuff::Common::colorString("WARNING in " + id() + ":")
                << " neither vector nor key 'lowerLeft' given, defaulting to 0.0!" << std::flush;
      lowerLefts = std::vector<ctype>(dim, ctype(0));
    }
    // get upper right
    std::vector<ctype> upperRights;
    if (settings.hasVector("upperRight")) {
      upperRights = settings.getVector("upperRight", ctype(1), dim);
      assert(upperRights.size() >= dim && "Given vector too short!");
    } else if (settings.hasKey("upperRight")) {
      const ctype upperRight = settings.get("upperRight", ctype(1));
      upperRights            = std::vector<ctype>(dim, upperRight);
    } else {
      std::cout << "\n" << Dune::Stuff::Common::colorString("WARNING in " + id() + ":")
                << " neither vector nor key 'upperRight' given, defaulting to 1.0!" << std::flush;
      upperRights = std::vector<ctype>(dim, ctype(1));
    }
    // get number of elements
    std::vector<unsigned int> tmpNumElements;
    if (settings.hasVector("numElements")) {
      tmpNumElements = settings.getVector("numElements", 1u, dim);
      assert(tmpNumElements.size() >= dim && "Given vector too short!");
    } else if (settings.hasKey("numElements")) {
      const unsigned int numElement = settings.get("numElements", 1u);
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
    return new ThisType(lowerLeft, upperRight, numElements);
  } // ... create(...)

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
}; // class GridProviderCube

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


} // namespace Stuff
} // namespace Dune

#endif // HAVE_DUNE_GRID

#endif // DUNE_STUFF_GRID_PROVIDER_CUBE_HH
