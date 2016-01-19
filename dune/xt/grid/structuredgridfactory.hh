// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014 - 2015)
//   Rene Milk       (2014 - 2015)
//   Sven Kaulmann   (2014)
//   Tobias Leibner  (2014)

#ifndef DUNE_STUFF_GRID_STRUCTURED_GRID_FACTORY_HH
#define DUNE_STUFF_GRID_STRUCTURED_GRID_FACTORY_HH

#include <memory>

// nothing here will compile w/o grid present
#if HAVE_DUNE_GRID

#include <dune/common/unused.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#if HAVE_DUNE_SPGRID
#include <dune/grid/spgrid.hh>
#endif

#include <dune/stuff/aliases.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/common/float_cmp.hh>

namespace Dune {

#if HAVE_DUNE_SPGRID

/** \brief Specialization of the StructuredGridFactory for SPGrid
 *
 *  This allows a SPGrid to be constructed using the
 *  StructuredGridFactory just like the unstructured Grids. Limitations:
 *  \li SPGrid does not support simplices
 */
template <class ct, int dim, template <int> class Refinement, class Comm>
class StructuredGridFactory<SPGrid<ct, dim, Refinement, Comm>>
{
  typedef SPGrid<ct, dim, Refinement, Comm> GridType;
  typedef typename GridType::ctype ctype;
  static const int dimworld = GridType::dimensionworld;

public:
  /** \brief Create a structured cube grid
   *
   *  \param lowerLeft  Lower left corner of the grid
   *  \param upperRight Upper right corner of the grid
   *  \param elements   Number of elements in each coordinate direction
   */
  static std::shared_ptr<GridType> createCubeGrid(const FieldVector<ctype, dimworld>& lowerLeft,
                                                  const FieldVector<ctype, dimworld>& upperRight,
                                                  const array<unsigned int, dim>& elements)
  {
    Dune::array<int, dim> cells;
    for (const auto i : DSC::valueRange(dim))
      cells[i] = elements[i];
    return std::make_shared<GridType>(lowerLeft, upperRight, cells);
  }

  /** \brief Create a structured cube grid
   *
   *  \param lowerLeft  Lower left corner of the grid
   *  \param upperRight Upper right corner of the grid
   *  \param elements   Number of elements in each coordinate direction
   *  \param overlap    Size of overlap in each coordinate direction
   */
  static std::shared_ptr<GridType>
  createCubeGrid(const FieldVector<ctype, dimworld>& lowerLeft, const FieldVector<ctype, dimworld>& upperRight,
                 const array<unsigned int, dim>& elements, const array<unsigned int, dim>& overlap,
                 Dune::MPIHelper::MPICommunicator communicator = Dune::MPIHelper::getCommunicator())
  {
    Dune::array<int, dim> cells;
    Dune::array<int, dim> over;
    for (const auto i : DSC::valueRange(dim)) {
      cells[i] = elements[i];
      over[i]  = overlap[i];
    }
    return std::make_shared<GridType>(lowerLeft, upperRight, cells, over, communicator);
  }

  /** \brief Create a structured simplex grid
   *
   *  \param lowerLeft  Lower left corner of the grid
   *  \param upperRight Upper right corner of the grid
   *  \param elements   Number of elements in each coordinate direction
   *
   *  \note Simplices are not supported in SGrid, so this functions
   *        unconditionally throws a GridError.
   */
  static shared_ptr<GridType> createSimplexGrid(const FieldVector<ctype, dimworld>& /*lowerLeft*/,
                                                const FieldVector<ctype, dimworld>& /*upperRight*/,
                                                const array<unsigned int, dim>& /*elements*/)
  {
    DUNE_THROW(GridError,
               className<StructuredGridFactory>() << "::createSimplexGrid(): Simplices are not supported "
                                                     "by SPGrid.");
  }
};

#endif // HAVE_DUNE_SPGRID

namespace Stuff {
namespace Grid {

template <class GridType>
static Dune::array<unsigned int, GridType::dimension> default_overlap()
{
  Dune::array<unsigned int, GridType::dimension> olvp;
  olvp.fill(1u);
  return olvp;
}

//! general case just forwards to Dune core implementation
template <class GridType>
class StructuredGridFactory : public Dune::StructuredGridFactory<GridType>
{
  typedef typename GridType::ctype ctype;

public:
  static std::shared_ptr<GridType>
  createCubeGrid(const Dune::FieldVector<ctype, GridType::dimension>& lowerLeft,
                 const Dune::FieldVector<ctype, GridType::dimension>& upperRight,
                 const Dune::array<unsigned int, GridType::dimension>& elements,
                 Dune::array<unsigned int, GridType::dimension> /*overlap*/ = default_overlap<GridType>(),
                 Dune::MPIHelper::MPICommunicator /*communicator*/          = Dune::MPIHelper::getCommunicator())
  {
    return Dune::StructuredGridFactory<GridType>::createCubeGrid(lowerLeft, upperRight, elements);
  }
};

#if HAVE_DUNE_SPGRID
template <class ct, int dim, template <int> class Refinement, class Comm>
class StructuredGridFactory<Dune::SPGrid<ct, dim, Refinement, Comm>>
    : public Dune::StructuredGridFactory<Dune::SPGrid<ct, dim, Refinement, Comm>>
{
  typedef Dune::SPGrid<ct, dim, Refinement, Comm> GridType;
  typedef typename GridType::ctype ctype;

public:
  static std::shared_ptr<GridType>
  createCubeGrid(const Dune::FieldVector<ctype, dim>& lowerLeft, const Dune::FieldVector<ctype, dim>& upperRight,
                 const Dune::array<unsigned int, dim>& elements,
                 Dune::array<unsigned int, dim> overlap        = default_overlap<GridType>(),
                 Dune::MPIHelper::MPICommunicator communicator = Dune::MPIHelper::getCommunicator())
  {
    return Dune::StructuredGridFactory<GridType>::createCubeGrid(
        lowerLeft, upperRight, elements, overlap, communicator);
  }
};

#endif

template <int dim, class Coords>
class StructuredGridFactory<Dune::YaspGrid<dim, Coords>>
    : public Dune::StructuredGridFactory<Dune::YaspGrid<dim, Coords>>
{
  typedef Dune::YaspGrid<dim, Coords> GridType;
  typedef typename GridType::ctype ctype;

public:
  //! creates a non-periodic Yasp grid
  static std::shared_ptr<GridType>
  createCubeGrid(const Dune::FieldVector<ctype, dim>& lowerLeft, const Dune::FieldVector<ctype, dim>& upperRight,
                 const Dune::array<unsigned int, dim>& elements_in,
                 Dune::array<unsigned int, dim> overlap        = default_overlap<GridType>(),
                 Dune::MPIHelper::MPICommunicator communicator = Dune::MPIHelper::getCommunicator())
  {
    const auto no_periodic_direction = std::bitset<dim>();
    std::array<int, dim> elements;
    std::copy(elements_in.begin(), elements_in.end(), elements.begin());
    auto overlap_check = overlap;
    overlap_check.fill(overlap[0]);
    for (auto i : DSC::valueRange(1, dim))
      if (overlap[i] != overlap[0])
        DUNE_THROW(Dune::InvalidStateException, "YaspGrid only supports uniform overlap");
    return std::make_shared<GridType>(lowerLeft, upperRight, elements, no_periodic_direction, overlap[0], communicator);
  }
};

} // namespace Grid
} // end of namespace Stuff
} // namespace Dune

#endif // HAVE_DUNE_GRID

#endif // DUNE_STUFF_GRID_STRUCTURED_GRID_FACTORY_HH
