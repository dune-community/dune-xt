// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Sven Kaulmann

#ifndef DUNE_STUFF_GRID_STRUCTURED_GRID_FACTORY_HH
#define DUNE_STUFF_GRID_STRUCTURED_GRID_FACTORY_HH

#include <dune/stuff/common/disable_warnings.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/stuff/common/reenable_warnings.hh>


namespace Dune {

#if HAVE_DUNE_SPGRID

#include <dune/grid/spgrid.hh>
#include <dune/stuff/aliases.hh>

/** \brief Specialization of the StructuredGridFactory for SPGrid
 *
 *  This allows a SPGrid to be constructed using the
 *  StructuredGridFactory just like the unstructured Grids. Limitations:
 *  \li SPGrid does not support simplices
 */
template <class ct, int dim, SPRefinementStrategy strategy, class Comm>
class StructuredGridFactory<SPGrid<ct, dim, strategy, Comm>>
{
  typedef SPGrid<ct, dim, strategy, Comm> GridType;
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
  static std::shared_ptr<GridType> createCubeGrid(const FieldVector<ctype, dimworld>& lowerLeft,
                                                  const FieldVector<ctype, dimworld>& upperRight,
                                                  const array<unsigned int, dim>& elements,
                                                  const array<unsigned int, dim>& overlap)
  {
    Dune::array<int, dim> cells;
    Dune::array<int, dim> over;
    for (const auto i : DSC::valueRange(dim)) {
      cells[i] = elements[i];
      over[i]  = overlap[i];
    }
    return std::make_shared<GridType>(lowerLeft, upperRight, cells, over);
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
  static shared_ptr<GridType> createSimplexGrid(const FieldVector<ctype, dimworld>& lowerLeft,
                                                const FieldVector<ctype, dimworld>& upperRight,
                                                const array<unsigned int, dim>& elements)
  {
    DUNE_THROW(GridError,
               className<StructuredGridFactory>() << "::createSimplexGrid(): Simplices are not supported "
                                                     "by SPGrid.");
  }
};


#endif // HAVE_DUNE_SPGRID


/** \brief Specialization of the StructuredGridFactory for SGrid< dim, dimWorld >
 *
 *  This allows a SGrid to be constructed using the
 *  StructuredGridFactory just like the unstructured Grids. Limitations:
 *  \li SGrid does not support simplices
 */
template <int dim, int dimworld>
class StructuredGridFactory<SGrid<dim, dimworld>>
{
  typedef SGrid<dim, dimworld> GridType;
  typedef typename GridType::ctype ctype;

public:
  /** \brief Create a structured cube grid
   *
   *  \param lowerLeft  Lower left corner of the grid
   *  \param upperRight Upper right corner of the grid
   *  \param elements   Number of elements in each coordinate direction
   */
  static shared_ptr<GridType> createCubeGrid(const FieldVector<ctype, dim>& lowerLeft,
                                             const FieldVector<ctype, dim>& upperRight,
                                             const array<unsigned int, dim>& elements)
  {
    FieldVector<int, dim> elements_;
    std::copy(elements.begin(), elements.end(), elements_.begin());

    return shared_ptr<GridType>(new GridType(elements_, lowerLeft, upperRight));
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
  static shared_ptr<GridType> createSimplexGrid(const FieldVector<ctype, dim>& lowerLeft,
                                                const FieldVector<ctype, dim>& upperRight,
                                                const array<unsigned int, dim>& elements)
  {
    DUNE_THROW(GridError,
               className<StructuredGridFactory>() << "::createSimplexGrid(): Simplices are not supported "
                                                     "by SGrid.");
  }
};

} // namespace Dune

#endif // DUNE_STUFF_GRID_STRUCTURED_GRID_FACTORY_HH
