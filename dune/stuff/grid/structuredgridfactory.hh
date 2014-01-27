// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_STUFF_GRID_STRUCTURED_GRID_FACTORY_HH
#define DUNE_STUFF_GRID_STRUCTURED_GRID_FACTORY_HH

#include <dune/grid/utility/structuredgridfactory.hh>

#if HAVE_DUNE_SPGRID

#include <dune/grid/spgrid.hh>
#include <dune/stuff/aliases.hh>

namespace Dune {

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

} // namespace Dune

#endif // HAVE_DUNE_SPGRID

#endif // DUNE_STUFF_GRID_STRUCTURED_GRID_FACTORY_HH
