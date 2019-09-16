// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   René Fritze     (2014 - 2016, 2018)
//   Sven Kaulmann   (2014)
//   Tobias Leibner  (2014, 2016, 2018)

#ifndef DUNE_XT_GRID_STRUCTUREDGRIDFACTORY_HH
#define DUNE_XT_GRID_STRUCTUREDGRIDFACTORY_HH

#include <memory>

#include <dune/common/unused.hh>

#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>

#if HAVE_DUNE_SPGRID
#  include <dune/grid/spgrid.hh>
#endif

#if HAVE_ALBERTA
#  include <dune/grid/albertagrid.hh>
#endif

#if HAVE_DUNE_ALUGRID
#  include <dune/alugrid/common/structuredgridfactory.hh>
#endif

#if HAVE_DUNE_UGGRID
#  include <dune/grid/uggrid.hh>
#endif

#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/common/logging.hh>

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
                                                  const std::array<unsigned int, dim>& elements)
  {
    std::array<int, dim> cells;
    for (const auto i : XT::Common::value_range(dim))
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
  createCubeGrid(const FieldVector<ctype, dimworld>& lowerLeft,
                 const FieldVector<ctype, dimworld>& upperRight,
                 const std::array<unsigned int, dim>& elements,
                 const std::array<unsigned int, dim>& overlap,
                 Dune::MPIHelper::MPICommunicator communicator = Dune::MPIHelper::getCommunicator())
  {
    std::array<int, dim> cells;
    std::array<int, dim> over;
    for (const auto i : XT::Common::value_range(dim)) {
      cells[i] = elements[i];
      over[i] = overlap[i];
    }
    return std::make_shared<GridType>(lowerLeft, upperRight, cells, over, communicator);
  }

  /** \brief Create a structured simplex grid
   *
   *  \param lowerLeft  Lower left corner of the grid
   *  \param upperRight Upper right corner of the grid
   *  \param elements   Number of elements in each coordinate direction
   *
   *  \note Simplices are not supported in SPGrid, so this functions
   *        unconditionally throws a GridError.
   */
  static shared_ptr<GridType> createSimplexGrid(const FieldVector<ctype, dimworld>& /*lowerLeft*/,
                                                const FieldVector<ctype, dimworld>& /*upperRight*/,
                                                const std::array<unsigned int, dim>& /*elements*/)
  {
    DUNE_THROW(GridError,
               className<StructuredGridFactory>() << "::createSimplexGrid(): Simplices are not supported "
                                                     "by SPGrid.");
  }
};

#endif // HAVE_DUNE_SPGRID

namespace XT {
namespace Grid {


template <class GridType>
static std::array<unsigned int, GridType::dimension> default_overlap()
{
  std::array<unsigned int, GridType::dimension> olvp;
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
                 const std::array<unsigned int, GridType::dimension>& elements,
                 std::array<unsigned int, GridType::dimension> /*overlap*/ = default_overlap<GridType>(),
                 Dune::MPIHelper::MPICommunicator /*mpi_comm*/ = Dune::MPIHelper::getCommunicator())
  {
    return Dune::StructuredGridFactory<GridType>::createCubeGrid(lowerLeft, upperRight, elements);
  }

  static std::shared_ptr<GridType>
  createSimplexGrid(const Dune::FieldVector<ctype, GridType::dimension>& lowerLeft,
                    const Dune::FieldVector<ctype, GridType::dimension>& upperRight,
                    const std::array<unsigned int, GridType::dimension>& elements,
                    Dune::MPIHelper::MPICommunicator /*mpi_comm*/ = Dune::MPIHelper::getCommunicator())
  {
    return Dune::StructuredGridFactory<GridType>::createSimplexGrid(lowerLeft, upperRight, elements);
  }
};

#if HAVE_DUNE_UGGRID
template <int dim>
class StructuredGridFactory<UGGrid<dim>> : public Dune::StructuredGridFactory<UGGrid<dim>>
{
  typedef UGGrid<dim> GridType;
  typedef typename GridType::ctype ctype;

public:
  static std::shared_ptr<GridType>
  createCubeGrid(const Dune::FieldVector<ctype, GridType::dimension>& lowerLeft,
                 const Dune::FieldVector<ctype, GridType::dimension>& upperRight,
                 const std::array<unsigned int, GridType::dimension>& elements,
                 std::array<unsigned int, GridType::dimension> = default_overlap<GridType>(),
                 Dune::MPIHelper::MPICommunicator mpi_comm = Dune::MPIHelper::getCommunicator())
  {
    if (Dune::MPIHelper::isFake)
      return Dune::StructuredGridFactory<GridType>::createCubeGrid(lowerLeft, upperRight, elements);
#  if HAVE_MPI
    if (mpi_comm == MPI_COMM_WORLD)
      return Dune::StructuredGridFactory<GridType>::createCubeGrid(lowerLeft, upperRight, elements);
#  endif
    DUNE_THROW(InvalidStateException, "Alberta construction cannot handle non-world communicators");
  }

  static std::shared_ptr<GridType>
  createSimplexGrid(const Dune::FieldVector<ctype, dim>& lowerLeft,
                    const Dune::FieldVector<ctype, dim>& upperRight,
                    const std::array<unsigned int, dim>& elements,
                    Dune::MPIHelper::MPICommunicator mpi_comm = Dune::MPIHelper::getCommunicator())
  {
    if (Dune::MPIHelper::isFake)
      return Dune::StructuredGridFactory<GridType>::createSimplexGrid(lowerLeft, upperRight, elements);
#  if HAVE_MPI
    if (mpi_comm == MPI_COMM_WORLD)
      return Dune::StructuredGridFactory<GridType>::createSimplexGrid(lowerLeft, upperRight, elements);
#  endif
    DUNE_THROW(InvalidStateException, "UGGRid construction cannot handle non-world communicators");
  }
};
#endif // HAVE_DUNE_UGGRID

#if HAVE_ALBERTA
template <int dim>
class StructuredGridFactory<AlbertaGrid<dim, dim>> : public Dune::StructuredGridFactory<AlbertaGrid<dim, dim>>
{
  typedef AlbertaGrid<dim, dim> GridType;
  typedef typename GridType::ctype ctype;

public:
  static std::shared_ptr<GridType>
  createCubeGrid(const Dune::FieldVector<ctype, GridType::dimension>& lowerLeft,
                 const Dune::FieldVector<ctype, GridType::dimension>& upperRight,
                 const std::array<unsigned int, GridType::dimension>& elements,
                 std::array<unsigned int, GridType::dimension> overlap = default_overlap<GridType>(),
                 Dune::MPIHelper::MPICommunicator mpi_comm = Dune::MPIHelper::getCommunicator())
  {
    bool warn = false;
    unsigned int sum = 0;
    for (auto i : Common::value_range(1, dim)) {
      warn = warn || (overlap[i] > 1);
      sum += overlap[i];
    }
    if (warn || sum >= dim)
      DXTC_LOG_INFO << "Ignoring non-default overlap for alberta cube creation";
    if (Dune::MPIHelper::isFake)
      return Dune::StructuredGridFactory<GridType>::createCubeGrid(lowerLeft, upperRight, elements);
#  if HAVE_MPI
    if (mpi_comm == MPI_COMM_WORLD)
      return Dune::StructuredGridFactory<GridType>::createCubeGrid(lowerLeft, upperRight, elements);
#  endif
    DUNE_THROW(InvalidStateException, "Alberta construction cannot handle non-world communicators");
  }

  static std::shared_ptr<GridType>
  createSimplexGrid(const Dune::FieldVector<ctype, dim>& lowerLeft,
                    const Dune::FieldVector<ctype, dim>& upperRight,
                    const std::array<unsigned int, dim>& elements,
                    Dune::MPIHelper::MPICommunicator mpi_comm = Dune::MPIHelper::getCommunicator())
  {
    if (Dune::MPIHelper::isFake)
      return Dune::StructuredGridFactory<GridType>::createSimplexGrid(lowerLeft, upperRight, elements);
#  if HAVE_MPI
    if (mpi_comm == MPI_COMM_WORLD)
      return Dune::StructuredGridFactory<GridType>::createSimplexGrid(lowerLeft, upperRight, elements);
#  endif
    DUNE_THROW(InvalidStateException, "Alberta construction cannot handle non-world communicators");
  }
};
#endif // HAVE_ALBERTA


#if HAVE_DUNE_ALUGRID
template <int dim_world, int dim, ALUGridRefinementType refineType, class Comm>
class StructuredGridFactory<ALUGrid<dim, dim_world, Dune::cube, refineType, Comm>>
  : public Dune::StructuredGridFactory<ALUGrid<dim, dim_world, Dune::cube, refineType, Comm>>
{
  typedef ALUGrid<dim, dim_world, Dune::cube, refineType, Comm> GridType;
  typedef typename GridType::ctype ctype;

public:
  static std::shared_ptr<GridType>
  createCubeGrid(const Dune::FieldVector<ctype, GridType::dimension>& lowerLeft,
                 const Dune::FieldVector<ctype, GridType::dimension>& upperRight,
                 const std::array<unsigned int, GridType::dimension>& elements,
                 std::array<unsigned int, GridType::dimension> = default_overlap<GridType>(),
                 Dune::MPIHelper::MPICommunicator mpi_comm = Dune::MPIHelper::getCommunicator())
  {
    if (Dune::MPIHelper::isFake)
      return Dune::StructuredGridFactory<GridType>::createCubeGrid(lowerLeft, upperRight, elements);
#  if HAVE_MPI
    if (mpi_comm == MPI_COMM_WORLD)
      return Dune::StructuredGridFactory<GridType>::createCubeGrid(lowerLeft, upperRight, elements);
#  endif
    DUNE_THROW(InvalidStateException, "ALUgrid Cube construction cannot handle non-world communicators");
  }
};
#endif // HAVE_ALBERTA

#if HAVE_DUNE_SPGRID
template <class ct, int dim, template <int> class Refinement, class Comm>
class StructuredGridFactory<Dune::SPGrid<ct, dim, Refinement, Comm>>
  : public Dune::StructuredGridFactory<Dune::SPGrid<ct, dim, Refinement, Comm>>
{
  typedef Dune::SPGrid<ct, dim, Refinement, Comm> GridType;
  typedef typename GridType::ctype ctype;

public:
  static std::shared_ptr<GridType>
  createCubeGrid(const Dune::FieldVector<ctype, dim>& lowerLeft,
                 const Dune::FieldVector<ctype, dim>& upperRight,
                 const std::array<unsigned int, dim>& elements,
                 std::array<unsigned int, dim> overlap = default_overlap<GridType>(),
                 Dune::MPIHelper::MPICommunicator communicator = Dune::MPIHelper::getCommunicator())
  {
    return Dune::StructuredGridFactory<GridType>::createCubeGrid(
        lowerLeft, upperRight, elements, overlap, communicator);
  }
};

#endif // HAVE_DUNE_SPGRID

template <int dim, class Coords>
class StructuredGridFactory<Dune::YaspGrid<dim, Coords>>
  : public Dune::StructuredGridFactory<Dune::YaspGrid<dim, Coords>>
{
  typedef Dune::YaspGrid<dim, Coords> GridType;
  typedef typename GridType::ctype ctype;

public:
  //! creates a non-periodic Yasp grid
  static std::shared_ptr<GridType>
  createCubeGrid(const Dune::FieldVector<ctype, dim>& lowerLeft,
                 const Dune::FieldVector<ctype, dim>& upperRight,
                 const std::array<unsigned int, dim>& elements_in,
                 std::array<unsigned int, dim> overlap = default_overlap<GridType>(),
                 Dune::MPIHelper::MPICommunicator communicator = Dune::MPIHelper::getCommunicator())
  {
    const auto no_periodic_direction = std::bitset<dim>();
    std::array<int, dim> elements;
    std::copy(elements_in.begin(), elements_in.end(), elements.begin());
    auto overlap_check = overlap;
    overlap_check.fill(overlap[0]);
    for (auto i : Common::value_range(1, dim))
      if (overlap[i] != overlap[0])
        DUNE_THROW(Dune::InvalidStateException, "YaspGrid only supports uniform overlap");
    return std::make_shared<GridType>(lowerLeft, upperRight, elements, no_periodic_direction, overlap[0], communicator);
  }
  static std::shared_ptr<GridType>
  createSimplexGrid(const Dune::FieldVector<ctype, GridType::dimension>&,
                    const Dune::FieldVector<ctype, GridType::dimension>&,
                    const std::array<unsigned int, GridType::dimension>&,
                    Dune::MPIHelper::MPICommunicator = Dune::MPIHelper::getCommunicator())
  {
    DUNE_THROW(NotImplemented, "YaspGrid does not provide Simplex Geometries");
  }
};

} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_STRUCTUREDGRIDFACTORY_HH
