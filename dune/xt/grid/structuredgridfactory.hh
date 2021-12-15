// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017, 2020)
//   René Fritze     (2014 - 2016, 2018 - 2020)
//   Sven Kaulmann   (2014)
//   Tobias Leibner  (2014, 2016, 2018, 2020)

#ifndef DUNE_XT_GRID_STRUCTUREDGRIDFACTORY_HH
#define DUNE_XT_GRID_STRUCTUREDGRIDFACTORY_HH

#include <memory>

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
#include <dune/xt/common/logging.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/common/unused.hh>

namespace Dune::XT::Grid {


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
  using ctype = typename GridType::ctype;

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
  using GridType = UGGrid<dim>;
  using ctype = typename GridType::ctype;

public:
  static std::shared_ptr<GridType>
  createCubeGrid(const Dune::FieldVector<ctype, GridType::dimension>& lowerLeft,
                 const Dune::FieldVector<ctype, GridType::dimension>& upperRight,
                 const std::array<unsigned int, GridType::dimension>& elements,
                 std::array<unsigned int, GridType::dimension> = default_overlap<GridType>(),
                 DXTC_MPI_ONLY Dune::MPIHelper::MPICommunicator mpi_comm = Dune::MPIHelper::getCommunicator())
  {
    if (Dune::MPIHelper::isFake) // NOLINT(readability-implicit-bool-conversion)
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
                    DXTC_MPI_ONLY Dune::MPIHelper::MPICommunicator mpi_comm = Dune::MPIHelper::getCommunicator())
  {
    if (Dune::MPIHelper::isFake) // NOLINT(readability-implicit-bool-conversion)
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
  using GridType = AlbertaGrid<dim, dim>;
  using ctype = typename GridType::ctype;

public:
  static std::shared_ptr<GridType>
  createCubeGrid(const Dune::FieldVector<ctype, GridType::dimension>& lowerLeft,
                 const Dune::FieldVector<ctype, GridType::dimension>& upperRight,
                 const std::array<unsigned int, GridType::dimension>& elements,
                 std::array<unsigned int, GridType::dimension> overlap = default_overlap<GridType>(),
                 DXTC_MPI_ONLY Dune::MPIHelper::MPICommunicator mpi_comm = Dune::MPIHelper::getCommunicator())
  {
    bool warn = false;
    unsigned int sum = 0;
    for (auto i : Common::value_range(1, dim)) {
      warn = warn || (overlap[i] > 1);
      sum += overlap[i];
    }
    if (warn || sum >= dim)
      DXTC_LOG_INFO << "Ignoring non-default overlap for alberta cube creation";
    if (Dune::MPIHelper::isFake) // NOLINT(readability-implicit-bool-conversion)
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
                    DXTC_MPI_ONLY Dune::MPIHelper::MPICommunicator mpi_comm = Dune::MPIHelper::getCommunicator())
  {
    if (Dune::MPIHelper::isFake) // NOLINT(readability-implicit-bool-conversion)
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
  using GridType = ALUGrid<dim, dim_world, Dune::cube, refineType, Comm>;
  using ctype = typename GridType::ctype;

public:
  static std::shared_ptr<GridType>
  createCubeGrid(const Dune::FieldVector<ctype, GridType::dimension>& lowerLeft,
                 const Dune::FieldVector<ctype, GridType::dimension>& upperRight,
                 const std::array<unsigned int, GridType::dimension>& elements,
                 std::array<unsigned int, GridType::dimension> = default_overlap<GridType>(),
                 DXTC_MPI_ONLY Dune::MPIHelper::MPICommunicator mpi_comm = Dune::MPIHelper::getCommunicator())
  {
    if (Dune::MPIHelper::isFake) // NOLINT(readability-implicit-bool-conversion)
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
  using GridType = Dune::SPGrid<ct, dim, Refinement, Comm>;
  using ctype = typename GridType::ctype;

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
  using GridType = Dune::YaspGrid<dim, Coords>;
  using ctype = typename GridType::ctype;

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

} // namespace Dune::XT::Grid

#endif // DUNE_XT_GRID_STRUCTUREDGRIDFACTORY_HH
