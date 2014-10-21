// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_COMMON_PARALLEL_HELPER_HH
#define DUNE_STUFF_COMMON_PARALLEL_HELPER_HH

#include <type_traits>

#include <dune/common/parallel/collectivecommunication.hh>

#if HAVE_DUNE_ISTL
#include <dune/istl/paamg/pinfo.hh>
#endif

namespace Dune {
namespace Stuff {


//! marker for sequential in MPI-enabled solver stuffs
struct SequentialCommunication
#if HAVE_DUNE_ISTL
    : public Dune::Amg::SequentialInformation
#endif
{
};


template <class GridCommImp>
struct UseParallelCommunication
{
#if HAVE_MPI && HAVE_DUNE_ISTL
  static constexpr bool value = std::is_same<GridCommImp, CollectiveCommunication<MPI_Comm>>::value;
#else
  static constexpr bool value = false;
#endif
};

/**
 * \brief calls MPI_Abort if enable-parallel, noop otherwise
 * \returns MPI_Abort if enable-parallel, 1 otherwise
 **/
int abort_all_mpi_processes();

} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_COMMON_PARALLEL_HELPER_HH
