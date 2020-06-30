// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014, 2016 - 2017)
//   René Fritze     (2014 - 2016, 2018 - 2019)
//   Tobias Leibner  (2014, 2019 - 2020)

#ifndef DUNE_XT_COMMON_PARALLEL_HELPER_HH
#define DUNE_XT_COMMON_PARALLEL_HELPER_HH

#include <type_traits>
#include <cassert>

#include <dune/common/version.hh>
#if DUNE_VERSION_GTE(DUNE_COMMON, 2, 7)
#  include <dune/common/parallel/communication.hh>
#else
#  include <dune/common/parallel/collectivecommunication.hh>
#endif

#include <dune/istl/paamg/pinfo.hh>

namespace Dune {
namespace XT {

//! marker for sequential in MPI-enabled solver stuffs
struct SequentialCommunication : public Dune::Amg::SequentialInformation
{};

template <class GridCommImp>
struct UseParallelCommunication
{
#if HAVE_MPI
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

} // namespace XT
} // namespace Dune

#endif // DUNE_XT_COMMON_PARALLEL_HELPER_HH
