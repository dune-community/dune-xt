// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2014, 2016 - 2018)
//   René Fritze     (2013, 2015 - 2016, 2018 - 2019)
//   Tobias Leibner  (2018 - 2020)

#ifndef DUNE_XT_COMMON_PARALLEL_COMMUNICATOR_HH
#define DUNE_XT_COMMON_PARALLEL_COMMUNICATOR_HH

#include <dune/common/parallel/communicator.hh>

#if !HAVE_MPI

// defines No_Comm + CollectiveCommunication
#  include <dune/common/parallel/communication.hh>

namespace Dune {


// dune-common is missing comparison operators for No_Comm
bool operator==(const No_Comm& /*lhs*/, const No_Comm& /*rhs*/);
bool operator!=(const No_Comm& lhs, const No_Comm& rhs);

// CollectiveCommunication<No_Comm> is also missing a comparison operator (CollectiveCommunication<MPI_Comm> does not
// have a comparison operator either, but it has a conversion operator to MPI_Comm which does have operator==)
bool operator==(const CollectiveCommunication<No_Comm>& /*lhs*/, const CollectiveCommunication<No_Comm>& /*rhs*/);
bool operator!=(const CollectiveCommunication<No_Comm>& lhs, const CollectiveCommunication<No_Comm>& rhs);


} // namespace Dune

#endif // !HAVE_MPI

#endif // DUNE_XT_COMMON_PARALLEL_COMMUNICATOR_HH
