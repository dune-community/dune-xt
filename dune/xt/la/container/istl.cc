// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)
//   Tobias Leibner  (2017)

#include <config.h>

#if HAVE_DUNE_ISTL

#include "istl.hh"


template class Dune::XT::LA::IstlDenseVector<double>;
template class Dune::XT::LA::IstlRowMajorSparseMatrix<double>;
// template std::ostream& operator<<(std::ostream&, const Dune::XT::LA::IstlRowMajorSparseMatrix<double>&);


#endif // HAVE_DUNE_ISTL
