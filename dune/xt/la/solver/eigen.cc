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

#if HAVE_EIGEN

#  include "eigen.hh"


template class Dune::XT::LA::Solver<Dune::XT::LA::EigenDenseMatrix<double>>;
// template void
// Dune::XT::LA::Solver<Dune::XT::LA::EigenDenseMatrix<double>>::apply(const Dune::XT::LA::EigenDenseVector<double>&,
//                                                                    Dune::XT::LA::EigenDenseVector<double>&) const;
// template void Dune::XT::LA::Solver<Dune::XT::LA::EigenDenseMatrix<double>>::apply(
//    const Dune::XT::LA::EigenDenseVector<double>&, Dune::XT::LA::EigenDenseVector<double>&, const std::string&) const;
// template void
// Dune::XT::LA::Solver<Dune::XT::LA::EigenDenseMatrix<double>>::apply(const Dune::XT::LA::EigenDenseVector<double>&,
//                                                                    Dune::XT::LA::EigenDenseVector<double>&,
//                                                                    const Dune::XT::Common::Configuration&) const;

template class Dune::XT::LA::Solver<Dune::XT::LA::EigenRowMajorSparseMatrix<double>>;
// template void Dune::XT::LA::Solver<Dune::XT::LA::EigenRowMajorSparseMatrix<double>>::apply(
//    const Dune::XT::LA::EigenDenseVector<double>&, Dune::XT::LA::EigenDenseVector<double>&) const;
// template void Dune::XT::LA::Solver<Dune::XT::LA::EigenRowMajorSparseMatrix<double>>::apply(
//    const Dune::XT::LA::EigenDenseVector<double>&, Dune::XT::LA::EigenDenseVector<double>&, const std::string&) const;
// template void Dune::XT::LA::Solver<Dune::XT::LA::EigenRowMajorSparseMatrix<double>>::apply(
//    const Dune::XT::LA::EigenDenseVector<double>&,
//    Dune::XT::LA::EigenDenseVector<double>&,
//    const Dune::XT::Common::Configuration&) const;


#endif // HAVE_EIGEN
