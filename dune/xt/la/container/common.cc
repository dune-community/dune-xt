#include <config.h>

#include "common.hh"


template class Dune::XT::LA::CommonDenseVector<double>;
template class Dune::XT::LA::CommonDenseMatrix<double>;
template class Dune::XT::LA::CommonSparseMatrix<double>;
// template void Dune::XT::LA::CommonSparseMatrix<double>::mv(const Dune::XT::LA::CommonDenseVector<double>&,
//                                                           Dune::XT::LA::CommonDenseVector<double>) const;
