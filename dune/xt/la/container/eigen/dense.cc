#include <config.h>

#if HAVE_EIGEN

#include "dense.hh"


template class Dune::XT::LA::EigenDenseVector<double>;
template class Dune::XT::LA::EigenMappedDenseVector<double>;
template class Dune::XT::LA::EigenDenseMatrix<double>;
// template void Dune::XT::LA::EigenDenseMatrix<double>::mv(const Dune::XT::LA::EigenDenseVector<double>&,
//                                                         Dune::XT::LA::EigenDenseVector<double>&) const;
// template void Dune::XT::LA::EigenDenseMatrix<double>::mv(const Dune::XT::LA::EigenMappedDenseVector<double>&,
//                                                         Dune::XT::LA::EigenMappedDenseVector<double>&) const;
// template void Dune::XT::LA::EigenDenseMatrix<double>::mv(const Dune::XT::LA::EigenDenseVector<double>&,
//                                                         Dune::XT::LA::EigenMappedDenseVector<double>&) const;
// template void Dune::XT::LA::EigenDenseMatrix<double>::mv(const Dune::XT::LA::EigenMappedDenseVector<double>&,
//                                                         Dune::XT::LA::EigenDenseVector<double>&) const;


#endif // HAVE_EIGEN
