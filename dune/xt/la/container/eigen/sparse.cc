#include <config.h>

#if HAVE_EIGEN

#include "sparse.hh"


template class Dune::XT::LA::EigenRowMajorSparseMatrix<double>;
// template void Dune::XT::LA::EigenRowMajorSparseMatrix<double>::mv(const Dune::XT::LA::EigenDenseVector<double>&,
//                                                                  Dune::XT::LA::EigenDenseVector<double>&) const;
// template void Dune::XT::LA::EigenRowMajorSparseMatrix<double>::mv(const
// Dune::XT::LA::EigenMappedDenseVector<double>&,
//                                                                  Dune::XT::LA::EigenMappedDenseVector<double>&)
//                                                                  const;
// template void Dune::XT::LA::EigenRowMajorSparseMatrix<double>::mv(const Dune::XT::LA::EigenDenseVector<double>&,
//                                                                  Dune::XT::LA::EigenMappedDenseVector<double>&)
//                                                                  const;
// template void Dune::XT::LA::EigenRowMajorSparseMatrix<double>::mv(const
// Dune::XT::LA::EigenMappedDenseVector<double>&,
//                                                                  Dune::XT::LA::EigenDenseVector<double>&) const;


#endif // HAVE_EIGEN
