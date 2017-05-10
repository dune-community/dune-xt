#include <config.h>

#if HAVE_EIGEN

#include "eigen.hh"


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
