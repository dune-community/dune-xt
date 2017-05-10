#include <config.h>

#if HAVE_DUNE_ISTL

#include "istl.hh"


template class Dune::XT::LA::IstlDenseVector<double>;
template class Dune::XT::LA::IstlRowMajorSparseMatrix<double>;
// template std::ostream& operator<<(std::ostream&, const Dune::XT::LA::IstlRowMajorSparseMatrix<double>&);


#endif // HAVE_DUNE_ISTL
