#include <config.h>

#if HAVE_DUNE_ISTL

#include "istl.hh"


template class Dune::XT::LA::Solver<Dune::XT::LA::IstlRowMajorSparseMatrix<double>>;


#endif // HAVE_DUNE_ISTL
