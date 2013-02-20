#ifndef DUNE_STUFF_LA_SOLVER_FASP_HH
#define DUNE_STUFF_LA_SOLVER_FASP_HH

#if HAVE_FASP
#if HAVE_EIGEN

extern "C" {
#include "fasp/fasp_functs.h"
}

#include <dune/stuff/la/container/eigen.hh>

#include "interface.hh"

namespace Dune {
namespace Stuff {
namespace LA {
namespace Solver {


template <class ElementImp>
class AmgFasp<Dune::Stuff::LA::Container::EigenRowMajorSparseMatrix<ElementImp>,
              Dune::Stuff::LA::Container::EigenDenseVector<ElementImp>>
    : public Interface<Dune::Stuff::LA::Container::EigenRowMajorSparseMatrix<ElementImp>,
                       Dune::Stuff::LA::Container::EigenDenseVector<ElementImp>>
{
public:
  typedef Interface<Dune::Stuff::LA::Container::EigenRowMajorSparseMatrix<ElementImp>,
                    Dune::Stuff::LA::Container::EigenDenseVector<ElementImp>> BaseType;

  typedef typename BaseType::MatrixType MatrixType;
  typedef typename BaseType::VectorType VectorType;
  typedef typename BaseType::ElementType ElementType;
  typedef typename BaseType::size_type size_type;

  virtual size_type apply(const MatrixType& _systemMatrix, const VectorType& _rhsVector, VectorType& solutionVector,
                          const size_type maxIter = 5000, const ElementType precision = 1e-12) const
  {
    MatrixType& systemMatrix = const_cast<MatrixType&>(_systemMatrix);
    VectorType& rhsVector    = const_cast<VectorType&>(_rhsVector);

    input_param inparam; // parameters from input files
    itsolver_param itparam; // parameters for itsolver
    AMG_param amgparam; // parameters for AMG
    ILU_param iluparam; // parameters for ILU
    Schwarz_param swzparam;
    fasp_param_init("input.dat", &inparam, &itparam, &amgparam, &iluparam, &swzparam);

    const int row = systemMatrix.rows();
    const int col = systemMatrix.cols();
    const int nnz = systemMatrix.backend().nonZeros();
    int* ja       = systemMatrix.backend().innerIndexPtr();
    int* ia       = systemMatrix.backend().outerIndexPtr();
    dCSRmat A;
    A.row = row;
    A.col = col;
    A.nnz = nnz;
    A.IA  = ia;
    A.JA  = ja;

    dvector b, x;
    b.row = rhsVector.size();
    b.val = rhsVector.backend().data();
    x.row = rhsVector.backend().rows();
    x.val = solutionVector.backend().data();
    return fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itparam, &amgparam);

  } // ... apply(...)

}; // class Fasp

} // namespace Solver
} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // HAVE_EIGEN
#endif // HAVE_FASP

#endif // DUNE_STUFF_LA_SOLVER_FASP_HH
