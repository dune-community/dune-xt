#ifndef DUNE_STUFF_LA_SOLVER_EIGEN_HH
#define DUNE_STUFF_LA_SOLVER_EIGEN_HH

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#elif defined(HAVE_CONFIG_H)
#include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#if HAVE_EIGEN

#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>

#include <dune/stuff/common/color.hh>
#include <dune/stuff/la/container/eigen.hh>

#include "interface.hh"

namespace Dune {
namespace Stuff {
namespace LA {
namespace Solver {

//! \attention Slow!
//! \todo Implement via Eigen::CG and identity preconditioner!
template <class ElementImp>
class Cg<Dune::Stuff::LA::Container::EigenRowMajorSparseMatrix<ElementImp>,
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

  Cg()
  {
    if (!Dune::Stuff::Common::Logger().created())
      Dune::Stuff::Common::Logger().create(Dune::Stuff::Common::LOG_CONSOLE | Dune::Stuff::Common::LOG_DEBUG);
    Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().debug();
    debug << "\n" << Dune::Stuff::Common::colorString("WARNING:") << " this solver is believed to be slow! "
          << std::flush;
  }

  virtual bool apply(const MatrixType& systemMatrix, const VectorType& rhsVector, VectorType& solutionVector,
                     const size_type maxIter = 5000, const ElementType precision = 1e-12) const
  {
    auto& x_i      = solutionVector.backend();
    const auto& b  = rhsVector.backend();
    const auto& A  = systemMatrix.backend();
    const int cols = A.cols();
    size_type iteration(1);
    ElementType rho(0), rho_prev(1), beta, alpha;
    const ElementType tolerance               = precision * precision * b.squaredNorm();
    typename VectorType::BackendType residuum = b - A * x_i;
    typename VectorType::BackendType correction_p(cols);
    typename VectorType::BackendType correction_q(cols);
    rho = residuum.squaredNorm();
    while (iteration <= maxIter) {
      if (iteration == 1) {
        correction_p = residuum;
      } else {
        beta         = rho / rho_prev;
        correction_p = residuum + beta * correction_p;
      }
      correction_q.noalias() = A * correction_p;
      alpha = rho / correction_p.dot(correction_q);
      x_i += alpha * correction_p;
      residuum -= alpha * correction_q;
      rho = residuum.squaredNorm();
      if (rho < tolerance)
        return true;
      rho_prev = rho;
    }
    return false;
  } // virtual bool apply(...)
}; // class Cg


template <class ElementImp>
class CgDiagonal<Dune::Stuff::LA::Container::EigenRowMajorSparseMatrix<ElementImp>,
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

  CgDiagonal()
  {
    if (!Dune::Stuff::Common::Logger().created())
      Dune::Stuff::Common::Logger().create(Dune::Stuff::Common::LOG_CONSOLE | Dune::Stuff::Common::LOG_DEBUG);
    Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().debug();
    debug << "\n" << Dune::Stuff::Common::colorString("WARNING:", Dune::Stuff::Common::Colors::red)
          << " this solver is believed to produce utterly wrong results! " << std::flush;
  }

  virtual bool apply(const MatrixType& systemMatrix, const VectorType& rhsVector, VectorType& solutionVector,
                     const size_type maxIter = 5000, const ElementType precision = 1e-12) const
  {
    typedef ::Eigen::ConjugateGradient<typename MatrixType::BackendType,
                                       ::Eigen::Lower,
                                       ::Eigen::DiagonalPreconditioner<ElementType>> EigenSolverType;
    EigenSolverType eigenSolver(systemMatrix.backend());
    eigenSolver.setMaxIterations(maxIter);
    eigenSolver.setTolerance(precision);
    solutionVector.backend() = eigenSolver.solve(rhsVector.backend());
    const ::Eigen::ComputationInfo info = eigenSolver.info();
    return (info == ::Eigen::Success);
  } // virtual bool apply(...)
}; // class CgDiagonal


template <class ElementImp>
class Bicgstab<Dune::Stuff::LA::Container::EigenRowMajorSparseMatrix<ElementImp>,
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

  Bicgstab()
  {
  }

  virtual bool apply(const MatrixType& systemMatrix, const VectorType& rhsVector, VectorType& solutionVector,
                     const size_type maxIter = 5000, const ElementType precision = 1e-12) const
  {
    typedef ::Eigen::BiCGSTAB<typename MatrixType::BackendType, ::Eigen::IdentityPreconditioner> EigenSolverType;
    EigenSolverType eigenSolver(systemMatrix.backend());
    eigenSolver.setMaxIterations(maxIter);
    eigenSolver.setTolerance(precision);
    solutionVector.backend() = eigenSolver.solve(rhsVector.backend());
    const ::Eigen::ComputationInfo info = eigenSolver.info();
    return (info == ::Eigen::Success);
  } // virtual bool apply(...)
}; // class Bicgstab


template <class ElementImp>
class BicgstabDiagonal<Dune::Stuff::LA::Container::EigenRowMajorSparseMatrix<ElementImp>,
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

  BicgstabDiagonal()
  {
  }

  virtual bool apply(const MatrixType& systemMatrix, const VectorType& rhsVector, VectorType& solutionVector,
                     const size_type maxIter = 5000, const ElementType precision = 1e-12) const
  {
    typedef ::Eigen::BiCGSTAB<typename MatrixType::BackendType, ::Eigen::DiagonalPreconditioner<ElementType>>
        EigenSolverType;
    EigenSolverType eigenSolver(systemMatrix.backend());
    eigenSolver.setMaxIterations(maxIter);
    eigenSolver.setTolerance(precision);
    solutionVector.backend() = eigenSolver.solve(rhsVector.backend());
    const ::Eigen::ComputationInfo info = eigenSolver.info();
    return (info == ::Eigen::Success);
  } // virtual bool apply(...)static_assert
}; // class BicgstabDiagonal


template <class ElementImp>
class BicgstabILUT<Dune::Stuff::LA::Container::EigenRowMajorSparseMatrix<ElementImp>,
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

  BicgstabILUT()
  {
  }

  virtual bool apply(const MatrixType& systemMatrix, const VectorType& rhsVector, VectorType& solutionVector,
                     const size_type maxIter = 5000, const ElementType precision = 1e-12) const
  {
    typedef ::Eigen::BiCGSTAB<typename MatrixType::BackendType, ::Eigen::IncompleteLUT<ElementType>> EigenSolverType;
    EigenSolverType eigenSolver(systemMatrix.backend());
    eigenSolver.setMaxIterations(maxIter);
    eigenSolver.setTolerance(precision);
    solutionVector.backend() = eigenSolver.solve(rhsVector.backend());
    const ::Eigen::ComputationInfo info = eigenSolver.info();
    return (info == ::Eigen::Success);
  } // virtual bool apply(...)
}; // class BicgstabILUT


template <class ElementImp>
class SimplicialLLT<Dune::Stuff::LA::Container::EigenRowMajorSparseMatrix<ElementImp>,
                    Dune::Stuff::LA::Container::EigenDenseVector<ElementImp>>
    : public SolverNotImplementedForThisMatrixVectorCombination<Dune::Stuff::LA::Container::
                                                                    EigenRowMajorSparseMatrix<ElementImp>,
                                                                Dune::Stuff::LA::Container::
                                                                    EigenDenseVector<ElementImp>>
{
public:
  typedef SolverNotImplementedForThisMatrixVectorCombination<Dune::Stuff::LA::Container::
                                                                 EigenRowMajorSparseMatrix<ElementImp>,
                                                             Dune::Stuff::LA::Container::EigenDenseVector<ElementImp>>
      BaseType;

  SimplicialLLT()
    : BaseType("\nERROR: only implemented for eigen matrices of type 'EigenColMajorSparseMatrix'!")
  {
  }
}; // class SimplicialLLT


// template< class ElementImp >
// class SimplicialLLT< Dune::Stuff::LA::Container::EigenColMajorSparseMatrix< ElementImp >,
//                    Dune::Stuff::LA::Container::EigenDenseVector< ElementImp > >
//  : public Interface< Dune::Stuff::LA::Container::EigenColMajorSparseMatrix< ElementImp >,
//                      Dune::Stuff::LA::Container::EigenDenseVector< ElementImp > >
//{
// public:
//  typedef Interface<  Dune::Stuff::LA::Container::EigenColMajorSparseMatrix< ElementImp >,
//                      Dune::Stuff::LA::Container::EigenDenseVector< ElementImp > >
//      BaseType;

//  typedef typename BaseType::MatrixType MatrixType;

//  typedef typename BaseType::VectorType VectorType;

//  typedef typename BaseType::ElementType ElementType;

//  typedef typename BaseType::size_type size_type;

//  SimplicialLLT()
//  {}

//  virtual bool apply(const MatrixType& systemMatrix,
//                     const VectorType& rhsVector,
//                     VectorType& solutionVector,
//                     const size_type /*maxIter*/ = 0,
//                     const ElementType /*precision*/ = 0) const
//  {
//    typedef ::Eigen::SimplicialLLT< typename MatrixType::BackendType, ::Eigen::Lower > EigenSolverType;
//    EigenSolverType eigenSolver(systemMatrix.base());
//    solutionVector.base() = eigenSolver.solve(rhsVector.base());
//    const ::Eigen::ComputationInfo info = eigenSolver.info();
//    return (info == ::Eigen::Success);
//  } // virtual bool apply(...)
//}; // class SimplicialLLT


template <class ElementImp>
class SimplicialLDLT<Dune::Stuff::LA::Container::EigenRowMajorSparseMatrix<ElementImp>,
                     Dune::Stuff::LA::Container::EigenDenseVector<ElementImp>>
    : public SolverNotImplementedForThisMatrixVectorCombination<Dune::Stuff::LA::Container::
                                                                    EigenRowMajorSparseMatrix<ElementImp>,
                                                                Dune::Stuff::LA::Container::
                                                                    EigenDenseVector<ElementImp>>
{
public:
  typedef SolverNotImplementedForThisMatrixVectorCombination<Dune::Stuff::LA::Container::
                                                                 EigenRowMajorSparseMatrix<ElementImp>,
                                                             Dune::Stuff::LA::Container::EigenDenseVector<ElementImp>>
      BaseType;

  SimplicialLDLT()
    : BaseType("\nERROR: only implemented for eigen matrices of type 'EigenColMajorSparseMatrix'!")
  {
  }
}; // class SimplicialLDLT


// template< class ElementImp >
// class SimplicialLDLT< Dune::Stuff::LA::Container::EigenColMajorSparseMatrix< ElementImp >,
//                      Dune::Stuff::LA::Container::EigenDenseVector< ElementImp > >
//  : public Interface< Dune::Stuff::LA::Container::EigenColMajorSparseMatrix< ElementImp >,
//                      Dune::Stuff::LA::Container::EigenDenseVector< ElementImp > >
//{
// public:
//  typedef Interface<  Dune::Stuff::LA::Container::EigenColMajorSparseMatrix< ElementImp >,
//                      Dune::Stuff::LA::Container::EigenDenseVector< ElementImp > >
//      BaseType;

//  typedef typename BaseType::MatrixType MatrixType;

//  typedef typename BaseType::VectorType VectorType;

//  typedef typename BaseType::ElementType ElementType;

//  typedef typename BaseType::size_type size_type;

//  SimplicialLLT()
//  {}

//  virtual bool apply(const MatrixType& systemMatrix,
//                     const VectorType& rhsVector,
//                     VectorType& solutionVector,
//                     const size_type /*maxIter*/ = 0,
//                     const ElementType /*precision*/ = 0) const
//  {
//    typedef ::Eigen::SimplicialLDLT< typename MatrixType::BackendType, ::Eigen::Lower > EigenSolverType;
//    EigenSolverType eigenSolver(systemMatrix.base());
//    solutionVector.base() = eigenSolver.solve(rhsVector.base());
//    const ::Eigen::ComputationInfo info = eigenSolver.info();
//    return (info == ::Eigen::Success);
//  } // virtual bool apply(...)
//}; // class SimplicialLDLT


} // namespace Solver
} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // HAVE_EIGEN

#endif // DUNE_STUFF_LA_SOLVER_EIGEN_HH
