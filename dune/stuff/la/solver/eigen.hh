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
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/aliases.hh>
#include <dune/stuff/la/container/eigen.hh>

#include "interface.hh"

namespace Dune {
namespace Stuff {
namespace LA {
namespace Solver {

template <class MatrixImp, class VectorImp>
struct IsEigenMV
{
  typedef typename std::is_base_of<Dune::Stuff::LA::Container::EigenMatrixInterface<typename MatrixImp::Traits>,
                                   MatrixImp>::type Mtype;
  typedef typename std::is_base_of<Dune::Stuff::LA::Container::EigenVectorInterface<typename VectorImp::Traits>,
                                   VectorImp>::type Vtype;
  static constexpr bool value = Mtype::value && Vtype::value;
};

//! \attention Slow!
//! \todo Implement via Eigen::CG and identity preconditioner!
template <class MatrixImp, class VectorImp>
class Cg<MatrixImp, VectorImp, typename std::enable_if<IsEigenMV<MatrixImp, VectorImp>::value>::type>
    : public Interface<MatrixImp, VectorImp>
{
public:
  typedef Interface<MatrixImp, VectorImp> BaseType;

  typedef typename BaseType::MatrixType MatrixType;
  typedef typename BaseType::VectorType VectorType;
  typedef typename BaseType::ElementType ElementType;
  typedef typename BaseType::size_type size_type;

  Cg()
  {
    DSC_LOG_DEBUG << "\n" << Dune::Stuff::Common::colorString("WARNING:") << " this solver is believed to be slow! "
                  << std::flush;
  }

  static Dune::ParameterTree createSampleDescription()
  {
    Dune::ParameterTree description;
    description["maxIter"]   = "5000";
    description["precision"] = "1e-12";
    return description;
  } // Dune::ParameterTree createSampleDescription()

  virtual size_type apply(const MatrixType& systemMatrix, const VectorType& rhsVector, VectorType& solutionVector,
                          const size_type maxIter = 5000, const ElementType precision = 1e-12,
                          const Dune::ParameterTree /*description*/ = Dune::ParameterTree()) const
  {
    auto& x_i      = solutionVector.backend();
    const auto& b  = rhsVector.backend();
    const auto& A  = systemMatrix.backend();
    const int cols = A.cols();
    size_type iteration(1);
    ElementType rho(0), rho_prev(1), beta, alpha;
    const ElementType tolerance = precision * precision * b.squaredNorm();
    typedef typename DSLC::EigenDenseVector<typename VectorType::ElementType>::BackendType RealEigenVector;
    RealEigenVector residuum = b - A * x_i;
    RealEigenVector correction_p(cols);
    RealEigenVector correction_q(cols);
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
        return 0;
      rho_prev = rho;
    }
    return 1;
  } // virtual bool apply(...)
}; // class Cg


template <class MatrixImp, class VectorImp>
class CgDiagonal<MatrixImp, VectorImp, typename std::enable_if<IsEigenMV<MatrixImp, VectorImp>::value>::type>
    : public Interface<MatrixImp, VectorImp>
{
public:
  typedef Interface<MatrixImp, VectorImp> BaseType;

  typedef typename BaseType::MatrixType MatrixType;
  typedef typename BaseType::VectorType VectorType;
  typedef typename BaseType::ElementType ElementType;
  typedef typename BaseType::size_type size_type;

  CgDiagonal()
  {
    DSC_LOG_DEBUG << "\n" << Dune::Stuff::Common::colorString("WARNING:", Dune::Stuff::Common::Colors::red)
                  << " this solver is believed to produce utterly wrong results! " << std::flush;
  }

  static Dune::ParameterTree createSampleDescription()
  {
    Dune::ParameterTree description;
    description["maxIter"]   = "5000";
    description["precision"] = "1e-12";
    return description;
  } // Dune::ParameterTree createSampleDescription()

  virtual size_type apply(const MatrixType& systemMatrix, const VectorType& rhsVector, VectorType& solutionVector,
                          const size_type maxIter = 5000, const ElementType precision = 1e-12,
                          const Dune::ParameterTree /*description*/ = Dune::ParameterTree()) const
  {
    typedef ::Eigen::ConjugateGradient<typename MatrixType::BackendType,
                                       ::Eigen::Lower,
                                       ::Eigen::DiagonalPreconditioner<ElementType>> EigenSolverType;
    EigenSolverType eigenSolver(systemMatrix.backend());
    eigenSolver.setMaxIterations(maxIter);
    eigenSolver.setTolerance(precision);
    solutionVector.backend() = eigenSolver.solve(rhsVector.backend());
    return BaseType::translateInfo(eigenSolver.info());
  } // virtual bool apply(...)
}; // class CgDiagonal


template <class MatrixImp, class VectorImp>
class Bicgstab<MatrixImp, VectorImp, typename std::enable_if<IsEigenMV<MatrixImp, VectorImp>::value>::type>
    : public Interface<MatrixImp, VectorImp>
{
public:
  typedef Interface<MatrixImp, VectorImp> BaseType;

  typedef typename BaseType::MatrixType MatrixType;
  typedef typename BaseType::VectorType VectorType;
  typedef typename BaseType::ElementType ElementType;
  typedef typename BaseType::size_type size_type;

  static Dune::ParameterTree createSampleDescription()
  {
    Dune::ParameterTree description;
    description["maxIter"]   = "5000";
    description["precision"] = "1e-12";
    return description;
  } // Dune::ParameterTree createSampleDescription()

  virtual size_type apply(const MatrixType& systemMatrix, const VectorType& rhsVector, VectorType& solutionVector,
                          const size_type maxIter = 5000, const ElementType precision = 1e-12,
                          const Dune::ParameterTree /*description*/ = Dune::ParameterTree()) const
  {
    typedef ::Eigen::BiCGSTAB<typename MatrixType::BackendType, ::Eigen::IdentityPreconditioner> EigenSolverType;
    EigenSolverType eigenSolver(systemMatrix.backend());
    eigenSolver.setMaxIterations(maxIter);
    eigenSolver.setTolerance(precision);
    solutionVector.backend() = eigenSolver.solve(rhsVector.backend());
    return BaseType::translateInfo(eigenSolver.info());
  } // virtual bool apply(...)
}; // class Bicgstab


template <class MatrixImp, class VectorImp>
class BicgstabDiagonal<MatrixImp, VectorImp, typename std::enable_if<IsEigenMV<MatrixImp, VectorImp>::value>::type>
    : public Interface<MatrixImp, VectorImp>
{
public:
  typedef Interface<MatrixImp, VectorImp> BaseType;

  typedef typename BaseType::MatrixType MatrixType;
  typedef typename BaseType::VectorType VectorType;
  typedef typename BaseType::ElementType ElementType;
  typedef typename BaseType::size_type size_type;

  static Dune::ParameterTree createSampleDescription()
  {
    Dune::ParameterTree description;
    description["maxIter"]   = "5000";
    description["precision"] = "1e-12";
    return description;
  } // Dune::ParameterTree createSampleDescription()

  virtual size_type apply(const MatrixType& systemMatrix, const VectorType& rhsVector, VectorType& solutionVector,
                          const size_type maxIter = 5000, const ElementType precision = 1e-12,
                          const Dune::ParameterTree /*description*/ = Dune::ParameterTree()) const
  {
    typedef ::Eigen::BiCGSTAB<typename MatrixType::BackendType, ::Eigen::DiagonalPreconditioner<ElementType>>
        EigenSolverType;
    EigenSolverType eigenSolver(systemMatrix.backend());
    eigenSolver.setMaxIterations(maxIter);
    eigenSolver.setTolerance(precision);
    solutionVector.backend() = eigenSolver.solve(rhsVector.backend());
    return BaseType::translateInfo(eigenSolver.info());
  } // virtual bool apply(...)static_assert
}; // class BicgstabDiagonal


template <class MatrixImp, class VectorImp>
class BicgstabILUT<MatrixImp, VectorImp, typename std::enable_if<IsEigenMV<MatrixImp, VectorImp>::value>::type>
    : public Interface<MatrixImp, VectorImp>
{
public:
  typedef Interface<MatrixImp, VectorImp> BaseType;

  typedef typename BaseType::MatrixType MatrixType;
  typedef typename BaseType::VectorType VectorType;
  typedef typename BaseType::ElementType ElementType;
  typedef typename BaseType::size_type size_type;

  static Dune::ParameterTree createSampleDescription()
  {
    Dune::ParameterTree description;
    description["maxIter"]   = "5000";
    description["precision"] = "1e-12";
    return description;
  } // Dune::ParameterTree createSampleDescription()

  virtual size_type apply(const MatrixType& systemMatrix, const VectorType& rhsVector, VectorType& solutionVector,
                          const size_type maxIter = 5000, const ElementType precision = 1e-12,
                          const Dune::ParameterTree /*description*/ = Dune::ParameterTree()) const
  {
    typedef ::Eigen::BiCGSTAB<typename MatrixType::BackendType, ::Eigen::IncompleteLUT<ElementType>> EigenSolverType;
    EigenSolverType solver(systemMatrix.backend());
    // configure solver and preconditioner
    double solverEps = DSC_CONFIG_GET("algorithm.solver.eps", 1.0e-8);
    solver.setTolerance(solverEps);
    const int solverVerbose = DSC_CONFIG_GET("algorithm.solver.verbose", 0);
    const int maxIterations = DSC_CONFIG_GET("algorithm.solver.maxIterations", solutionVector.size());
    solver.setMaxIterations(maxIterations);
    double dropTol =
        DSC_CONFIG_GET("algorithm.solver.preconditioner.dropTol", Eigen::NumTraits<double>::dummy_precision());
    int fillFactor = DSC_CONFIG_GET("algorithm.solver.preconditioner.fillFactor", 10);
    solver.preconditioner().setDroptol(dropTol);
    solver.preconditioner().setFillfactor(fillFactor);


    solutionVector.backend() = solver.solve(rhsVector.backend());
    return BaseType::translateInfo(solver.info());
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

  static Dune::ParameterTree createSampleDescription()
  {
    return Dune::ParameterTree();
  } // Dune::ParameterTree createSampleDescription()
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

  static Dune::ParameterTree createSampleDescription()
  {
    return Dune::ParameterTree();
  } // Dune::ParameterTree createSampleDescription()
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
