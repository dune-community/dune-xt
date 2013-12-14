// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff/
// Copyright Holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_LA_SOLVER_STUFF_HH
#define DUNE_STUFF_LA_SOLVER_STUFF_HH

//#if HAVE_DUNE_ISTL

#include <type_traits>

#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/paamg/amg.hh>

#include <dune/stuff/la/container/istl.hh>

#include "interface.hh"

namespace Dune {
namespace Stuff {
namespace LA {

template <class M, class V>
struct IsIstlMV
{
  typedef typename std::is_same<IstlRowMajorSparseMatrix<typename M::ScalarType>, M>::type Mtype;
  typedef typename std::is_same<IstlDenseVector<typename V::ScalarType>, V>::type Vtype;
  static constexpr bool value = Mtype::value && Vtype::value;
};


// template< class MatrixImp, class VectorImp >
// class CgSolver< MatrixImp, VectorImp, typename std::enable_if< IsIstlMV< MatrixImp, VectorImp >::value >::type >
//  : public SolverInterface< MatrixImp, VectorImp >
//{
//  typedef SolverInterface< MatrixImp, VectorImp > BaseType;
// public:
//  typedef typename BaseType::MatrixType MatrixType;
//  typedef typename BaseType::VectorType VectorType;
//  typedef typename BaseType::ScalarType ScalarType;

//  virtual size_t apply(const MatrixType& lhs,
//                       const VectorType& rhs,
//                       VectorType& solution,
//                       const Dune::ParameterTree description = BaseType::defaultIterativeSettings()) const
//  {
//    assert(false);
//  } // ... apply(...)
//}; // class CgSolver


// template< class MatrixImp, class VectorImp >
// class CgDiagonalSolver< MatrixImp, VectorImp, typename std::enable_if< IsIstlMV< MatrixImp, VectorImp >::value
// >::type >
//  : public SolverInterface< MatrixImp, VectorImp >
//{
//  typedef SolverInterface< MatrixImp, VectorImp > BaseType;
// public:
//  typedef typename BaseType::MatrixType MatrixType;
//  typedef typename BaseType::VectorType VectorType;
//  typedef typename BaseType::ScalarType ScalarType;

//  virtual size_t apply(const MatrixType& lhs,
//                       const VectorType& rhs,
//                       VectorType& solution,
//                       const Dune::ParameterTree description = BaseType::defaultIterativeSettings()) const
//  {
//    assert(false);
//  } // ... apply(...)
//}; // class CgDiagonalSolver


// template< class MatrixImp, class VectorImp >
// class BicgstabSolver< MatrixImp, VectorImp, typename std::enable_if< IsIstlMV< MatrixImp, VectorImp >::value >::type
// >
//  : public SolverInterface< MatrixImp, VectorImp >
//{
//  typedef SolverInterface< MatrixImp, VectorImp > BaseType;
// public:
//  typedef typename BaseType::MatrixType MatrixType;
//  typedef typename BaseType::VectorType VectorType;
//  typedef typename BaseType::ScalarType ScalarType;

//  virtual size_t apply(const MatrixType& lhs,
//                       const VectorType& rhs,
//                       VectorType& solution,
//                       const Dune::ParameterTree description = BaseType::defaultIterativeSettings()) const
//  {
//    assert(false);
//  } // ... apply(...)
//}; // class BicgstabSolver


// template< class MatrixImp, class VectorImp >
// class BicgstabDiagonalSolver< MatrixImp, VectorImp,
//                              typename std::enable_if< IsIstlMV< MatrixImp, VectorImp >::value >::type >
//  : public SolverInterface< MatrixImp, VectorImp >
//{
//  typedef SolverInterface< MatrixImp, VectorImp > BaseType;
// public:
//  typedef typename BaseType::MatrixType MatrixType;
//  typedef typename BaseType::VectorType VectorType;
//  typedef typename BaseType::ScalarType ScalarType;

//  virtual size_t apply(const MatrixType& lhs,
//                       const VectorType& rhs,
//                       VectorType& solution,
//                       const Dune::ParameterTree description = BaseType::defaultIterativeSettings()) const
//  {
//    assert(false);
//  } // ... apply(...)
//}; // class BicgstabDiagonalSolver


template <class MatrixImp, class VectorImp>
class BicgstabILUTSolver<MatrixImp, VectorImp, typename std::enable_if<IsIstlMV<MatrixImp, VectorImp>::value>::type>
    : public SolverInterface<MatrixImp, VectorImp>
{
  typedef SolverInterface<MatrixImp, VectorImp> BaseType;

public:
  typedef typename BaseType::MatrixType MatrixType;
  typedef typename BaseType::VectorType VectorType;
  typedef typename BaseType::ScalarType ScalarType;

  static Dune::ParameterTree defaultSettings()
  {
    Dune::ParameterTree description                 = BaseType::defaultIterativeSettings();
    description["precision"]                        = "1.0e-8";
    description["preconditioner.iterations"]        = "10";
    description["preconditioner.relaxation_factor"] = "1.0";
    return description;
  } // Dune::ParameterTree defaultSettings()

  virtual size_t apply(const MatrixType& lhs, const VectorType& rhs, VectorType& solution,
                       const Dune::ParameterTree description = BaseType::defaultIterativeSettings()) const
  {
    typedef MatrixAdapter<typename MatrixType::BackendType,
                          typename VectorType::BackendType,
                          typename VectorType::BackendType> MatrixOperatorType;
    MatrixOperatorType matrix_operator(*(lhs.backend_));
    typedef SeqILUn<typename MatrixType::BackendType,
                    typename VectorType::BackendType,
                    typename VectorType::BackendType> PreconditionerType;
    PreconditionerType preconditioner(
        *(lhs.backend_),
        description.get("preconditioner.iterations",
                        defaultSettings().template get<size_t>("preconditioner.iterations")),
        description.get("preconditioner.relaxation_factor",
                        defaultSettings().template get<double>("preconditioner.relaxation_factor")));
    typedef BiCGSTABSolver<typename VectorType::BackendType> SolverType;
    SolverType solver(matrix_operator,
                      preconditioner,
                      description.get("precision", defaultSettings().template get<double>("precision")),
                      description.get("maxIter", defaultSettings().template get<size_t>("maxIter")),
                      description.get("verbose", defaultSettings().template get<int>("verbose")));
    Dune::InverseOperatorResult stat;
    solver.apply(*(solution.backend_), *(rhs.backend_), stat);
    if (stat.converged)
      return 0;
    else
      return 3;
  } // ... apply(...)
}; // class BicgstabILUTSolver


// template< class MatrixImp, class VectorImp >
// class DirectSuperLuSolver< MatrixImp, VectorImp,
//                           typename std::enable_if< IsIstlMV< MatrixImp, VectorImp >::value >::type >
//  : public SolverInterface< MatrixImp, VectorImp >
//{
//  typedef SolverInterface< MatrixImp, VectorImp > BaseType;
// public:
//  typedef typename BaseType::MatrixType MatrixType;
//  typedef typename BaseType::VectorType VectorType;
//  typedef typename BaseType::ScalarType ScalarType;

//  DirectSuperLuSolver()
//  {}

//  virtual size_t apply(const MatrixType& lhs,
//                       const VectorType& rhs,
//                       VectorType& solution,
//                       const Dune::ParameterTree description = BaseType::defaultIterativeSettings()) const
//  {
//    assert(false);
//  } // ... apply(...)
//}; // class DirectSuperLuSolver


template <class MatrixImp, class VectorImp>
class AmgSolver<MatrixImp, VectorImp, typename std::enable_if<IsIstlMV<MatrixImp, VectorImp>::value>::type>
    : public SolverInterface<MatrixImp, VectorImp>
{
  typedef SolverInterface<MatrixImp, VectorImp> BaseType;

public:
  typedef typename BaseType::MatrixType MatrixType;
  typedef typename BaseType::VectorType VectorType;
  typedef typename BaseType::ScalarType ScalarType;

  static Dune::ParameterTree defaultSettings()
  {
    return BaseType::defaultIterativeSettings();
  } // Dune::ParameterTree defaultSettings()

  virtual size_t apply(const MatrixType& lhs, const VectorType& rhs, VectorType& solution,
                       const Dune::ParameterTree description = BaseType::defaultIterativeSettings()) const
  {
    typedef MatrixAdapter<typename MatrixType::BackendType,
                          typename VectorType::BackendType,
                          typename VectorType::BackendType> MatrixOperatorType;
    MatrixOperatorType matrix_operator(*(lhs.backend_));
    typedef SeqILU0<typename MatrixType::BackendType,
                    typename VectorType::BackendType,
                    typename VectorType::BackendType> PreconditionerType;
    typedef Dune::Amg::AMG<MatrixOperatorType, typename VectorType::BackendType, PreconditionerType> SmootherType;
    Dune::SeqScalarProduct<typename VectorType::BackendType> scalar_product;
    typedef typename Dune::Amg::SmootherTraits<PreconditionerType>::Arguments SmootherArgs;
    SmootherArgs smootherArgs;
    smootherArgs.iterations       = 1;
    smootherArgs.relaxationFactor = 1;
    Dune::Amg::Parameters params(15, 2000);
    params.setDefaultValuesAnisotropic(2); // <- dim
    typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<typename MatrixType::BackendType,
                                                                      Dune::Amg::FirstDiagonal>> AmgCriterion;
    AmgCriterion amg_criterion(params);
    SmootherType smoother(matrix_operator, amg_criterion, smootherArgs);
    typedef BiCGSTABSolver<typename VectorType::BackendType> SolverType;
    SolverType solver(matrix_operator,
                      scalar_product,
                      smoother,
                      description.get("precision", defaultSettings().template get<double>("precision")),
                      description.get("maxIter", defaultSettings().template get<size_t>("maxIter")),
                      description.get("verbose", defaultSettings().template get<int>("verbose")));
    Dune::InverseOperatorResult stat;
    solver.apply(*(solution.backend_), *(rhs.backend_), stat);
    if (stat.converged)
      return 0;
    else
      return 3;
  } // ... apply(...)
}; // class AmgSolver


} // namespace LA
} // namespace Stuff
} // namespace Dune

//#endif // HAVE_DUNE_ISTL

#endif // DUNE_STUFF_LA_SOLVER_STUFF_HH
