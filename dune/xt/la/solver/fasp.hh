// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2014, 2016 - 2017)
//   Rene Milk       (2013, 2015 - 2016, 2018)
//   Tobias Leibner  (2017)

#ifndef DUNE_XT_LA_SOLVER_FASP_HH
#define DUNE_XT_LA_SOLVER_FASP_HH

#if HAVE_FASP
#  if HAVE_EIGEN

extern "C"
{
#    include "fasp_functs.h"
}

#    include <dune/xt/la/container/eigen.hh>

#    include "interface.hh"

namespace Dune {
namespace XT {
namespace LA {

template <class ElementImp>
class AmgSolver<Dune::XT::LA::EigenRowMajorSparseMatrix<ElementImp>, Dune::XT::LA::EigenDenseVector<ElementImp>>
  : public SolverInterface<Dune::XT::LA::EigenRowMajorSparseMatrix<ElementImp>,
                           Dune::XT::LA::EigenDenseVector<ElementImp>>
{
public:
  typedef SolverInterface<Dune::XT::LA::EigenRowMajorSparseMatrix<ElementImp>,
                          Dune::XT::LA::EigenDenseVector<ElementImp>>
      BaseType;

  typedef typename BaseType::MatrixType MatrixType;
  typedef typename BaseType::VectorType VectorType;
  typedef typename BaseType::ScalarType ScalarType;

  static Dune::ParameterTree defaultSettings()
  {
    Dune::ParameterTree description = BaseType::defaultIterativeSettings();
    // these parameters were taken from the init.dat that Ludmil gave me...
    description["input_param.print_level"] = "3";
    description["input_param.output_type"] = "0";
    description["input_param.workdir"] = "./";
    description["input_param.problem_num"] = "14 ";
    description["input_param.solver_type"] = "1";
    description["input_param.precond_type"] = "2";
    description["input_param.stop_type"] = "1";
    description["input_param.itsolver_tol"] = "1.000000e-10";
    description["input_param.itsolver_maxit"] = "1000";
    description["input_param.restart"] = "20";
    description["input_param.ILU_type"] = "1";
    description["input_param.ILU_lfil"] = "0";
    description["input_param.ILU_droptol"] = "1.000000e-01";
    description["input_param.ILU_relax"] = "9.000000e-01";
    description["input_param.ILU_permtol"] = "1.000000e-03";
    description["input_param.Schwarz_mmsize"] = "200";
    description["input_param.Schwarz_maxlvl"] = "2";
    description["input_param.Schwarz_type"] = "1";
    description["input_param.AMG_type"] = "1";
    description["input_param.AMG_levels"] = "20";
    description["input_param.AMG_cycle_type"] = "1";
    description["input_param.AMG_smoother"] = "2";
    description["input_param.AMG_relaxation"] = "1.100000e+00";
    description["input_param.AMG_polynomial_degree"] = "3";
    description["input_param.AMG_presmooth_iter"] = "2";
    description["input_param.AMG_postsmooth_iter"] = "2";
    description["input_param.AMG_coarse_dof"] = "100";
    description["input_param.AMG_tol"] = "1.000000e-08";
    description["input_param.AMG_maxit"] = "1";
    description["input_param.AMG_ILU_levels"] = "0";
    description["input_param.AMG_coarse_scaling"] = "0";
    description["input_param.AMG_amli_degree"] = "2";
    description["input_param.AMG_nl_amli_krylov_type"] = "6";
    description["input_param.AMG_schwarz_levels"] = "0";
    description["input_param.AMG_coarsening_type"] = "1";
    description["input_param.AMG_interpolation_type"] = "1";
    description["input_param.AMG_strong_threshold"] = "3.000000e-01";
    description["input_param.AMG_truncation_threshold"] = "4.000000e-01";
    description["input_param.AMG_max_row_sum"] = "9.000000e-01";
    description["input_param.AMG_strong_coupled"] = "8.000000e-02";
    description["input_param.AMG_max_aggregation"] = "20";
    description["input_param.AMG_tentative_smooth"] = "6.700000e-01";
    description["input_param.AMG_smooth_filter"] = "0";
    description["itsolver_param.itsolver_type"] = "1";
    description["itsolver_param.precond_type"] = "2";
    description["itsolver_param.stop_type"] = "1";
    description["itsolver_param.maxit"] = "1000";
    description["itsolver_param.tol"] = "1.000000e-10";
    description["itsolver_param.restart"] = "20";
    description["itsolver_param.print_level"] = "3";
    description["AMG_param.AMG_type"] = "1";
    description["AMG_param.print_level"] = "3";
    description["AMG_param.maxit"] = "1";
    description["AMG_param.tol"] = "1.000000e-08";
    description["AMG_param.max_levels"] = "20";
    description["AMG_param.coarse_dof"] = "100";
    description["AMG_param.cycle_type"] = "1";
    description["AMG_param.smoother"] = "2";
    description["AMG_param.smooth_order"] = "1";
    description["AMG_param.presmooth_iter"] = "2";
    description["AMG_param.postsmooth_iter"] = "2";
    description["AMG_param.relaxation"] = "1.100000e+00";
    description["AMG_param.polynomial_degree"] = "3";
    description["AMG_param.coarse_scaling"] = "0";
    description["AMG_param.amli_degree"] = "2";
    description["AMG_param.amli_coef"] = "1.100000e+00";
    description["AMG_param.nl_amli_krylov_type"] = "6";
    description["AMG_param.coarsening_type"] = "1";
    description["AMG_param.interpolation_type"] = "1";
    description["AMG_param.strong_threshold"] = "3.000000e-01";
    description["AMG_param.max_row_sum"] = "9.000000e-01";
    description["AMG_param.truncation_threshold"] = "4.000000e-01";
    description["AMG_param.strong_coupled"] = "8.000000e-02";
    description["AMG_param.max_aggregation"] = "20";
    description["AMG_param.tentative_smooth"] = "6.700000e-01";
    description["AMG_param.smooth_filter"] = "0";
    description["AMG_param.ILU_levels"] = "0";
    description["AMG_param.ILU_type"] = "1";
    description["AMG_param.ILU_lfil"] = "0";
    description["AMG_param.ILU_droptol"] = "1.000000e-01";
    description["AMG_param.ILU_relax"] = "9.000000e-01";
    description["AMG_param.ILU_permtol"] = "1.000000e-03";
    description["AMG_param.schwarz_levels"] = "0";
    description["AMG_param.schwarz_mmsize"] = "200";
    description["AMG_param.schwarz_maxlvl"] = "2";
    description["AMG_param.schwarz_type"] = "1";
    description["ILU_param.print_level"] = "3";
    description["ILU_param.ILU_type"] = "1";
    description["ILU_param.ILU_lfil"] = "0";
    description["ILU_param.ILU_droptol"] = "1.000000e-01";
    description["ILU_param.ILU_relax"] = "9.000000e-01";
    description["ILU_param.ILU_permtol"] = "1.000000e-03";
    description["Schwarz_param.print_level"] = "3";
    description["Schwarz_param.schwarz_type"] = "1";
    description["Schwarz_param.schwarz_maxlvl"] = "2";
    description["Schwarz_param.schwarz_mmsize"] = "200";
    return description;
  } // Dune::ParameterTree defaultSettings()

  /**
   *  \attention  There is a const_cast inside, in order to forward non-const pointers to fasp. I hope they do not
   *              touch the matrix, but who knows...
   */
  virtual size_t apply(const MatrixType& _systemMatrix,
                       const VectorType& _rhsVector,
                       VectorType& solutionVector,
                       const Dune::ParameterTree description = defaultSettings()) const
  {
    const size_t maxIter = description.get<size_t>("maxIter");
    const ScalarType precision = description.get<ScalarType>("precision");
    // init system matrix and right hand side
    MatrixType& systemMatrix = const_cast<MatrixType&>(_systemMatrix);
    VectorType& rhsVector = const_cast<VectorType&>(_rhsVector);
    dCSRmat A;
    A.row = systemMatrix.rows();
    A.col = systemMatrix.cols();
    A.nnz = systemMatrix.backend().nonZeros();
    A.IA = systemMatrix.backend().outerIndexPtr();
    A.JA = systemMatrix.backend().innerIndexPtr();
    A.val = systemMatrix.backend().valuePtr();
    dvector f, x;
    f.row = rhsVector.size();
    f.val = rhsVector.backend().data();
    x.row = rhsVector.backend().rows();
    x.val = solutionVector.backend().data();
    // init parameters
    input_param inparam = initInputParams(maxIter, precision, description);
    itsolver_param itparam = initItsolverParams(maxIter, precision, description);
    AMG_param amgparam = initAMGParams(1, precision, description); // the 1 is on purpose!
    ILU_param iluparam = initIluParams(maxIter, precision, description);
    Schwarz_param swzparam = initSchwarzParams(maxIter, precision, description);
    // call fasp (this is taken from the fasp example test)
    int status = -1;
    // Preconditioned Krylov methods
    if (inparam.solver_type >= 1 && inparam.solver_type <= 20) {
      // Using no preconditioner for Krylov iterative methods
      if (inparam.precond_type == PREC_NULL) {
        status = fasp_solver_dcsr_krylov(&A, &f, &x, &itparam);
      }
      // Using diag(A) as preconditioner for Krylov iterative methods
      else if (inparam.precond_type == PREC_DIAG) {
        status = fasp_solver_dcsr_krylov_diag(&A, &f, &x, &itparam);
      }
      // Using AMG as preconditioner for Krylov iterative methods
      else if (inparam.precond_type == PREC_AMG || inparam.precond_type == PREC_FMG) {
        if (inparam.print_level > PRINT_NONE)
          fasp_param_amg_print(&amgparam);
        status = fasp_solver_dcsr_krylov_amg(&A, &f, &x, &itparam, &amgparam);
      }
      // Using ILU as preconditioner for Krylov iterative methods Q: Need to change!
      else if (inparam.precond_type == PREC_ILU) {
        if (inparam.print_level > PRINT_NONE)
          fasp_param_ilu_print(&iluparam);
        status = fasp_solver_dcsr_krylov_ilu(&A, &f, &x, &itparam, &iluparam);
      }
      // Using Schwarz as preconditioner for Krylov iterative methods
      else if (inparam.precond_type == PREC_SCHWARZ) {
        if (inparam.print_level > PRINT_NONE)
          fasp_param_schwarz_print(&swzparam);
        status = fasp_solver_dcsr_krylov_schwarz(&A, &f, &x, &itparam, &swzparam);
      } else {
        printf("### ERROR: Wrong preconditioner type %d!!!\n", inparam.precond_type);
        status = ERROR_SOLVER_PRECTYPE;
      }
    }
    // AMG as the iterative solver
    else if (inparam.solver_type == SOLVER_AMG) {
      if (inparam.print_level > PRINT_NONE)
        fasp_param_amg_print(&amgparam);
      fasp_solver_amg(&A, &f, &x, &amgparam);
    }
    // Full AMG as the iterative solver
    else if (inparam.solver_type == SOLVER_FMG) {
      if (inparam.print_level > PRINT_NONE)
        fasp_param_amg_print(&amgparam);
      fasp_solver_famg(&A, &f, &x, &amgparam);
    } else {
      DUNE_THROW(Dune::RangeError, "### ERROR: Wrong solver type: " << inparam.solver_type << "!");
      status = ERROR_SOLVER_TYPE;
    }
    if (status > 0)
      return 0;
    else
      return 3;
  } // ... apply(...)

private:
  input_param
  initInputParams(const size_t& maxIter, const ScalarType& precision, const Dune::ParameterTree& description) const
  {
    input_param inputParam;
    inputParam.print_level = description.get<int>("input_param.print_level", 0);
    inputParam.output_type = description.get<int>("input_param.output_type", 0);
    //    inputParam.workdir = description.get< char >("input_param.workdir", '.');
    inputParam.problem_num = 14;
    inputParam.solver_type = description.get<int>("input_param.solver_type", 1);
    inputParam.precond_type = description.get<int>("input_param.precond_type", 2);
    inputParam.stop_type = description.get<int>("input_param.stop_type", 1);
    inputParam.itsolver_tol = precision;
    inputParam.itsolver_maxit = maxIter;
    inputParam.restart = description.get<int>("input_param.restart", 20);
    inputParam.ILU_type = description.get<int>("input_param.ILU_type", 1);
    inputParam.ILU_lfil = description.get<int>("input_param.ILU_lfil", 0);
    inputParam.ILU_droptol = description.get<double>("input_param.ILU_droptol", 1.000000e-01);
    inputParam.ILU_relax = description.get<double>("input_param.ILU_relax", 9.000000e-01);
    inputParam.ILU_permtol = description.get<double>("input_param.ILU_permtol", 1.000000e-03);
    inputParam.Schwarz_mmsize = description.get<int>("input_param.Schwarz_mmsize", 200);
    inputParam.Schwarz_maxlvl = description.get<int>("input_param.Schwarz_maxlvl", 2);
    inputParam.Schwarz_type = description.get<int>("input_param.Schwarz_type", 1);
    inputParam.AMG_type = description.get<int>("input_param.AMG_type", 1);
    inputParam.AMG_levels = description.get<int>("input_param.AMG_levels", 20);
    inputParam.AMG_cycle_type = description.get<int>("input_param.AMG_cycle_type", 1);
    inputParam.AMG_smoother = description.get<int>("input_param.AMG_smoother", 2);
    inputParam.AMG_relaxation = description.get<double>("input_param.AMG_relaxation", 1.100000e+00);
    inputParam.AMG_polynomial_degree = description.get<int>("input_param.AMG_polynomial_degree", 3);
    inputParam.AMG_presmooth_iter = description.get<int>("input_param.AMG_presmooth_iter", 2);
    inputParam.AMG_postsmooth_iter = description.get<int>("input_param.AMG_postsmooth_iter", 2);
    inputParam.AMG_coarse_dof = description.get<int>("input_param.AMG_coarse_dof", 100);
    inputParam.AMG_tol = description.get<double>("input_param.AMG_tol", 1.000000e-08);
    inputParam.AMG_maxit = description.get<int>("input_param.AMG_maxit", 1);
    inputParam.AMG_ILU_levels = description.get<int>("input_param.AMG_ILU_levels", 0);
    inputParam.AMG_coarse_scaling = description.get<int>("input_param.AMG_coarse_scaling", 0);
    inputParam.AMG_amli_degree = description.get<int>("input_param.AMG_amli_degree", 2);
    inputParam.AMG_nl_amli_krylov_type = description.get<int>("input_param.AMG_nl_amli_krylov_type", 6);
    inputParam.AMG_schwarz_levels = description.get<int>("input_param.AMG_schwarz_levels", 0);
    inputParam.AMG_coarsening_type = description.get<int>("input_param.AMG_coarsening_type", 1);
    inputParam.AMG_interpolation_type = description.get<int>("input_param.AMG_interpolation_type", 1);
    inputParam.AMG_strong_threshold = description.get<double>("input_param.AMG_strong_threshold", 3.000000e-01);
    inputParam.AMG_truncation_threshold = description.get<double>("input_param.AMG_truncation_threshold", 4.000000e-01);
    inputParam.AMG_max_row_sum = description.get<double>("input_param.AMG_max_row_sum", 9.000000e-01);
    inputParam.AMG_strong_coupled = description.get<double>("input_param.AMG_strong_coupled", 8.000000e-02);
    inputParam.AMG_max_aggregation = description.get<int>("input_param.AMG_max_aggregation", 20);
    inputParam.AMG_tentative_smooth = description.get<double>("input_param.AMG_tentative_smooth", 6.700000e-01);
    inputParam.AMG_smooth_filter = description.get<int>("input_param.AMG_smooth_filter", 0);
    return inputParam;
  } // ... initInputParams(...)

  itsolver_param
  initItsolverParams(const size_t& maxIter, const ScalarType& precision, const Dune::ParameterTree& description) const
  {
    itsolver_param itsolverParams;
    itsolverParams.itsolver_type = description.get<int>("itsolver_param.itsolver_type", 1);
    itsolverParams.precond_type = description.get<int>("itsolver_param.precond_type", 2);
    itsolverParams.stop_type = description.get<int>("itsolver_param.stop_type", 1);
    itsolverParams.maxit = maxIter;
    itsolverParams.tol = precision;
    itsolverParams.restart = description.get<int>("itsolver_param.restart", 20);
    itsolverParams.print_level = description.get<int>("itsolver_param.print_level", 0);
    return itsolverParams;
  } // ... initItsolverParams(...)

  AMG_param
  initAMGParams(const size_t& maxIter, const ScalarType& precision, const Dune::ParameterTree& description) const
  {
    AMG_param amgParams;
    amgParams.AMG_type = description.get<int>("AMG_param.AMG_type", 1);
    amgParams.print_level = description.get<int>("AMG_param.print_level", 0);
    amgParams.maxit = maxIter;
    amgParams.tol = precision;
    amgParams.max_levels = description.get<int>("AMG_param.max_levels", 20);
    amgParams.coarse_dof = description.get<int>("AMG_param.coarse_dof", 100);
    amgParams.cycle_type = description.get<int>("AMG_param.cycle_type", 1);
    amgParams.smoother = description.get<int>("AMG_param.smoother", 2);
    amgParams.smooth_order = description.get<int>("AMG_param.smooth_order", 1);
    amgParams.presmooth_iter = description.get<int>("AMG_param.presmooth_iter", 2);
    amgParams.postsmooth_iter = description.get<int>("AMG_param.postsmooth_iter", 2);
    amgParams.relaxation = description.get<double>("AMG_param.relaxation", 1.1);
    amgParams.polynomial_degree = description.get<int>("AMG_param.polynomial_degree", 3);
    amgParams.coarse_scaling = description.get<int>("AMG_param.coarse_scaling", 0);
    amgParams.amli_degree = description.get<int>("AMG_param.amli_degree", 2);
    double tmp = description.get<double>("AMG_param.amli_coef", 1.1);
    amgParams.amli_coef = &tmp;
    amgParams.nl_amli_krylov_type = description.get<int>("AMG_param.nl_amli_krylov_type", 6);
    amgParams.coarsening_type = description.get<int>("AMG_param.coarsening_type", 1);
    amgParams.interpolation_type = description.get<int>("AMG_param.interpolation_type", 1);
    amgParams.strong_threshold = description.get<double>("AMG_param.strong_threshold", 3.e-1);
    amgParams.max_row_sum = description.get<double>("AMG_param.max_row_sum", 9.0e-1);
    amgParams.truncation_threshold = description.get<double>("AMG_param.truncation_threshold", 1.0e-1);
    amgParams.strong_coupled = description.get<double>("AMG_param.strong_coupled", 8.0e-2);
    amgParams.max_aggregation = description.get<int>("AMG_param.max_aggregation", 20);
    amgParams.tentative_smooth = description.get<double>("AMG_param.tentative_smooth", 6.7e-1);
    amgParams.smooth_filter = description.get<int>("AMG_param.smooth_filter", 0);
    amgParams.ILU_levels = description.get<int>("AMG_param.ILU_levels", 0);
    amgParams.ILU_type = description.get<int>("AMG_param.ILU_type", 1);
    amgParams.ILU_lfil = description.get<int>("AMG_param.ILU_lfil", 0);
    amgParams.ILU_droptol = description.get<double>("AMG_param.ILU_droptol", 1.0e-1);
    amgParams.ILU_relax = description.get<double>("AMG_param.ILU_relax", 9.0e-1);
    amgParams.ILU_permtol = description.get<double>("AMG_param.ILU_permtol", 1.0e-3);
    amgParams.schwarz_levels = description.get<int>("AMG_param.schwarz_levels", 0);
    amgParams.schwarz_mmsize = description.get<int>("AMG_param.schwarz_mmsize", 200);
    amgParams.schwarz_maxlvl = description.get<int>("AMG_param.schwarz_maxlvl", 2);
    amgParams.schwarz_type = description.get<int>("AMG_param.schwarz_type", 1);
    return amgParams;
  } // ... initAMGParams(...)

  ILU_param initIluParams(const size_t& /*maxIter*/,
                          const ScalarType& /*precision*/,
                          const Dune::ParameterTree& description) const
  {
    ILU_param iluParams;
    iluParams.print_level = description.get<int>("ILU_param.print_level", 0);
    iluParams.ILU_type = description.get<int>("ILU_param.ILU_type", 1);
    iluParams.ILU_lfil = description.get<int>("ILU_param.ILU_lfil", 0);
    iluParams.ILU_droptol = description.get<double>("ILU_param.ILU_droptol", 1.0e-1);
    iluParams.ILU_relax = description.get<double>("ILU_param.ILU_relax", 9.0e-1);
    iluParams.ILU_permtol = description.get<double>("ILU_param.ILU_permtol", 1.0e-3);
    return iluParams;
  } // ... initIluParams(...)

  Schwarz_param initSchwarzParams(const size_t& /*maxIter*/,
                                  const ScalarType& /*precision*/,
                                  const Dune::ParameterTree& description) const
  {
    Schwarz_param schwarzParams;
    schwarzParams.print_level = description.get<int>("schwarzParams.print_level", 0);
    schwarzParams.schwarz_type = description.get<int>("schwarzParams.schwarz_type", 1);
    schwarzParams.schwarz_maxlvl = description.get<int>("schwarzParams.schwarz_maxlvl", 2);
    schwarzParams.schwarz_mmsize = description.get<int>("schwarzParams.schwarz_mmsize", 200);
    return schwarzParams;
  } // ... initSchwarzParams(...)
}; // class AmgSolver

} // namespace LA
} // namespace XT
} // namespace Dune

#  endif // HAVE_EIGEN
#endif // HAVE_FASP

#endif // DUNE_XT_LA_SOLVER_FASP_HH
