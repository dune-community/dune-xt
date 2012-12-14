#ifndef DUNE_STUFF_LA_SOLVER_EIGEN_HH
#define DUNE_STUFF_LA_SOLVER_EIGEN_HH

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#elif defined(HAVE_CONFIG_H)
#include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#if HAVE_EIGEN

#include <dune/stuff/common/color.hh>

#include "eigen/interface.hh"
#include "eigen/iterative.hh"
//#include "eigen/direct.hh"

namespace Dune {
namespace Stuff {
namespace LA {
namespace Solver {
namespace Eigen {


template <class MatrixType>
Interface<MatrixType>* create(const std::string type = "eigen.iterative.bicgstab.diagonal")
{
  if (type == "eigen.iterative.bicgstab.diagonal") {
    return new BicgstabDiagonal<MatrixType>;
  } else if (type == "eigen.iterative.bicgstab.incompletelut") {
    return new BicgstabIncompleteLUT<MatrixType>;
  } else if (type == "eigen.iterative.cg.diagonal") {
    std::cout << "\n" << Dune::Stuff::Common::highlightString("WARNING:", 9) << " you selected the solver '" << type
              << "' which is believed to produce strange results!" << std::endl;
    return new CgDiagonal<MatrixType>;
  } else if (type == "eigen.iterative.cg.incompletelut") {
    std::cout << "\n" << Dune::Stuff::Common::highlightString("WARNING:", 9) << " you selected the solver '" << type
              << "' which is believed to not converge!" << std::endl;
    return new CgIncompleteLUT<MatrixType>;
    //  } else if (type == "eigen.direct.simplicialllt") {
    //    return new SimplicialLLT< MatrixType >;
    //  } else if (type == "eigen.direct.simplicialldlt") {
    //    return new SimplicialLDLT< MatrixType >;
  } else
    DUNE_THROW(Dune::RangeError, "\nERROR: unknown linear solver '" << type << "' requested!");
} // Interface< ElementType >* create(const std::string type = "eigen.bicgstab.incompletelut")


} // namespace Eigen
} // namespace Solver
} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // HAVE_EIGEN

#endif // DUNE_STUFF_LA_SOLVER_EIGEN_HH
