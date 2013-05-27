#ifndef DUNE_STUFF_LA_SOLVER_HH
#define DUNE_STUFF_LA_SOLVER_HH

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#elif defined(HAVE_CONFIG_H)
#include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#include <string>
#include <vector>

#include <dune/stuff/common/color.hh>
#include <dune/stuff/common/parameter/tree.hh>

#include "solver/interface.hh"
#include "solver/eigen.hh"
#include "solver/fasp.hh"

namespace Dune {
namespace Stuff {
namespace LA {


static std::vector<std::string> solverTypes()
{
  return
  {
    "bicgstab.ilut", "bicgstab.diagonal", "bicgstab", "cg", "cg.diagonal", "simplicialllt", "simplicialldlt"
#if HAVE_FASP
        ,
        "amg.fasp"
#endif // FAVE_FASP
  };
} // std::vector< std::string > types()


template <class MatrixType, class VectorType>
Dune::ParameterTree solverDefaultSettings(const std::string type)
{
  if (type == "bicgstab.ilut") {
    return Dune::Stuff::LA::BicgstabILUTSolver<MatrixType, VectorType>::defaultSettings();
#if HAVE_FASP
  } else if (type == "amg.fasp") {
    return Dune::Stuff::LA::AmgFaspSolver<MatrixType, VectorType>::defaultSettings();
#endif // HAVE_FASP
  }
  return SolverInterface<MatrixType, VectorType>::defaultIterativeSettings();
} // Dune::ParameterTree defaultSettings(const std::string type)


template <class MatrixType, class VectorType>
SolverInterface<MatrixType, VectorType>* createSolver(const std::string type = solverTypes()[0])
{
  if (type == "cg") {
    return new Dune::Stuff::LA::CgSolver<MatrixType, VectorType>();
  } else if (type == "cg.diagonal") {
    return new Dune::Stuff::LA::CgDiagonalSolver<MatrixType, VectorType>();
  } else if (type == "bicgstab") {
    return new Dune::Stuff::LA::BicgstabSolver<MatrixType, VectorType>();
  } else if (type == "bicgstab.diagonal") {
    return new Dune::Stuff::LA::BicgstabDiagonalSolver<MatrixType, VectorType>();
  } else if (type == "bicgstab.ilut") {
    return new Dune::Stuff::LA::BicgstabILUTSolver<MatrixType, VectorType>();
  } else if (type == "simplicialllt") {
    return new Dune::Stuff::LA::SimplicialLLTSolver<MatrixType, VectorType>();
  } else if (type == "simplicialldlt") {
    return new Dune::Stuff::LA::SimplicialLDLTSolver<MatrixType, VectorType>();
#if HAVE_FASP
  } else if (type == "amg.fasp") {
    return new Dune::Stuff::LA::AmgFaspSolver<MatrixType, VectorType>();
#endif // HAVE_FASP
  } else
    DUNE_THROW(Dune::RangeError,
               "\n" << Dune::Stuff::Common::colorString("ERROR:", Dune::Stuff::Common::Colors::red)
                    << " unknown linear solver '"
                    << type
                    << "' requested!");
} // SolverInterface< ElementType >* create(const std::string type = "eigen.bicgstab.incompletelut")


} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_LA_SOLVER_HH
