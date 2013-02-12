#ifndef DUNE_STUFF_LA_SOLVER_HH
#define DUNE_STUFF_LA_SOLVER_HH

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#elif defined(HAVE_CONFIG_H)
#include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#include <dune/common/shared_ptr.hh>
#include <dune/stuff/common/color.hh>
#include <dune/stuff/common/parameter/tree.hh>

#include "solver/interface.hh"
#include "solver/eigen.hh"

namespace Dune {
namespace Stuff {
namespace LA {
namespace Solver {


Dune::ParameterTree createSampleDescription()
{
  Dune::ParameterTree description;
  description["maxIter"]   = "5000";
  description["precision"] = "1e-12";
  return description;
} // Dune::ParameterTree createSampleDescription(const std::string type)


template <class MatrixType, class VectorType>
Dune::shared_ptr<Interface<MatrixType, VectorType>> create(const std::string type = "bicgstab")
{
  if (type == "cg") {
    return Dune::make_shared<Dune::Stuff::LA::Solver::Cg<MatrixType, VectorType>>();
  } else if (type == "cg.diagonal") {
    return Dune::make_shared<Dune::Stuff::LA::Solver::CgDiagonal<MatrixType, VectorType>>();
  } else if (type == "bicgstab") {
    return Dune::make_shared<Dune::Stuff::LA::Solver::Bicgstab<MatrixType, VectorType>>();
  } else if (type == "bicgstab.diagonal") {
    return Dune::make_shared<Dune::Stuff::LA::Solver::BicgstabDiagonal<MatrixType, VectorType>>();
  } else if (type == "bicgstab.ilut") {
    return Dune::make_shared<Dune::Stuff::LA::Solver::BicgstabILUT<MatrixType, VectorType>>();
  } else if (type == "simplicialllt") {
    return Dune::make_shared<Dune::Stuff::LA::Solver::SimplicialLLT<MatrixType, VectorType>>();
  } else if (type == "simplicialldlt") {
    return Dune::make_shared<Dune::Stuff::LA::Solver::SimplicialLDLT<MatrixType, VectorType>>();
  } else
    DUNE_THROW(Dune::RangeError,
               "\n" << Dune::Stuff::Common::colorString("ERROR:", Dune::Stuff::Common::Colors::red)
                    << " unknown linear solver '"
                    << type
                    << "' requested!");
} // Interface< ElementType >* create(const std::string type = "eigen.bicgstab.incompletelut")


} // namespace Solver
} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_LA_SOLVER_HH
