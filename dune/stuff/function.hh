#ifndef DUNE_STUFF_FUNCTION_HH
#define DUNE_STUFF_FUNCTION_HH

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#elif defined(HAVE_CONFIG_H)
#include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

#include <string>
#include <vector>

#include <dune/common/shared_ptr.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/exceptions.hh>

#include <dune/stuff/common/color.hh>

namespace Dune {
namespace Stuff {

// some forwards, includes are below
template <class D, int d, class R, int r>
class FunctionInterface;
template <class D, int d, class R, int r>
class FunctionCheckerboard;
template <class D, int d, class R, int r>
class FunctionExpression;
template <class D, int d, class R, int r>
class FunctionSpe10Model1;
template <class D, int d, class R, int r>
class FunctionAffineParametricDefault;
template <class D, int d, class R, int r>
class FunctionAffineParametricCheckerboard;


std::vector<std::string> availableFunctions()
{
  return {"function.expression",
          "function.checkerboard",
          "function.spe10.model1",
          "function.affineparametric.default",
          "function.affineparametric.checkerboard"};
} // std::vector< std::string > types()


template <class D, int d, class R, int r>
Dune::ParameterTree createSampleFunctionDescription(const std::string type)
{
  if (type == "function.checkerboard") {
    typedef FunctionCheckerboard<D, d, R, r> FunctionType;
    return FunctionType::createSampleDescription();
  } else if (type == "function.expression") {
    typedef FunctionExpression<D, d, R, r> FunctionType;
    return FunctionType::createSampleDescription();
  } else if (type == "function.spe10.model1") {
    typedef FunctionSpe10Model1<D, d, R, r> FunctionType;
    return FunctionType::createSampleDescription();
  } else if (type == "function.affineparametric.default") {
    typedef FunctionAffineParametricDefault<D, d, R, r> FunctionType;
    return FunctionType::createSampleDescription();
  } else if (type == "function.affineparametric.checkerboard") {
    typedef FunctionAffineParametricCheckerboard<D, d, R, r> FunctionType;
    return FunctionType::createSampleDescription();
  } else
    DUNE_THROW(Dune::RangeError,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " unknown function '" << type
                    << "' requested!");
} // ... createSampleFunctionDescription(...)


template <class D, int d, class R, int r>
FunctionInterface<D, d, R, r>* createFunction(const std::string type,
                                              const Dune::ParameterTree description = Dune::ParameterTree())
{
  if (type == "function.checkerboard") {
    return FunctionCheckerboard<D, d, R, r>::create(description);
  } else if (type == "function.expression") {
    return FunctionExpression<D, d, R, r>::create(description);
  } else if (type == "function.spe10.model1") {
    return FunctionSpe10Model1<D, d, R, r>::create(description);
  } else if (type == "function.affineparametric.default") {
    return FunctionAffineParametricDefault<D, d, R, r>::create(description);
  } else if (type == "function.affineparametric.checkerboard") {
    return FunctionAffineParametricCheckerboard<D, d, R, r>::create(description);
  } else
    DUNE_THROW(Dune::RangeError,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " unknown function '" << type
                    << "' requested!");
} // ... createFunction(...)


} // namespace Stuff
} // namespace Dune

#include "function/interface.hh"
#include "function/checkerboard.hh"
#include "function/expression.hh"
#include "function/spe10.hh"
#include "function/affineparametric/default.hh"
#include "function/affineparametric/checkerboard.hh"

#endif // DUNE_STUFF_FUNCTION_HH
