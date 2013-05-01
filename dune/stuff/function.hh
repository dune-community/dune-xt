#ifndef DUNE_STUFF_FUNCTION_HH
#define DUNE_STUFF_FUNCTION_HH

#include <string>
#include <vector>

#include <dune/common/parametertree.hh>
#include <dune/common/exceptions.hh>

#include <dune/stuff/common/color.hh>

namespace Dune {
namespace Stuff {


// function forwards, includes are below
template <class D, int d, class R, int rR, int rC>
class FunctionInterface;
template <class D, int d, class R, int rR, int rC>
class FunctionCheckerboard;
template <class D, int d, class R, int rR, int rC>
class FunctionExpression;
template <class D, int d, class R, int rR, int rC>
class FunctionConstant;
template <class D, int d, class R, int rR, int rC>
class FunctionSpe10Model1;
// affine parametric forwards, includes are below
template <class D, int d, class R, int rR, int rC>
class AffineParametricFunctionCheckerboard;


/**
 *  \brief      This is a way to create new functions fulfilling the FunctionInterface using their string identifier.
 *
 *  \attention  This class will not compile for all dimensions. The errors should give you a hint which specializations
 *              are needed below.
 */
template <class D, int d, class R, int rR, int rC>
class Functions
{
public:
  static std::vector<std::string> available()
  {
    return {"function.expression", "function.checkerboard", "function.constant", "function.spe10.model1"};
  } // ... available(...)

  static Dune::ParameterTree createSampleDescription(const std::string type)
  {
    if (type == "function.checkerboard")
      return FunctionCheckerboard<D, d, R, rR, rC>::createSampleDescription();
    else if (type == "function.expression")
      return FunctionExpression<D, d, R, rR, rC>::createSampleDescription();
    else if (type == "function.constant")
      return FunctionConstant<D, d, R, rR, rC>::createSampleDescription();
    else if (type == "function.spe10.model1")
      return FunctionSpe10Model1<D, d, R, rR, rC>::createSampleDescription();
    else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " unknown function '" << type
                      << "' requested!");
  } // ... createSampleDescription(...)

  static FunctionInterface<D, d, R, rR, rC>* create(const std::string type,
                                                    const Dune::ParameterTree description = Dune::ParameterTree())
  {
    if (type == "function.checkerboard")
      return FunctionCheckerboard<D, d, R, rR, rC>::create(description);
    else if (type == "function.expression")
      return FunctionExpression<D, d, R, rR, rC>::create(description);
    else if (type == "function.constant")
      return FunctionConstant<D, d, R, rR, rC>::create(description);
    else if (type == "function.spe10.model1")
      return FunctionSpe10Model1<D, d, R, rR, rC>::createSampleDescription();
    else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " unknown function '" << type
                      << "' requested!");
  } // ... create(...)
}; // class Functions


/**
 *  \brief      This is a way to create new affine parametric functions fulfilling the AffineParametricFunctionInterface
 *              using their string identifier.
 *
 *  \attention  This class will not compile for all dimensions. The errors should give you a hint which specializations
 *              are needed below.
 */
template <class D, int d, class R, int rR, int rC>
class AffineParametricFunctions
{
public:
  static std::vector<std::string> available()
  {
    return {"function.affineparametric.checkerboard"};
  } // ... available(...)

  static Dune::ParameterTree createSampleDescription(const std::string type)
  {
    if (type == "function.affineparametric.checkerboard")
      return AffineParametricFunctionCheckerboard<D, d, R, rR, rC>::createSampleDescription();
    else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " unknown function '" << type
                      << "' requested!");
  } // ... createSampleDescription(...)

  static FunctionInterface<D, d, R, rR, rC>* create(const std::string type,
                                                    const Dune::ParameterTree description = Dune::ParameterTree())
  {
    if (type == "function.affineparametric.checkerboard")
      return AffineParametricFunctionCheckerboard<D, d, R, rR, rC>::create(description);
    else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " unknown function '" << type
                      << "' requested!");
  } // ... create(...)
}; // class AffineParametricFunctions


} // namespace Stuff
} // namespace Dune

#include "function/interface.hh"
#include "function/checkerboard.hh"
#include "function/expression.hh"
#include "function/spe10.hh"
#include "function/constant.hh"
#include "function/affineparametric/checkerboard.hh"
//#include "function/affineparametric/default.hh"

#endif // DUNE_STUFF_FUNCTION_HH
