#ifndef DUNE_STUFF_FUNCTION_HH
#define DUNE_STUFF_FUNCTION_HH

#include <string>
#include <vector>

#include <dune/common/parametertree.hh>
#include <dune/common/exceptions.hh>

#include <dune/stuff/common/color.hh>

namespace Dune {
namespace Stuff {


// forwards, includes are below
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


/**
 *  \brief      This is a way to create new functions fulfilling the FunctionInterface using their string identifier.
 *
 *  \attention  This class will not compile for all dimensions. The errors should give you a hint which specializations
 *              are needed below.
 */
template <class D, int d, class R, int rR, int rC = 1>
class Functions
{
public:
  static std::vector<std::string> available()
  {
    return {"function.expression", "function.checkerboard", "function.constant", "function.spe10.model1"};
  } // ... available(...)

  static Dune::ParameterTree defaultSettings(const std::string type = available()[0])
  {
    if (type == "function.checkerboard")
      return FunctionCheckerboard<D, d, R, rR, rC>::defaultSettings();
    else if (type == "function.expression")
      return FunctionExpression<D, d, R, rR, rC>::defaultSettings();
    else if (type == "function.constant")
      return FunctionConstant<D, d, R, rR, rC>::defaultSettings();
    else if (type == "function.spe10.model1")
      return FunctionSpe10Model1<D, d, R, rR, rC>::defaultSettings();
    else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " unknown function '" << type
                      << "' requested!");
  } // ... defaultSettings(...)

  static FunctionInterface<D, d, R, rR, rC>* create(const std::string type = available()[0],
                                                    const Dune::ParameterTree settings = defaultSettings())
  {
    if (type == "function.checkerboard")
      return FunctionCheckerboard<D, d, R, rR, rC>::create(settings);
    else if (type == "function.expression")
      return FunctionExpression<D, d, R, rR, rC>::create(settings);
    else if (type == "function.constant")
      return FunctionConstant<D, d, R, rR, rC>::create(settings);
    else if (type == "function.spe10.model1")
      return FunctionSpe10Model1<D, d, R, rR, rC>::create(settings);
    else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " unknown function '" << type
                      << "' requested!");
  } // ... create(...)
}; // class Functions


///**
// *  \brief      This is a way to create new affine parametric functions fulfilling the
// AffineParametricFunctionInterface
// *              using their string identifier.
// *
// *  \attention  This class will not compile for all dimensions. The errors should give you a hint which
// specializations
// *              are needed below.
// */
// template< class D, int d, class R, int rR, int rC = 1>
// class AffineParametricFunctions
//{
// public:
//  static std::vector< std::string > available()
//  {
//    return {
//          "function.affineparametric.checkerboard"
//        , "function.affineparametric.default"
//    };
//  } // ... available(...)

//  static Dune::ParameterTree defaultSettings(const std::string type = available()[0])
//  {
//    if (type == "function.affineparametric.checkerboard")
//      return AffineParametricFunctionCheckerboard< D, d, R, rR, rC >::defaultSettings();
//    else if (type == "function.affineparametric.default")
//      return AffineParametricFunctionDefault< D, d, R, rR, rC >::defaultSettings();
//    else
//      DUNE_THROW(Dune::RangeError,
//                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
//                 << " unknown function '" << type << "' requested!");
//  } // ... defaultSettings(...)

//  static AffineParametricFunctionInterface< D, d, R, rR, rC >* create(const std::string type = available()[0],
//                                                                      const Dune::ParameterTree settings =
//                                                                      defaultSettings())
//  {
//    if (type == "function.affineparametric.checkerboard")
//      return AffineParametricFunctionCheckerboard< D, d, R, rR, rC >::create(settings);
//    else if (type == "function.affineparametric.default")
//      return AffineParametricFunctionDefault< D, d, R, rR, rC >::create(settings);
//    else
//      DUNE_THROW(Dune::RangeError,
//                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
//                 << " unknown function '" << type << "' requested!");
//  } // ... create(...)
//}; // class AffineParametricFunctions


// template< class D, int d, class R, int rR, int rC = 1 >
// class GenericStationaryFunctions
//{
//  typedef Functions< D, d, R, rR, rC > Funcs;
//  typedef AffineParametricFunctions< D, d, R, rR, rC > AffineFuncs;

// public:
//  static std::vector< std::string > available()
//  {
//    return {Funcs::available(), AffineFuncs::available() };
//  }

//  static Dune::ParameterTree defaultSettings(const std::string type = available()[0])
//  {
//    const auto funcs = Funcs::available();
//    const auto affineFuncs = AffineFuncs::available();
//    bool found = true;
//    for (auto& func : funcs)
//      if (type == func)
//        found = true;
//    if (found)
//      return Funcs::defaultSettings(type);
//    for (auto& affineFunc : affineFuncs)
//      if (type == affineFunc)
//        found = true;
//    if (found)
//      return AffineFuncs::defaultSettings(type);
//    else
//      DUNE_THROW(Dune::RangeError,
//                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
//                 << " unknown function '" << type << "' requested!");
//  }

//  static GenericStationaryFunctionInterface< D, d, R, rR, rC >* create(const std::string type = available()[0],
//                                                                       const Dune::ParameterTree settings =
//                                                                       defaultSettings())
//  {
//    const auto funcs = Funcs::available();
//    const auto affineFuncs = AffineFuncs::available();
//    bool found = false;
//    for (auto& func : funcs)
//      if (type == func)
//        found = true;
//    if (found)
//      return Funcs::create(type, settings);
//    for (auto& affineFunc : affineFuncs)
//      if (type == affineFunc)
//        found = true;
//    if (found)
//      return AffineFuncs::create(type, settings);
//    else
//      DUNE_THROW(Dune::RangeError,
//                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
//                 << " unknown function '" << type << "' requested!");
//  }
//}; // class GenericStationaryFunctions


} // namespace Stuff
} // namespace Dune

#include "function/interface.hh"
#include "function/checkerboard.hh"
#include "function/expression.hh"
#include "function/spe10.hh"
#include "function/constant.hh"
//#include "function/affineparametric/checkerboard.hh"
//#include "function/affineparametric/default.hh"

#endif // DUNE_STUFF_FUNCTION_HH
