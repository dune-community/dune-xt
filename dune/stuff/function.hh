#ifndef DUNE_STUFF_FUNCTION_HH
#define DUNE_STUFF_FUNCTION_HH

#include <string>
#include <vector>

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


template <class D, int d, class R, int r>
class Functions
{
public:
  static std::vector<std::string> available()
  {
    return {"function.expression",
            "function.checkerboard",
            "function.affineparametric.default",
            "function.affineparametric.checkerboard"};
  } // ... available(...)

  static Dune::ParameterTree createSampleDescription(const std::string type)
  {
    if (type == "function.checkerboard")
      return FunctionCheckerboard<D, d, R, r>::createSampleDescription();
    else if (type == "function.expression")
      return FunctionExpression<D, d, R, r>::createSampleDescription();
    else if (type == "function.affineparametric.default")
      return FunctionAffineParametricDefault<D, d, R, r>::createSampleDescription();
    else if (type == "function.affineparametric.checkerboard")
      return FunctionAffineParametricCheckerboard<D, d, R, r>::createSampleDescription();
    else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " unknown function '" << type
                      << "' requested!");
  } // ... createSampleDescription(...)

  static FunctionInterface<D, d, R, r>* create(const std::string type,
                                               const Dune::ParameterTree description = Dune::ParameterTree())
  {
    if (type == "function.checkerboard")
      return FunctionCheckerboard<D, d, R, r>::create(description);
    else if (type == "function.expression")
      return FunctionExpression<D, d, R, r>::create(description);
    else if (type == "function.affineparametric.default")
      return FunctionAffineParametricDefault<D, d, R, r>::create(description);
    else if (type == "function.affineparametric.checkerboard")
      return FunctionAffineParametricCheckerboard<D, d, R, r>::create(description);
    else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " unknown function '" << type
                      << "' requested!");
  } // ... create(...)
}; // class Functions


// spezialization for dims 2 -> 1
template <class D, class R>
class Functions<D, 2, R, 1>
{
public:
  static std::vector<std::string> available()
  {
    return {"function.expression",
            "function.checkerboard",
            "function.spe10.model1",
            "function.affineparametric.default",
            "function.affineparametric.checkerboard"};
  } // ... available(...)

  static Dune::ParameterTree createSampleDescription(const std::string type)
  {
    if (type == "function.checkerboard")
      return FunctionCheckerboard<D, 2, R, 1>::createSampleDescription();
    else if (type == "function.expression")
      return FunctionExpression<D, 2, R, 1>::createSampleDescription();
    else if (type == "function.spe10.model1")
      return FunctionSpe10Model1<D, 2, R, 1>::createSampleDescription();
    else if (type == "function.affineparametric.default")
      return FunctionAffineParametricDefault<D, 2, R, 1>::createSampleDescription();
    else if (type == "function.affineparametric.checkerboard")
      return FunctionAffineParametricCheckerboard<D, 2, R, 1>::createSampleDescription();
    else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " unknown function '" << type
                      << "' requested!");
  } // ... createSampleDescription(...)

  static FunctionInterface<D, 2, R, 1>* create(const std::string type,
                                               const Dune::ParameterTree description = Dune::ParameterTree())
  {
    if (type == "function.checkerboard")
      return FunctionCheckerboard<D, 2, R, 1>::create(description);
    else if (type == "function.expression")
      return FunctionExpression<D, 2, R, 1>::create(description);
    else if (type == "function.spe10.model1")
      return FunctionSpe10Model1<D, 2, R, 1>::create(description);
    else if (type == "function.affineparametric.default")
      return FunctionAffineParametricDefault<D, 2, R, 1>::create(description);
    else if (type == "function.affineparametric.checkerboard")
      return FunctionAffineParametricCheckerboard<D, 2, R, 1>::create(description);
    else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " unknown function '" << type
                      << "' requested!");
  } // ... create(...)
}; // class Functions< D, 2, R, 1 >


} // namespace Stuff
} // namespace Dune

#include "function/interface.hh"
#include "function/checkerboard.hh"
#include "function/expression.hh"
#include "function/spe10.hh"
#include "function/affineparametric/default.hh"
#include "function/affineparametric/checkerboard.hh"

#endif // DUNE_STUFF_FUNCTION_HH
