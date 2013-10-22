// This file is part of the dune-stuff project:
//   http://users.dune-project.org/projects/dune-stuff/
// Copyright Holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_FUNCTIONS_HH
#define DUNE_STUFF_FUNCTIONS_HH

#include <string>
#include <vector>

#include <dune/common/parametertree.hh>
#include <dune/common/exceptions.hh>

#include <dune/stuff/common/color.hh>

#include "functions/interfaces.hh"
#include "functions/checkerboard.hh"
#include "functions/expression.hh"
#include "functions/constant.hh"
#include "functions/spe10.hh"

namespace Dune {
namespace Stuff {


/**
 *  \brief  This is a way to create new functions fulfilling the LocalizableFunctionInterface using their string
 *          identifier.
 *
 *  \note If you implement a new function, you should add it to the libdunestuff by editing the functions.cc file as
 *        follows. Inside the DUNE_STUFF_FUNCTIONS_CC_GENERATE_SPECIALIZATION macro, for each function, you will find
 *        three lines like
 *          FUNCTION_Expression_ ## ddim ## _ ## rdim ## _ ## rCdim ## _STATIC_ID(E, D, ddim, R, rdim, rCdim) \
 *        Yust copy each of those and replace Expression by the class name of your Function. In addition, below the
 *        DUNE_STUFF_FUNCTIONS_CC_GENERATE_SPECIALIZATION macro you will find a very long list of macro definitions like
 *          FUNCTION_Expression_1_1_1_ ... (E, D, ddim, R, rdim, rCdim)
 *        for each function with ... replaced by STATIC_ID, SETTINGS and CREATE. You will have to copy those 81 lines of
 *        code for your function and insert which combination of dimensions you did implement and which not. This might
 *        seem totally ridiculous but it can be acchieved quite fast by a few find and replace commands of your favorite
 *        text editor. Patches welcome if you have a better idea!
 */
template <class E, class D, int d, class R, int r, int rC = 1>
class Functions
{
public:
  static std::vector<std::string> available() = delete;

  static Dune::ParameterTree defaultSettings(const std::string /*type = available()[0]*/) = delete;

  static LocalizableFunctionInterface<E, D, d, R, r, rC>*
  create(const std::string /*type = available()[0]*/, const ParameterTree /*settings = defaultSettings()*/) = delete;
}; // class Functions

#define DUNE_STUFF_FUNCTIONS_GENERATE_DIM_RANGE(ddim)                                                                  \
  DUNE_STUFF_FUNCTIONS_GENERATE_DIM_RANGE_COLS(ddim, 1)                                                                \
  DUNE_STUFF_FUNCTIONS_GENERATE_DIM_RANGE_COLS(ddim, 2)                                                                \
  DUNE_STUFF_FUNCTIONS_GENERATE_DIM_RANGE_COLS(ddim, 3)

#define DUNE_STUFF_FUNCTIONS_GENERATE_DIM_RANGE_COLS(ddim, rdim)                                                       \
  DUNE_STUFF_FUNCTIONS_GENERATE_CLASS(ddim, rdim, 1)                                                                   \
  DUNE_STUFF_FUNCTIONS_GENERATE_CLASS(ddim, rdim, 2)                                                                   \
  DUNE_STUFF_FUNCTIONS_GENERATE_CLASS(ddim, rdim, 3)

#define DUNE_STUFF_FUNCTIONS_GENERATE_CLASS(ddim, rdim, rCdim)                                                         \
  template <class E, class D, class R>                                                                                 \
  class Functions<E, D, ddim, R, rdim, rCdim>                                                                          \
  {                                                                                                                    \
    static const int d  = ddim;                                                                                        \
    static const int r  = rdim;                                                                                        \
    static const int rC = rCdim;                                                                                       \
                                                                                                                       \
  public:                                                                                                              \
    static std::vector<std::string> available();                                                                       \
                                                                                                                       \
    static Dune::ParameterTree defaultSettings(const std::string type = available()[0]);                               \
                                                                                                                       \
    static LocalizableFunctionInterface<E, D, ddim, R, rdim, rCdim>*                                                   \
    create(const std::string type = available()[0], const ParameterTree settings = defaultSettings());                 \
  }; // class Functions< ..., ddim, ..., rdim, rCdim >

DUNE_STUFF_FUNCTIONS_GENERATE_DIM_RANGE(1)
DUNE_STUFF_FUNCTIONS_GENERATE_DIM_RANGE(2)
DUNE_STUFF_FUNCTIONS_GENERATE_DIM_RANGE(3)

#undef DUNE_STUFF_FUNCTIONS_GENERATE_CLASS
#undef DUNE_STUFF_FUNCTIONS_GENERATE_DIM_RANGE_COLS
#undef DUNE_STUFF_FUNCTIONS_GENERATE_DIM_RANGE

} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTIONS_HH
