// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
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
 *  \note If you implement a new function, you should add it to the libdunestuff. Inside the
 *        DUNE_STUFF_FUNCTIONS_GENERATE_CLASS macro you will find three lines like (for each function) like
 *          DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_ ## ddim ## _ ## rdim ## _ ## rCdim ## _STATIC_ID(E, D, ddim, R,
 * rdim, rCdim) \
 *        Yust copy each of those and replace Expression by the class name of your Function. In addition, below the
 *        DUNE_STUFF_FUNCTIONS_GENERATE_CLASS macro you will find a very long list of macro definitions like
 *          DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_1_1_1_ ... (E, D, ddim, R, rdim, rCdim)
 *        for each function with ... replaced by STATIC_ID, SETTINGS and CREATE. You will have to copy those 81 lines of
 *        code for your function and insert which combination of dimensions you did implement and which not. This might
 *        seem totally ridiculous but it can be acchieved quite fast by a few find and replace commands of your favorite
 *        text editor. Patches welcome if you have a better idea!
 */
template <class E, class D, int d, class R, int r, int rC = 1>
class FunctionsProvider
{
public:
  static std::vector<std::string> available() = delete;

  static Dune::ParameterTree defaultSettings(const std::string /*type = available()[0]*/) = delete;

  static LocalizableFunctionInterface<E, D, d, R, r, rC>*
  create(const std::string /*type = available()[0]*/, const ParameterTree /*settings = defaultSettings()*/) = delete;
}; // class FunctionsProvider

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
  class FunctionsProvider<E, D, ddim, R, rdim, rCdim>                                                                  \
  {                                                                                                                    \
    static const int d  = ddim;                                                                                        \
    static const int r  = rdim;                                                                                        \
    static const int rC = rCdim;                                                                                       \
                                                                                                                       \
  public:                                                                                                              \
    static std::vector<std::string> available()                                                                        \
    {                                                                                                                  \
      return {Functions::Constant<E, D, d, R, r, rC>::static_id()                                                      \
                  DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_##ddim##_##rdim##_##rCdim##_STATIC_ID(                   \
                      E, D, ddim, R, rdim, rCdim)                                                                      \
                      DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_##ddim##_##rdim##_##rCdim##_STATIC_ID(             \
                          E, D, ddim, R, rdim, rCdim)                                                                  \
                          DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_##ddim##_##rdim##_##rCdim##_STATIC_ID(          \
                              E, D, ddim, R, rdim, rCdim)};                                                            \
    } /* ... available(...) */                                                                                         \
                                                                                                                       \
    static Dune::ParameterTree defaultSettings(const std::string type = available()[0])                                \
    {                                                                                                                  \
      if (type.compare(Functions::Constant<E, D, d, R, r, rC>::static_id()) == 0)                                      \
        return Functions::Constant<E, D, d, R, r, rC>::defaultSettings();                                              \
      DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_##ddim##_##rdim##_##rCdim##_SETTINGS(E, D, ddim, R, rdim, rCdim)     \
          DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_##ddim##_##rdim##_##rCdim##_SETTINGS(                          \
              E, D, ddim, R, rdim, rCdim)                                                                              \
              DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_##ddim##_##rdim##_##rCdim##_SETTINGS(                       \
                  E, D, ddim, R, rdim, rCdim) else DUNE_THROW(Dune::RangeError,                                        \
                                                              "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")    \
                                                                   << " unknown function '"                            \
                                                                   << type                                             \
                                                                   << "' requested!");                                 \
    } /* ... defaultSettings(...) */                                                                                   \
                                                                                                                       \
    static LocalizableFunctionInterface<E, D, ddim, R, rdim, rCdim>*                                                   \
    create(const std::string type = available()[0], const ParameterTree settings = defaultSettings())                  \
    {                                                                                                                  \
      if (type.compare(Functions::Constant<E, D, d, R, r, rC>::static_id()) == 0)                                      \
        return Functions::Constant<E, D, d, R, r, rC>::create(settings);                                               \
      DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_##ddim##_##rdim##_##rCdim##_CREATE(                                  \
          E, D, ddim, R, rdim, rCdim, settings)                                                                        \
          DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_##ddim##_##rdim##_##rCdim##_CREATE(                            \
              E, D, ddim, R, rdim, rCdim, settings)                                                                    \
              DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_##ddim##_##rdim##_##rCdim##_CREATE(                         \
                  E, D, ddim, R, rdim, rCdim, settings) else DUNE_THROW(Dune::RangeError,                              \
                                                                        "\n" << Dune::Stuff::Common::colorStringRed(   \
                                                                                    "ERROR:")                          \
                                                                             << " unknown function '"                  \
                                                                             << type                                   \
                                                                             << "' requested!");                       \
    } /* ... create(...) */                                                                                            \
  }; // class FunctionsProvider< ..., ddim, ..., rdim, rCdim >
// DUNE_STUFF_FUNCTIONS_GENERATE_CLASS

// ====================
// ==== Expression ====
// ====================
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_1_1_1_STATIC_ID(E, D, ddim, R, rdim, rCdim)                        \
  , Functions::Expression<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_1_1_2_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_1_1_3_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_1_2_1_STATIC_ID(E, D, ddim, R, rdim, rCdim)                        \
  , Functions::Expression<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_1_2_2_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_1_2_3_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_1_3_1_STATIC_ID(E, D, ddim, R, rdim, rCdim)                        \
  , Functions::Expression<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_1_3_2_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_1_3_3_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_2_1_1_STATIC_ID(E, D, ddim, R, rdim, rCdim)                        \
  , Functions::Expression<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_2_1_2_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_2_1_3_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_2_2_1_STATIC_ID(E, D, ddim, R, rdim, rCdim)                        \
  , Functions::Expression<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_2_2_2_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_2_2_3_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_2_3_1_STATIC_ID(E, D, ddim, R, rdim, rCdim)                        \
  , Functions::Expression<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_2_3_2_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_2_3_3_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_3_1_1_STATIC_ID(E, D, ddim, R, rdim, rCdim)                        \
  , Functions::Expression<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_3_1_2_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_3_1_3_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_3_2_1_STATIC_ID(E, D, ddim, R, rdim, rCdim)                        \
  , Functions::Expression<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_3_2_2_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_3_2_3_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_3_3_1_STATIC_ID(E, D, ddim, R, rdim, rCdim)                        \
  , Functions::Expression<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_3_3_2_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_3_3_3_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/

#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_1_1_1_SETTINGS(E, D, ddim, R, rdim, rCdim)                         \
  else if (type.compare(Functions::Expression<E, D, d, R, r, rC>::static_id())                                         \
           == 0) return Functions::Expression<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_1_1_2_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_1_1_3_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_1_2_1_SETTINGS(E, D, ddim, R, rdim, rCdim)                         \
  else if (type.compare(Functions::Expression<E, D, d, R, r, rC>::static_id())                                         \
           == 0) return Functions::Expression<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_1_2_2_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_1_2_3_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_1_3_1_SETTINGS(E, D, ddim, R, rdim, rCdim)                         \
  else if (type.compare(Functions::Expression<E, D, d, R, r, rC>::static_id())                                         \
           == 0) return Functions::Expression<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_1_3_2_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_1_3_3_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_2_1_1_SETTINGS(E, D, ddim, R, rdim, rCdim)                         \
  else if (type.compare(Functions::Expression<E, D, d, R, r, rC>::static_id())                                         \
           == 0) return Functions::Expression<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_2_1_2_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_2_1_3_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_2_2_1_SETTINGS(E, D, ddim, R, rdim, rCdim)                         \
  else if (type.compare(Functions::Expression<E, D, d, R, r, rC>::static_id())                                         \
           == 0) return Functions::Expression<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_2_2_2_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_2_2_3_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_2_3_1_SETTINGS(E, D, ddim, R, rdim, rCdim)                         \
  else if (type.compare(Functions::Expression<E, D, d, R, r, rC>::static_id())                                         \
           == 0) return Functions::Expression<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_2_3_2_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_2_3_3_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_3_1_1_SETTINGS(E, D, ddim, R, rdim, rCdim)                         \
  else if (type.compare(Functions::Expression<E, D, d, R, r, rC>::static_id())                                         \
           == 0) return Functions::Expression<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_3_1_2_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_3_1_3_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_3_2_1_SETTINGS(E, D, ddim, R, rdim, rCdim)                         \
  else if (type.compare(Functions::Expression<E, D, d, R, r, rC>::static_id())                                         \
           == 0) return Functions::Expression<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_3_2_2_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_3_2_3_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_3_3_1_SETTINGS(E, D, ddim, R, rdim, rCdim)                         \
  else if (type.compare(Functions::Expression<E, D, d, R, r, rC>::static_id())                                         \
           == 0) return Functions::Expression<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_3_3_2_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_3_3_3_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/

#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_1_1_1_CREATE(E, D, ddim, R, rdim, rCdim, settings)                 \
  else if (type.compare(Functions::Expression<E, D, d, R, r, rC>::static_id())                                         \
           == 0) return Functions::Expression<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_1_1_2_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_1_1_3_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_1_2_1_CREATE(E, D, ddim, R, rdim, rCdim, settings)                 \
  else if (type.compare(Functions::Expression<E, D, d, R, r, rC>::static_id())                                         \
           == 0) return Functions::Expression<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_1_2_2_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_1_2_3_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_1_3_1_CREATE(E, D, ddim, R, rdim, rCdim, settings)                 \
  else if (type.compare(Functions::Expression<E, D, d, R, r, rC>::static_id())                                         \
           == 0) return Functions::Expression<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_1_3_2_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_1_3_3_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_2_1_1_CREATE(E, D, ddim, R, rdim, rCdim, settings)                 \
  else if (type.compare(Functions::Expression<E, D, d, R, r, rC>::static_id())                                         \
           == 0) return Functions::Expression<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_2_1_2_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_2_1_3_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_2_2_1_CREATE(E, D, ddim, R, rdim, rCdim, settings)                 \
  else if (type.compare(Functions::Expression<E, D, d, R, r, rC>::static_id())                                         \
           == 0) return Functions::Expression<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_2_2_2_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_2_2_3_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_2_3_1_CREATE(E, D, ddim, R, rdim, rCdim, settings)                 \
  else if (type.compare(Functions::Expression<E, D, d, R, r, rC>::static_id())                                         \
           == 0) return Functions::Expression<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_2_3_2_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_2_3_3_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_3_1_1_CREATE(E, D, ddim, R, rdim, rCdim, settings)                 \
  else if (type.compare(Functions::Expression<E, D, d, R, r, rC>::static_id())                                         \
           == 0) return Functions::Expression<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_3_1_2_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_3_1_3_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_3_2_1_CREATE(E, D, ddim, R, rdim, rCdim, settings)                 \
  else if (type.compare(Functions::Expression<E, D, d, R, r, rC>::static_id())                                         \
           == 0) return Functions::Expression<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_3_2_2_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_3_2_3_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_3_3_1_CREATE(E, D, ddim, R, rdim, rCdim, settings)                 \
  else if (type.compare(Functions::Expression<E, D, d, R, r, rC>::static_id())                                         \
           == 0) return Functions::Expression<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_3_3_2_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Expression_3_3_3_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/

// ======================
// ==== Checkerboard ====
// ======================
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_1_1_1_STATIC_ID(E, D, ddim, R, rdim, rCdim)                      \
  , Functions::Checkerboard<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_1_1_2_STATIC_ID(E, D, ddim, R, rdim, rCdim)                      \
  , Functions::Checkerboard<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_1_1_3_STATIC_ID(E, D, ddim, R, rdim, rCdim)                      \
  , Functions::Checkerboard<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_1_2_1_STATIC_ID(E, D, ddim, R, rdim, rCdim)                      \
  , Functions::Checkerboard<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_1_2_2_STATIC_ID(E, D, ddim, R, rdim, rCdim)                      \
  , Functions::Checkerboard<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_1_2_3_STATIC_ID(E, D, ddim, R, rdim, rCdim)                      \
  , Functions::Checkerboard<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_1_3_1_STATIC_ID(E, D, ddim, R, rdim, rCdim)                      \
  , Functions::Checkerboard<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_1_3_2_STATIC_ID(E, D, ddim, R, rdim, rCdim)                      \
  , Functions::Checkerboard<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_1_3_3_STATIC_ID(E, D, ddim, R, rdim, rCdim)                      \
  , Functions::Checkerboard<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_2_1_1_STATIC_ID(E, D, ddim, R, rdim, rCdim)                      \
  , Functions::Checkerboard<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_2_1_2_STATIC_ID(E, D, ddim, R, rdim, rCdim)                      \
  , Functions::Checkerboard<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_2_1_3_STATIC_ID(E, D, ddim, R, rdim, rCdim)                      \
  , Functions::Checkerboard<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_2_2_1_STATIC_ID(E, D, ddim, R, rdim, rCdim)                      \
  , Functions::Checkerboard<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_2_2_2_STATIC_ID(E, D, ddim, R, rdim, rCdim)                      \
  , Functions::Checkerboard<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_2_2_3_STATIC_ID(E, D, ddim, R, rdim, rCdim)                      \
  , Functions::Checkerboard<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_2_3_1_STATIC_ID(E, D, ddim, R, rdim, rCdim)                      \
  , Functions::Checkerboard<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_2_3_2_STATIC_ID(E, D, ddim, R, rdim, rCdim)                      \
  , Functions::Checkerboard<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_2_3_3_STATIC_ID(E, D, ddim, R, rdim, rCdim)                      \
  , Functions::Checkerboard<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_3_1_1_STATIC_ID(E, D, ddim, R, rdim, rCdim)                      \
  , Functions::Checkerboard<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_3_1_2_STATIC_ID(E, D, ddim, R, rdim, rCdim)                      \
  , Functions::Checkerboard<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_3_1_3_STATIC_ID(E, D, ddim, R, rdim, rCdim)                      \
  , Functions::Checkerboard<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_3_2_1_STATIC_ID(E, D, ddim, R, rdim, rCdim)                      \
  , Functions::Checkerboard<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_3_2_2_STATIC_ID(E, D, ddim, R, rdim, rCdim)                      \
  , Functions::Checkerboard<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_3_2_3_STATIC_ID(E, D, ddim, R, rdim, rCdim)                      \
  , Functions::Checkerboard<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_3_3_1_STATIC_ID(E, D, ddim, R, rdim, rCdim)                      \
  , Functions::Checkerboard<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_3_3_2_STATIC_ID(E, D, ddim, R, rdim, rCdim)                      \
  , Functions::Checkerboard<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_3_3_3_STATIC_ID(E, D, ddim, R, rdim, rCdim)                      \
  , Functions::Checkerboard<E, D, ddim, R, rdim, rCdim>::static_id()

#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_1_1_1_SETTINGS(E, D, ddim, R, rdim, rCdim)                       \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_1_1_2_SETTINGS(E, D, ddim, R, rdim, rCdim)                       \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_1_1_3_SETTINGS(E, D, ddim, R, rdim, rCdim)                       \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_1_2_1_SETTINGS(E, D, ddim, R, rdim, rCdim)                       \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_1_2_2_SETTINGS(E, D, ddim, R, rdim, rCdim)                       \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_1_2_3_SETTINGS(E, D, ddim, R, rdim, rCdim)                       \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_1_3_1_SETTINGS(E, D, ddim, R, rdim, rCdim)                       \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_1_3_2_SETTINGS(E, D, ddim, R, rdim, rCdim)                       \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_1_3_3_SETTINGS(E, D, ddim, R, rdim, rCdim)                       \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_2_1_1_SETTINGS(E, D, ddim, R, rdim, rCdim)                       \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_2_1_2_SETTINGS(E, D, ddim, R, rdim, rCdim)                       \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_2_1_3_SETTINGS(E, D, ddim, R, rdim, rCdim)                       \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_2_2_1_SETTINGS(E, D, ddim, R, rdim, rCdim)                       \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_2_2_2_SETTINGS(E, D, ddim, R, rdim, rCdim)                       \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_2_2_3_SETTINGS(E, D, ddim, R, rdim, rCdim)                       \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_2_3_1_SETTINGS(E, D, ddim, R, rdim, rCdim)                       \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_2_3_2_SETTINGS(E, D, ddim, R, rdim, rCdim)                       \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_2_3_3_SETTINGS(E, D, ddim, R, rdim, rCdim)                       \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_3_1_1_SETTINGS(E, D, ddim, R, rdim, rCdim)                       \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_3_1_2_SETTINGS(E, D, ddim, R, rdim, rCdim)                       \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_3_1_3_SETTINGS(E, D, ddim, R, rdim, rCdim)                       \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_3_2_1_SETTINGS(E, D, ddim, R, rdim, rCdim)                       \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_3_2_2_SETTINGS(E, D, ddim, R, rdim, rCdim)                       \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_3_2_3_SETTINGS(E, D, ddim, R, rdim, rCdim)                       \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_3_3_1_SETTINGS(E, D, ddim, R, rdim, rCdim)                       \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_3_3_2_SETTINGS(E, D, ddim, R, rdim, rCdim)                       \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_3_3_3_SETTINGS(E, D, ddim, R, rdim, rCdim)                       \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::defaultSettings();

#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_1_1_1_CREATE(E, D, ddim, R, rdim, rCdim, settings)               \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_1_1_2_CREATE(E, D, ddim, R, rdim, rCdim, settings)               \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_1_1_3_CREATE(E, D, ddim, R, rdim, rCdim, settings)               \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_1_2_1_CREATE(E, D, ddim, R, rdim, rCdim, settings)               \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_1_2_2_CREATE(E, D, ddim, R, rdim, rCdim, settings)               \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_1_2_3_CREATE(E, D, ddim, R, rdim, rCdim, settings)               \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_1_3_1_CREATE(E, D, ddim, R, rdim, rCdim, settings)               \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_1_3_2_CREATE(E, D, ddim, R, rdim, rCdim, settings)               \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_1_3_3_CREATE(E, D, ddim, R, rdim, rCdim, settings)               \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_2_1_1_CREATE(E, D, ddim, R, rdim, rCdim, settings)               \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_2_1_2_CREATE(E, D, ddim, R, rdim, rCdim, settings)               \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_2_1_3_CREATE(E, D, ddim, R, rdim, rCdim, settings)               \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_2_2_1_CREATE(E, D, ddim, R, rdim, rCdim, settings)               \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_2_2_2_CREATE(E, D, ddim, R, rdim, rCdim, settings)               \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_2_2_3_CREATE(E, D, ddim, R, rdim, rCdim, settings)               \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_2_3_1_CREATE(E, D, ddim, R, rdim, rCdim, settings)               \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_2_3_2_CREATE(E, D, ddim, R, rdim, rCdim, settings)               \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_2_3_3_CREATE(E, D, ddim, R, rdim, rCdim, settings)               \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_3_1_1_CREATE(E, D, ddim, R, rdim, rCdim, settings)               \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_3_1_2_CREATE(E, D, ddim, R, rdim, rCdim, settings)               \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_3_1_3_CREATE(E, D, ddim, R, rdim, rCdim, settings)               \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_3_2_1_CREATE(E, D, ddim, R, rdim, rCdim, settings)               \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_3_2_2_CREATE(E, D, ddim, R, rdim, rCdim, settings)               \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_3_2_3_CREATE(E, D, ddim, R, rdim, rCdim, settings)               \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_3_3_1_CREATE(E, D, ddim, R, rdim, rCdim, settings)               \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_3_3_2_CREATE(E, D, ddim, R, rdim, rCdim, settings)               \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Checkerboard_3_3_3_CREATE(E, D, ddim, R, rdim, rCdim, settings)               \
  else if (type.compare(Functions::Checkerboard<E, D, d, R, r, rC>::static_id())                                       \
           == 0) return Functions::Checkerboard<E, D, d, R, r, rC>::create(settings);

// =====================
// ==== Spe10Model1 ====
// =====================
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_1_1_1_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_1_1_2_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_1_1_3_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_1_2_1_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_1_2_2_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_1_2_3_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_1_3_1_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_1_3_2_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_1_3_3_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_2_1_1_STATIC_ID(E, D, ddim, R, rdim, rCdim)                       \
  , Functions::Spe10Model1<E, D, ddim, R, rdim, rCdim>::static_id()
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_2_1_2_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_2_1_3_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_2_2_1_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_2_2_2_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_2_2_3_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_2_3_1_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_2_3_2_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_2_3_3_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_3_1_1_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_3_1_2_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_3_1_3_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_3_2_1_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_3_2_2_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_3_2_3_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_3_3_1_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_3_3_2_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_3_3_3_STATIC_ID(E, D, ddim, R, rdim, rCdim) /**/

#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_1_1_1_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_1_1_2_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_1_1_3_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_1_2_1_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_1_2_2_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_1_2_3_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_1_3_1_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_1_3_2_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_1_3_3_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_2_1_1_SETTINGS(E, D, ddim, R, rdim, rCdim)                        \
  else if (type.compare(Functions::Spe10Model1<E, D, d, R, r, rC>::static_id())                                        \
           == 0) return Functions::Spe10Model1<E, D, d, R, r, rC>::defaultSettings();
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_2_1_2_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_2_1_3_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_2_2_1_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_2_2_2_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_2_2_3_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_2_3_1_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_2_3_2_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_2_3_3_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_3_1_1_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_3_1_2_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_3_1_3_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_3_2_1_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_3_2_2_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_3_2_3_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_3_3_1_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_3_3_2_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_3_3_3_SETTINGS(E, D, ddim, R, rdim, rCdim) /**/

#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_1_1_1_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_1_1_2_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_1_1_3_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_1_2_1_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_1_2_2_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_1_2_3_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_1_3_1_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_1_3_2_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_1_3_3_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_2_1_1_CREATE(E, D, ddim, R, rdim, rCdim, settings)                \
  else if (type.compare(Functions::Spe10Model1<E, D, d, R, r, rC>::static_id())                                        \
           == 0) return Functions::Spe10Model1<E, D, d, R, r, rC>::create(settings);
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_2_1_2_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_2_1_3_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_2_2_1_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_2_2_2_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_2_2_3_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_2_3_1_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_2_3_2_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_2_3_3_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_3_1_1_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_3_1_2_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_3_1_3_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_3_2_1_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_3_2_2_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_3_2_3_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_3_3_1_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_3_3_2_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/
#define DUNE_STUFF_FUNCTIONS_HH_FUNCTION_Spe10Model1_3_3_3_CREATE(E, D, ddim, R, rdim, rCdim, settings) /**/

DUNE_STUFF_FUNCTIONS_GENERATE_DIM_RANGE(1)
DUNE_STUFF_FUNCTIONS_GENERATE_DIM_RANGE(2)
DUNE_STUFF_FUNCTIONS_GENERATE_DIM_RANGE(3)

#undef DUNE_STUFF_FUNCTIONS_GENERATE_CLASS
#undef DUNE_STUFF_FUNCTIONS_GENERATE_DIM_RANGE_COLS
#undef DUNE_STUFF_FUNCTIONS_GENERATE_DIM_RANGE

} // namespace Stuff
} // namespace Dune


#ifdef DUNE_STUFF_FUNCTIONS_TO_LIB
#define DUNE_STUFF_FUNCTIONS_LIST_DIMRANGE(etype, ddim)                                                                \
  DUNE_STUFF_FUNCTIONS_LIST_DIMRANGECOLS(etype, ddim, 1)                                                               \
  DUNE_STUFF_FUNCTIONS_LIST_DIMRANGECOLS(etype, ddim, 2)                                                               \
  DUNE_STUFF_FUNCTIONS_LIST_DIMRANGECOLS(etype, ddim, 3)

#define DUNE_STUFF_FUNCTIONS_LIST_DIMRANGECOLS(etype, ddim, rdim)                                                      \
  DUNE_STUFF_FUNCTIONS_LIST_DOMAINFIELDTYPES(etype, ddim, rdim, 1)                                                     \
  DUNE_STUFF_FUNCTIONS_LIST_DOMAINFIELDTYPES(etype, ddim, rdim, 2)                                                     \
  DUNE_STUFF_FUNCTIONS_LIST_DOMAINFIELDTYPES(etype, ddim, rdim, 3)

#define DUNE_STUFF_FUNCTIONS_LIST_DOMAINFIELDTYPES(etype, ddim, rdim, rcdim)                                           \
  DUNE_STUFF_FUNCTIONS_LIST_RANGEFIELDTYPES(etype, double, ddim, rdim, rcdim)

#define DUNE_STUFF_FUNCTIONS_LIST_RANGEFIELDTYPES(etype, dftype, ddim, rdim, rcdim)                                    \
  DUNE_STUFF_FUNCTIONS_LAST_EXPANSION(etype, dftype, ddim, double, rdim, rcdim)                                        \
  DUNE_STUFF_FUNCTIONS_LAST_EXPANSION(etype, dftype, ddim, long double, rdim, rcdim)

#define DUNE_STUFF_FUNCTIONS_LAST_EXPANSION(etype, dftype, ddim, rftype, rdim, rcdim)                                  \
  extern template class Dune::Stuff::Functions<etype, dftype, ddim, rftype, rdim, rcdim>;

#if HAVE_DUNE_GRID

DUNE_STUFF_FUNCTIONS_LIST_DIMRANGE(DuneStuffFunctionsInterfacesSGrid1dEntityType, 1)
DUNE_STUFF_FUNCTIONS_LIST_DIMRANGE(DuneStuffFunctionsInterfacesSGrid2dEntityType, 2)
DUNE_STUFF_FUNCTIONS_LIST_DIMRANGE(DuneStuffFunctionsInterfacesSGrid3dEntityType, 3)

DUNE_STUFF_FUNCTIONS_LIST_DIMRANGE(DuneStuffFunctionsInterfacesYaspGrid1dEntityType, 1)
DUNE_STUFF_FUNCTIONS_LIST_DIMRANGE(DuneStuffFunctionsInterfacesYaspGrid2dEntityType, 2)
DUNE_STUFF_FUNCTIONS_LIST_DIMRANGE(DuneStuffFunctionsInterfacesYaspGrid3dEntityType, 3)

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H

DUNE_STUFF_FUNCTIONS_LIST_DIMRANGE(DuneStuffFunctionsInterfacesAluSimplexGrid2dEntityType, 2)
DUNE_STUFF_FUNCTIONS_LIST_DIMRANGE(DuneStuffFunctionsInterfacesAluSimplexGrid3dEntityType, 3)
DUNE_STUFF_FUNCTIONS_LIST_DIMRANGE(DuneStuffFunctionsInterfacesAluCubeGrid3dEntityType, 3)

#endif // HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
#endif // HAVE_DUNE_GRID

#undef DUNE_STUFF_FUNCTIONS_LAST_EXPANSION
#undef DUNE_STUFF_FUNCTIONS_LIST_RANGEFIELDTYPES
#undef DUNE_STUFF_FUNCTIONS_LIST_DOMAINFIELDTYPES
#undef DUNE_STUFF_FUNCTIONS_LIST_DIMRANGECOLS
#undef DUNE_STUFF_FUNCTIONS_LIST_DIMRANGE

#endif // DUNE_STUFF_FUNCTIONS_TO_LIB
#endif // DUNE_STUFF_FUNCTIONS_HH
