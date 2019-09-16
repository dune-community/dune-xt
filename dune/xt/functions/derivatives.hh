// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_XT_FUNCTIONS_DERIVATIVES_HH
#define DUNE_XT_FUNCTIONS_DERIVATIVES_HH

#include <dune/xt/functions/base/derivatives-of-element-functions.hh>
#include <dune/xt/functions/interfaces/element-functions.hh>

namespace Dune {
namespace XT {
namespace Functions {


template <class E, class R>
DivergenceElementFunction<ElementFunctionInterface<E, E::dimension, 1, R>>
    divergence(ElementFunctionInterface<E, E::dimension, 1, R>& func)
{
  return DivergenceElementFunction<ElementFunctionInterface<E, E::dimension, 1, R>>(func);
}


template <class E, class R>
GradientElementFunction<ElementFunctionInterface<E, 1, 1, R>> gradient(ElementFunctionInterface<E, 1, 1, R>& func)
{
  return GradientElementFunction<ElementFunctionInterface<E, 1, 1, R>>(func);
}


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_DERIVATIVES_HH
