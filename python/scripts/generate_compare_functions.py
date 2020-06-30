#!/usr/bin/env python
#
# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017)
#   René Fritze     (2016, 2018 - 2019)
#   Tobias Leibner  (2020)
# ~~~

import sys
import os
from string import Template

common_tpl = '''
template <Style style,
            class FirstType,
            class SecondType,
            class ToleranceType = typename MT<FirstType>::Eps>
  typename std::enable_if<internal::cmp_type_check<FirstType, SecondType, ToleranceType>::value, bool>::type ${id}(
      const FirstType& first,
      const SecondType& second,
      const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& rtol =
          Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType, style>::value(),
      const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& atol =
          Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType, style>::value())
  {
    return internal::Call<FirstType, SecondType, typename Dune::FloatCmp::EpsilonType<ToleranceType>::Type, style>::
        ${id}(first, second, rtol, atol);
  }

  template <class FirstType, class SecondType, class ToleranceType = typename MT<FirstType>::Eps>
  typename std::enable_if<internal::cmp_type_check<FirstType, SecondType, ToleranceType>::value, bool>::type ${id}(
      const FirstType& first,
      const SecondType& second,
      const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& rtol =
          Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType>::value(),
      const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& atol =
          Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType>::value())
  {
    return ${id}<Style::defaultStyle>(first, second, rtol, atol);
  }
'''

test_tpl = '''
  template <Dune::XT::Common::FloatCmp::Style style,
            class FirstType,
            class SecondType,
            class ToleranceType = typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps>
  typename std::
      enable_if<Dune::XT::Common::FloatCmp::internal::cmp_type_check<FirstType, SecondType, ToleranceType>::value,
                void>::type
          DXTC_EXPECT_FLOAT_${ID}(const FirstType& expected,
                                 const SecondType& actual,
                                 const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& rtol =
                                     Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType, style>::value(),
                                 const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& atol =
                                     Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType, style>::value())
  {
    const auto ${id} = Dune::XT::Common::FloatCmp::internal::
        Call<FirstType, SecondType, typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps, style>::${id};
    EXPECT_PRED4(${id}, expected, actual, rtol, atol);
  }

  template <class FirstType, class SecondType, class ToleranceType = typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps>
  typename std::enable_if<Dune::XT::Common::FloatCmp::internal::cmp_type_check<FirstType, SecondType, ToleranceType>::
                              value,
                          void>::type
      DXTC_EXPECT_FLOAT_${ID}(
          const FirstType& first,
          const SecondType& second,
          const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& rtol =
              Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType,
                                                         Dune::XT::Common::FloatCmp::Style::defaultStyle>::value(),
          const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& atol =
              Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType,
                                                         Dune::XT::Common::FloatCmp::Style::defaultStyle>::value())
  {
    DXTC_EXPECT_FLOAT_${ID}<Dune::XT::Common::FloatCmp::Style::defaultStyle>(first, second, rtol, atol);
  }
'''

bindir = os.path.dirname(os.path.abspath(__file__))
fn_test = os.path.join(bindir, '..', 'dune', 'xt', 'common', 'test', 'float_cmp_generated.hxx')
fn_common = os.path.join(bindir, '..', 'dune', 'xt', 'common', 'float_cmp_generated.hxx')
cmps = ['eq', 'ne', 'gt', 'lt', 'ge', 'le']
with open(fn_test, 'wt') as test_header, open(fn_common, 'wt') as common_header:
    for name in cmps:
        common_header.write(Template(common_tpl).substitute(id=name))
        test_header.write(Template(test_tpl).substitute(id=name, ID=name.upper()))
