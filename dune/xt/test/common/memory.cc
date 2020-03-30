// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   René Fritze     (2016, 2018 - 2019)

#include <dune/xt/test/main.hxx>

#include <dune/common/dynmatrix.hh>
#include <dune/common/tupleutility.hh>

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/fvector.hh>

using namespace Dune::XT::Common;


struct ScopeTest : public testing::Test
{
  using T = int;
  static constexpr T constant = 1;

  template <class P>
  void deref(P& p)
  {
    EXPECT_NE(p, nullptr);
    auto g = *p;
    EXPECT_EQ(constant, g);
  }

  template <class P>
  void access(const P& p)
  {
    auto g = p.access();
    EXPECT_EQ(constant, g);
  }

  template <template <class F> class Provider, class P>
  void scope(P& p)
  {
    Provider<T> shared_provider(p);
    deref(p);
  }

  void check_shared()
  {
    auto shared = std::make_shared<T>(constant);
    scope<ConstStorageProvider>(shared);
    deref(shared);
    scope<StorageProvider>(shared);
    deref(shared);
  }

  void check_const()
  {
    using CSP = ConstStorageProvider<T>;
    access(CSP{new T(constant)});
    CSP{};
    T e{constant};
    CSP{e};
  }
};
constexpr typename ScopeTest::T ScopeTest::constant;

TEST_F(ScopeTest, All)
{
  this->check_shared();
  this->check_const();
}
