// This file is part of the dune-xt-common project:
//   https://github.com/dune-community/dune-xt-common
// Copyright 2009-2018 dune-xt-common developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   René Fritze    (2018)
//   Tobias Leibner (2018)

// This one has to come first (includes the config.h)!
#include <dune/xt/common/test/main.hxx>

#include <dune/xt/common/crtp.hh>

template <class Traits>
class TestInterface;

struct FailTraits
{
  using BackendType = double;
  using derived_type = TestInterface<FailTraits>;
};

struct TestImp;

struct ImpTraits
{
  using BackendType = double;
  using derived_type = TestImp;
};

template <class Traits>
class TestInterface : public Dune::XT::Common::CRTPInterface<TestInterface<Traits>, Traits>
{
public:
  typedef typename Traits::BackendType BackendType;

  inline BackendType& backend()
  {
    CHECK_CRTP(this->as_imp().backend());
    return this->as_imp().backend();
  }

  inline const BackendType& backend() const
  {
    CHECK_CRTP(this->as_imp().backend());
    return this->as_imp().backend();
  }
}; // class ProvidesBackend


struct TestImp : public TestInterface<ImpTraits>
{
  using BackendType = ImpTraits::BackendType;
  BackendType foo_;

  inline BackendType& backend()
  {
    return foo_;
  }

  inline const BackendType& backend() const
  {
    return foo_;
  }
};

GTEST_TEST(crtp, fail)
{
  TestInterface<FailTraits> test_iface;
  EXPECT_THROW(
      {
        auto& st = test_iface.backend();
        static_cast<void>(st);
      },
      Dune::XT::Common::Exceptions::CRTP_check_failed);
}

GTEST_TEST(crtp, success)
{
  TestImp test_imp;
  auto& st = test_imp.backend();
  static_cast<void>(st);
}
