// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014, 2016 - 2017)
//   René Fritze     (2012 - 2016, 2018 - 2019)
//   Tobias Leibner  (2014, 2016, 2020)

#include <tuple>

#include <dune/xt/test/main.hxx>

#include <dune/common/dynmatrix.hh>
#include <dune/common/tupleutility.hh>
#include <dune/common/unused.hh>

#include <dune/xt/common/math.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/common/tuple.hh>

class A
{
public:
  typedef int MyFancyType;
};

class B
{
public:
  typedef int MyFancyType;
};

typedef std::tuple<A, B, B> MyTuple;

struct print_value
{
  template <class U, class V>
  static void run(std::ostream& out)
  {
    out << "(" << U::value << "," << V::value << ")" << std::endl;
  }

  template <class U, class V>
  static void run()
  {
    run<U, V>(std::cout);
  }
};

GTEST_TEST(TypeTransform, All)
{
  typedef TUPLE_TYPEDEFS_2_TUPLE(MyTuple, MyFancyType) MySubTupleType DUNE_UNUSED;
}

GTEST_TEST(Product, All)
{
  using namespace Dune::XT::Common::TupleProduct;
  typedef boost::mpl::vector<Int<1>, Int<2>> u_types;
  typedef boost::mpl::vector<Int<3>, Int<4>> v_types;
  typedef Combine<u_types, v_types, print_value>::Generate<> base_generator_type;
  base_generator_type::Run();
  base_generator_type::Run(std::cerr);
}

template <class arg1, class arg2>
struct tplA
{
  using type = arg1;
};

template <class arg1, class arg2>
struct tplB
{
  using type = arg2;
};

template <class Tuple>
void type_call()
{
  if constexpr (std::tuple_size_v<Tuple> == 3) {
    using TupleElement = typename Dune::XT::Common::list_content<Tuple>::template head<int, int>;
    using Type = typename TupleElement::type;
    static_assert(std::is_same<TupleElement, tplA<int, int>>::value, "");
    static_assert(std::is_same<Type, int>::value, "");
  } else if constexpr (std::tuple_size_v<Tuple> == 2) {
    using TupleElement = typename Dune::XT::Common::list_content<Tuple>::template head<int, int>;
    using Type = typename TupleElement::type;
    static_assert(std::is_same<TupleElement, tplB<int, int>>::value, "");
    static_assert(std::is_same<Type, int>::value, "");
  } else if constexpr (std::tuple_size_v<Tuple> == 1) {
    using TupleElement = typename Dune::XT::Common::list_content<Tuple>::template head<int, int>;
    static_assert(std::is_same<TupleElement, int>::value, "");
  }

  using Tail = typename Tuple::template tail_type<int, int>;
  type_call<Tail>();
}

template <>
void type_call<Dune::XT::Common::null_template_tuple>()
{
  // recursion terminator
}

GTEST_TEST(TemplateTuple, All)
{
  using namespace Dune::XT::Common;
  using tt = template_tuple<tplwrap<tplA>, tplwrap<tplB>, int>;
  type_call<tt>();
}
