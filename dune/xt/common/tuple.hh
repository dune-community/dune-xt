// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012, 2014 - 2017)
//   René Fritze     (2010 - 2013, 2015 - 2016, 2018 - 2019)
//   Stefan Girke    (2012)
//   Tobias Leibner  (2014, 2016 - 2018, 2020)

#ifndef DUNE_XT_COMMON_TUPLE_HH
#define DUNE_XT_COMMON_TUPLE_HH

#include <utility>

#include <dune/common/typetraits.hh>


#include <boost/tuple/tuple.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/deref.hpp>
#include <boost/mpl/begin.hpp>
#include <boost/mpl/end.hpp>
#include <boost/mpl/next.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/type_traits.hpp>

#define TMAX(t_, no_) (std::tuple_size<t_>::value >= (no_ + 1) ? no_ : 0)
#define TELE(t_, s_, no_) typename std::tuple_element<TMAX(t_, no_), t_>::type::s_

template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8, class T9, size_t Size>
struct RightTrimTuple
{
  typedef std::tuple<T1, T2, T3, T4, T5, T6, T7, T8, T9> type;
};

template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8, class T9>
struct RightTrimTuple<T1, T2, T3, T4, T5, T6, T7, T8, T9, 9>
{
  typedef std::tuple<T1, T2, T3, T4, T5, T6, T7, T8, T9> type;
};

template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8, class T9>
struct RightTrimTuple<T1, T2, T3, T4, T5, T6, T7, T8, T9, 8>
{
  typedef std::tuple<T1, T2, T3, T4, T5, T6, T7, T8> type;
};

template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8, class T9>
struct RightTrimTuple<T1, T2, T3, T4, T5, T6, T7, T8, T9, 7>
{
  typedef std::tuple<T1, T2, T3, T4, T5, T6, T7> type;
};

template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8, class T9>
struct RightTrimTuple<T1, T2, T3, T4, T5, T6, T7, T8, T9, 6>
{
  typedef std::tuple<T1, T2, T3, T4, T5, T6> type;
};

template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8, class T9>
struct RightTrimTuple<T1, T2, T3, T4, T5, T6, T7, T8, T9, 5>
{
  typedef std::tuple<T1, T2, T3, T4, T5> type;
};

template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8, class T9>
struct RightTrimTuple<T1, T2, T3, T4, T5, T6, T7, T8, T9, 4>
{
  typedef std::tuple<T1, T2, T3, T4> type;
};

template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8, class T9>
struct RightTrimTuple<T1, T2, T3, T4, T5, T6, T7, T8, T9, 3>
{
  typedef std::tuple<T1, T2, T3> type;
};

template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8, class T9>
struct RightTrimTuple<T1, T2, T3, T4, T5, T6, T7, T8, T9, 2>
{
  typedef std::tuple<T1, T2> type;
};

template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8, class T9>
struct RightTrimTuple<T1, T2, T3, T4, T5, T6, T7, T8, T9, 1>
{
  using type = std::tuple<T1>;
};

#define TUPLE_TYPEDEFS_2_TUPLE(t_, s_)                                                                                 \
  typename RightTrimTuple<TELE(t_, s_, 0),                                                                             \
                          TELE(t_, s_, 1),                                                                             \
                          TELE(t_, s_, 2),                                                                             \
                          TELE(t_, s_, 3),                                                                             \
                          TELE(t_, s_, 4),                                                                             \
                          TELE(t_, s_, 5),                                                                             \
                          TELE(t_, s_, 6),                                                                             \
                          TELE(t_, s_, 7),                                                                             \
                          TELE(t_, s_, 8),                                                                             \
                          std::tuple_size<t_>::value>::type

/**
 * @def TUPLE_TYPEDEFS_2_TUPLE( t_, s_ )
 *
 * @brief extracts types of the elements of a std::tuple and
 * creates a new std::tuple with these extracted types.
 *
 * first argument: the tuple type
 * second argument: the typedef of the elements of the tuple
 *
 * example:
 * Let two classes be given:
 *
 * @code
 * class A
 * {
 *   public:
 *   typedef int MyFancyType;
 * };
 *
 * class B
 * {
 *   public:
 *   typedef int MyFancyType;
 * };
 * @endcode
 *
 * Define a tuple, i.e.
 *
 * @code
 * typedef std::tuple< A, B, B >
 *   MyTuple;
 * @endcode
 *
 * Instead of writing
 *
 * @code
 * typedef std::tuple< typename tuple_element<0,MyTuple>::type::MyFancyType,
 *                      typename tuple_element<1,MyTuple>::type::MyFancyType,
 *                      typename tuple_element<2,MyTuple>::type::MyFancyType >
 *   MySubTupleType;
 * @endcode
 *
 * just write
 *
 * @code
 * typedef TUPLE_TYPEDEFS_2_TUPLE( MyTuple, MyFancyType )
 *   MySubTupleType;
 * @endcode
 *
 * without (general) restriction to the size
 * of the tuple @c MyTuple. Enjoy it!
 */

namespace Dune::XT::Common {

// reduced from
// http://stackoverflow.com/questions/1492204/is-it-possible-to-generate-types-with-all-combinations-of-template-arguments
namespace TupleProduct {
using boost::is_same;
using boost::mpl::begin;
using boost::mpl::deref;
using boost::mpl::end;
using boost::mpl::if_;
using boost::mpl::next;

struct end_of_recursion_tag
{
  template <class... Args>
  static void Run(Args&&... /*args*/)
  {}
};

template <class UTypes, // Forward Sequence, e.g. boost::mpl::vector
          class VTypes, // Forward Sequence, e.g. boost::mpl::vector
          class TestFunc // class type that has a nested templated run() member function
          >
struct Combine
{
  // forward declaration
  template <class UIterator, class VIterator>
  struct Generate;

  // this class implements recursion body
  template <class UIterator, class VIterator>
  struct Next
  {
    // u_begin is not necessary ;)
    // it would be cheaper not to pre-declare all of them since we force evaluation
    // however this dramatically increase the readability
    using v_begin = typename begin<VTypes>::type;

    using u_end = typename end<UTypes>::type;
    using v_end = typename end<VTypes>::type;

    using u_next = typename next<UIterator>::type;
    using v_next = typename next<VIterator>::type;

    typedef
        typename if_<boost::is_same<v_next, v_end>,
                     typename if_<boost::is_same<u_next, u_end>, end_of_recursion_tag, Generate<u_next, v_begin>>::type,
                     Generate<UIterator, v_next>>::type type;
  };

  //  this class run test on generated types in thos round and go to next*/
  template <class UIterator = typename begin<UTypes>::type, class VIterator = typename begin<VTypes>::type>
  struct Generate
  {
    //  generate <<next>> target type
    typedef typename Next<UIterator, VIterator>::type next_type;

    static void Run()
    {
      // test on the generated types of this round of recursion
      TestFunc::template run<typename deref<UIterator>::type, typename deref<VIterator>::type>();

      // go to the next round of recursion
      next_type::Run();
    }

    template <class... Args>
    static void Run(Args&&... args)
    {
      // test on the generated types of this round of recursion
      TestFunc::template run<typename deref<UIterator>::type, typename deref<VIterator>::type>(
          std::forward<Args>(args)...);

      // go to the next round of recursion
      next_type::Run(std::forward<Args>(args)...);
    }
  };
};

} // namespace TupleProduct

//! from https://stackoverflow.com/questions/17424477/implementation-c14-make-integer-sequence/17426611
template <size_t... Is>
struct index_sequence
{
  using type = index_sequence;
};


namespace internal {


template <class S1, class S2>
struct Concat;
template <class S1, class S2>
using Concat_t = typename Concat<S1, S2>::type;

template <size_t... I1, size_t... I2>
struct Concat<index_sequence<I1...>, index_sequence<I2...>>
{
  using type = index_sequence<I1..., (sizeof...(I1) + I2)...>;
};


} // namespace internal


template <size_t N>
struct make_index_sequence;
template <size_t N>
using make_index_sequence_t = typename make_index_sequence<N>::type;

template <size_t N>
struct make_index_sequence
{
  using type = internal::Concat_t<make_index_sequence_t<N / 2>, make_index_sequence_t<N - N / 2>>;
};

// break conditions
template <>
struct make_index_sequence<1>
{
  using type = index_sequence<0>;
};
template <>
struct make_index_sequence<0>
{
  using type = index_sequence<>;
};

//! T_aliased< T, Index > is always the type T, no matter what index is
template <typename T, std::size_t index>
using T_aliased = T;

//! make_identical_tuple< T, N >::type is a std::tuple< T, ... , T > with a length of N
template <typename T, std::size_t N, typename I = typename make_index_sequence<N>::type>
struct make_identical_tuple;

template <typename T, std::size_t N, std::size_t... Indices>
struct make_identical_tuple<T, N, index_sequence<Indices...>>
{
  using type = std::tuple<T_aliased<T, Indices>...>;

  static type create(const T& t)
  {
    return type(T_aliased<T, Indices>{t}...);
  }
};

//! implementation of std::tuple_element with O(log(N)) template recursion depth
//! from https://stackoverflow.com/questions/18593057/stdtuple-element-need-deep-template-instantination
namespace internal {
template <std::size_t>
struct Any
{
  Any(...) {}
};

template <typename T>
struct wrapper
{
  using type = T;
};

template <std::size_t... Is>
struct get_nth_helper
{
  template <typename T>
  static auto deduce(Any<Is>..., wrapper<T>, ...) -> wrapper<T>;
};

template <std::size_t... Is, typename... Ts>
auto deduce_seq(index_sequence<Is...>, wrapper<Ts>... pp) -> decltype(get_nth_helper<Is...>::deduce(pp...));
} // namespace internal

template <std::size_t N, class Tuple>
struct tuple_element;

template <std::size_t N, class... Ts>
struct tuple_element<N, std::tuple<Ts...>>
{
  using wrapped_type = decltype(internal::deduce_seq(make_index_sequence_t<N>{}, internal::wrapper<Ts>()...));
  using type = typename wrapped_type::type;
};

/** These classes allow recursively visiting a typelist that may contain templates names
 *  So if you ever wanted to iterate over "AllDirichletBoundaryInfo, AllNeumannBoundaryInfo, etc"
 *  and then input another type as tpl arg, sthis is for you.
 *
 * **/

//! use this to wrap template classes into the template_tuple
template <template <class...> class B>
struct tplwrap
{};

template <class... Ms>
class template_tuple;

namespace internal {
template <class B, class...>
struct from_tplwrap
{
  using type = B;
};

template <template <class...> class B, class... T>
struct from_tplwrap<tplwrap<B>, T...>
{
  using type = B<T...>;
};


template <class... Templates>
struct type_help
{
  template <class... Parameters>
  using type = boost::tuple<typename from_tplwrap<Templates, Parameters...>::type...>;
};

template <class... Args>
struct type_help<boost::tuple<Args...>>
{
  template <class...>
  using type = boost::tuple<Args...>;
};

template <class... Args>
struct tail_hlp
{
  using type = template_tuple<Args...>;
};

template <class... Args>
struct tail_hlp<boost::tuple<Args...>>
{
  using type = template_tuple<Args...>;
};

template <class A>
struct head_hlp
{
  using type = A;
};

template <class A, class B>
struct head_hlp<boost::tuples::cons<A, B>>
{
  using type = A;
};
} // namespace internal

template <class... WrapperOrTypes>
class template_tuple
{

public:
  template <typename... T>
  using type = typename internal::type_help<WrapperOrTypes...>::template type<T...>;
  template <typename... T>
  using tail_type = typename internal::tail_hlp<typename type<T...>::tail_type>::type;
  template <typename... T>
  using head_type = typename internal::head_hlp<typename type<T...>::head_type>::type;
};

using null_template_tuple = template_tuple<boost::tuples::null_type>;

} // namespace Dune::XT::Common

#endif // DUNE_XT_COMMON_TUPLE_HH
