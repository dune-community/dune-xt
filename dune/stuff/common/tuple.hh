#ifndef DUNE_STUFF_TUPLE_HH
#define DUNE_STUFF_TUPLE_HH

#include <dune/common/tuples.hh>
// #include <dune/fem/misc/femtuples.hh>


// ========== easy sub tuple extraction ===========
#define TMAX(t_, no_) (tuple_size<t_>::value <= (no_ + 1) ? 0 : no_)
#define TELE(t_, s_, no_) typename tuple_element<TMAX(t_, no_), t_>::type::s_
#define TCOND(t_, no_) (tuple_size<t_>::value <= no_)
#define SUBTUPLE(t_, s_)                                                                                                                                                                                                                                                                                                                                  \
  typename SelectType<TCOND(t_, 1),                                                                                                                                                                                                                                                                                                                       \
                      Dune::tuple<TELE(t_, s_, 0)>,                                                                                                                                                                                                                                                                                                       \
                      typename SelectType<TCOND(t_, 2),                                                                                                                                                                                                                                                                                                   \
                                          Dune::tuple<TELE(t_, s_, 0), TELE(t_, s_, 1)>,                                                                                                                                                                                                                                                                  \
                                          typename SelectType<TCOND(t_, 3),                                                                                                                                                                                                                                                                               \
                                                              Dune::tuple<TELE(t_, s_, 0),                                                                                                                                                                                                                                                                \
                                                                          TELE(t_, s_, 1),                                                                                                                                                                                                                                                                \
                                                                          TELE(t_, s_, 2)>,                                                                                                                                                                                                                                                               \
                                                              typename SelectType<TCOND(t_, 4),                                                                                                                                                                                                                                                           \
                                                                                  Dune::tuple<TELE(t_, s_, 0),                                                                                                                                                                                                                                            \
                                                                                              TELE(t_, s_, 1),                                                                                                                                                                                                                                            \
                                                                                              TELE(t_, s_, 2),                                                                                                                                                                                                                                            \
                                                                                              TELE(t_, s_, 3)>,                                                                                                                                                                                                                                           \
                                                                                  typename SelectType<TCOND(t_, 5),                                                                                                                                                                                                                                       \
                                                                                                      Dune::                                                                                                                                                                                                                                              \
                                                                                                          tuple<TELE(                                                                                                                                                                                                                                     \
                                                                                                                    t_,                                                                                                                                                                                                                                   \
                                                                                                                    s_,                                                                                                                                                                                                                                   \
                                                                                                                    0),                                                                                                                                                                                                                                   \
                                                                                                                TELE(                                                                                                                                                                                                                                     \
                                                                                                                    t_,                                                                                                                                                                                                                                   \
                                                                                                                    s_,                                                                                                                                                                                                                                   \
                                                                                                                    1),                                                                                                                                                                                                                                   \
                                                                                                                TELE(                                                                                                                                                                                                                                     \
                                                                                                                    t_,                                                                                                                                                                                                                                   \
                                                                                                                    s_,                                                                                                                                                                                                                                   \
                                                                                                                    2),                                                                                                                                                                                                                                   \
                                                                                                                TELE(                                                                                                                                                                                                                                     \
                                                                                                                    t_,                                                                                                                                                                                                                                   \
                                                                                                                    s_,                                                                                                                                                                                                                                   \
                                                                                                                    3),                                                                                                                                                                                                                                   \
                                                                                                                TELE(                                                                                                                                                                                                                                     \
                                                                                                                    t_,                                                                                                                                                                                                                                   \
                                                                                                                    s_,                                                                                                                                                                                                                                   \
                                                                                                                    4)>,                                                                                                                                                                                                                                  \
                                                                                                      typename SelectType<TCOND(                                                                                                                                                                                                                          \
                                                                                                                              t_,                                                                                                                                                                                                                         \
                                                                                                                              6),                                                                                                                                                                                                                         \
                                                                                                                          Dune::                                                                                                                                                                                                                          \
                                                                                                                              tuple<TELE(                                                                                                                                                                                                                 \
                                                                                                                                        t_,                                                                                                                                                                                                               \
                                                                                                                                        s_,                                                                                                                                                                                                               \
                                                                                                                                        0),                                                                                                                                                                                                               \
                                                                                                                                    TELE(                                                                                                                                                                                                                 \
                                                                                                                                        t_,                                                                                                                                                                                                               \
                                                                                                                                        s_,                                                                                                                                                                                                               \
                                                                                                                                        1),                                                                                                                                                                                                               \
                                                                                                                                    TELE(                                                                                                                                                                                                                 \
                                                                                                                                        t_,                                                                                                                                                                                                               \
                                                                                                                                        s_,                                                                                                                                                                                                               \
                                                                                                                                        2),                                                                                                                                                                                                               \
                                                                                                                                    TELE(                                                                                                                                                                                                                 \
                                                                                                                                        t_,                                                                                                                                                                                                               \
                                                                                                                                        s_,                                                                                                                                                                                                               \
                                                                                                                                        3),                                                                                                                                                                                                               \
                                                                                                                                    TELE(                                                                                                                                                                                                                 \
                                                                                                                                        t_,                                                                                                                                                                                                               \
                                                                                                                                        s_,                                                                                                                                                                                                               \
                                                                                                                                        4),                                                                                                                                                                                                               \
                                                                                                                                    TELE(                                                                                                                                                                                                                 \
                                                                                                                                        t_,                                                                                                                                                                                                               \
                                                                                                                                        s_,                                                                                                                                                                                                               \
                                                                                                                                        5)>,                                                                                                                                                                                                              \
                                                                                                                          typename SelectType<TCOND(                                                                                                                                                                                                      \
                                                                                                                                                  t_,                                                                                                                                                                                                     \
                                                                                                                                                  7),                                                                                                                                                                                                     \
                                                                                                                                              Dune::tuple<TELE(t_, s_, 0), TELE(t_, s_, 1), TELE(t_, s_, 2), TELE(t_, s_, 3), TELE(t_, s_, 4), TELE(t_, s_, 5), TELE(t_, s_, 6)>,                                                                         \
                                                                                                                                              typename SelectType<TCOND(                                                                                                                                                                                  \
                                                                                                                                                                      t_,                                                                                                                                                                                 \
                                                                                                                                                                      8),                                                                                                                                                                                 \
                                                                                                                                                                  Dune::tuple<TELE(t_, s_, 0), TELE(t_, s_, 1), TELE(t_, s_, 2), TELE(t_, s_, 3), TELE(t_, s_, 4), TELE(t_, s_, 5), TELE(t_, s_, 6), TELE(t_, s_, 7)>,                                    \
                                                                                                                                                                  typename SelectType<TCOND(                                                                                                                                                              \
                                                                                                                                                                                          t_,                                                                                                                                                             \
                                                                                                                                                                                          9),                                                                                                                                                             \
                                                                                                                                                                                      Dune::                                                                                                                                                              \
                                                                                                                                                                                          tuple<TELE(t_, s_, 0), TELE(t_, s_, 1), TELE(t_, s_, 2), TELE(t_, s_, 3), TELE(t_, s_, 4), TELE(t_, s_, 5), TELE(t_, s_, 6), TELE(t_, s_, 7), TELE(t_, s_, 8)>, \
                                                                                                                                                                                      Dune::                                                                                                                                                              \
                                                                                                                                                                                          tuple<TELE(                                                                                                                                                     \
                                                                                                                                                                                                    t_,                                                                                                                                                   \
                                                                                                                                                                                                    s_,                                                                                                                                                   \
                                                                                                                                                                                                    0),                                                                                                                                                   \
                                                                                                                                                                                                TELE(                                                                                                                                                     \
                                                                                                                                                                                                    t_,                                                                                                                                                   \
                                                                                                                                                                                                    s_,                                                                                                                                                   \
                                                                                                                                                                                                    1),                                                                                                                                                   \
                                                                                                                                                                                                TELE(                                                                                                                                                     \
                                                                                                                                                                                                    t_,                                                                                                                                                   \
                                                                                                                                                                                                    s_,                                                                                                                                                   \
                                                                                                                                                                                                    2),                                                                                                                                                   \
                                                                                                                                                                                                TELE(                                                                                                                                                     \
                                                                                                                                                                                                    t_, s_, 3),                                                                                                                                           \
                                                                                                                                                                                                TELE(                                                                                                                                                     \
                                                                                                                                                                                                    t_,                                                                                                                                                   \
                                                                                                                                                                                                    s_,                                                                                                                                                   \
                                                                                                                                                                                                    4),                                                                                                                                                   \
                                                                                                                                                                                                TELE(                                                                                                                                                     \
                                                                                                                                                                                                    t_,                                                                                                                                                   \
                                                                                                                                                                                                    s_,                                                                                                                                                   \
                                                                                                                                                                                                    5),                                                                                                                                                   \
                                                                                                                                                                                                TELE(                                                                                                                                                     \
                                                                                                                                                                                                    t_,                                                                                                                                                   \
                                                                                                                                                                                                    s_,                                                                                                                                                   \
                                                                                                                                                                                                    6),                                                                                                                                                   \
                                                                                                                                                                                                TELE(                                                                                                                                                     \
                                                                                                                                                                                                    t_,                                                                                                                                                   \
                                                                                                                                                                                                    s_,                                                                                                                                                   \
                                                                                                                                                                                                    7),                                                                                                                                                   \
                                                                                                                                                                                                TELE(                                                                                                                                                     \
                                                                                                                                                                                                    t_,                                                                                                                                                   \
                                                                                                                                                                                                    s_,                                                                                                                                                   \
                                                                                                                                                                                                    8)>>::Type>::Type>::                                                                                                                                  \
                                                                                                                              Type>::                                                                                                                                                                                                                     \
                                                                                                          Type>::                                                                                                                                                                                                                                         \
                                                                                      Type>::Type>::Type>::Type>::Type
// "SUBTUPLE" extracts types of the elements of a tuple and
// creates a new tuple with these extracted types
//
// example:
//
// Let two classes be given:
//
// class A{ public: typedef int MyFancyType; };
// class B{ public: typedef int MyFancyType; };
//
// Define a tuple, i.e.
//
// typedef Dune::tuple< A, B, B >
//   MyTuple;
//
// Instead of writing
//
// typedef Dune::tuple< typename tuple_element<0,MyTuple>::type::MyFancyType,
//                      typename tuple_element<1,MyTuple>::type::MyFancyType,
//                      typename tuple_element<2,MyTuple>::type::MyFancyType >
//   MySubTupleType;
//
// just write
//
// typedef SUBTUPLE( MyTuple, MyFancyType )
//   MySubTupleType;
//
// without (general) restriction to the size
// of the tuple MyTuple. Enjoy it!
// ===============================================


namespace Dune {
namespace Stuff {
namespace Common {

//! split our function wrapper into sep pointers for Datawriter input
template <class T1, class T2 = T1, class T3 = T2, class T4 = T3>
struct TupleSerializer
{
  typedef Dune::
      tuple<const typename T1::DiscreteVelocityFunctionType*, const typename T1::DiscretePressureFunctionType*,
            const typename T2::DiscreteVelocityFunctionType*, const typename T2::DiscretePressureFunctionType*,
            const typename T3::DiscreteVelocityFunctionType*, const typename T3::DiscretePressureFunctionType*,
            const typename T4::DiscreteVelocityFunctionType*, const typename T4::DiscretePressureFunctionType*>
          TupleType;

  static TupleType& getTuple(T1& t1, T2& t2, T3& t3, T4& t4)
  {
    // yay for dangling pointers, but using a local static here fubared sequential runs with diff grid
    TupleType* t = new TupleType(&(t1.discreteVelocity()),
                                 &(t1.discretePressure()),
                                 &(t2.discreteVelocity()),
                                 &(t2.discretePressure()),
                                 &(t3.discreteVelocity()),
                                 &(t3.discretePressure()),
                                 &(t4.discreteVelocity()),
                                 &(t4.discretePressure()));

    return *t;
  } // getTuple

  static TupleType& getTuple(T1& t1, T2& t2, T3& t3)
  {
    // yay for dangling pointers, but using a local static here fubared sequential runs with diff grid
    TupleType* t = new TupleType(&(t1.discreteVelocity()),
                                 &(t1.discretePressure()),
                                 &(t2.discreteVelocity()),
                                 &(t2.discretePressure()),
                                 &(t3.discreteVelocity()),
                                 &(t3.discretePressure()),
                                 nullptr,
                                 nullptr);

    return *t;
  } // getTuple

  static TupleType& getTuple(T1& t1, T2& t2)
  {
    // yay for dangling pointers, but using a local static here fubared sequential runs with diff grid
    TupleType* t = new TupleType(&(t1.discreteVelocity()),
                                 &(t1.discretePressure()),
                                 &(t2.discreteVelocity()),
                                 &(t2.discretePressure()),
                                 nullptr,
                                 nullptr,
                                 nullptr,
                                 nullptr);

    return *t;
  } // getTuple

  static TupleType& getTuple(T1& t1)
  {
    // yay for dangling pointers, but using a local static here fubared sequential runs with diff grid
    TupleType* t = new TupleType(
        &(t1.discreteVelocity()), &(t1.discretePressure()), nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);

    return *t;
  } // getTuple
};

template <class T1, class T2 = T1, class T3 = T1, class T4 = T1, class T5 = T1, class T6 = T1, class T7 = T1,
          class T8 = T1, class T9 = T1>
struct FullTuple : public Dune::tuple<T1, T2, T3, T4, T5, T6, T7, T8, T9>
{
  FullTuple(T1 t1, T2 t2, T3 t3, T4 t4, T5 t5, T6 t6, T7 t7, T8 t8, T9 t9)
    : Dune::tuple<T1, T2, T3, T4, T5, T6, T7, T8, T9>(t1, t2, t3, t4, t5, t6, t7, t8, t9)
  {
  }
};

} // namespace Common
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_TUPLE_HH

/** Copyright (c) 2012, Rene Milk
   * All rights reserved.
   *
   * Redistribution and use in source and binary forms, with or without
   * modification, are permitted provided that the following conditions are met:
   *
   * 1. Redistributions of source code must retain the above copyright notice, this
   *    list of conditions and the following disclaimer.
   * 2. Redistributions in binary form must reproduce the above copyright notice,
   *    this list of conditions and the following disclaimer in the documentation
   *    and/or other materials provided with the distribution.
   *
   * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
   * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
   * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
   * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
   * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
   * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
   *
   * The views and conclusions contained in the software and documentation are those
   * of the authors and should not be interpreted as representing official policies,
   * either expressed or implied, of the FreeBSD Project.
   **/
