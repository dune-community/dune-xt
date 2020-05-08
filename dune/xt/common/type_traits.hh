// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2018)
//   René Fritze     (2012 - 2019)
//   Tobias Leibner  (2016, 2018 - 2020)

#ifndef DUNE_XT_COMMON_TYPE_TRAITS_HH
#define DUNE_XT_COMMON_TYPE_TRAITS_HH

#ifdef __GNUC__
#  include <cxxabi.h>
#endif

#include <complex>
#include <memory>
#include <type_traits>
#include <iostream>

#include <dune/common/bigunsignedint.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/unused.hh>

#include <dune/xt/common/exceptions.hh>

/** use this to define Typename specializations in the GLOBAL namespace ONLY **/
#define DUNE_XT_COMMON_TYPENAME(NAME)                                                                                  \
  namespace Dune {                                                                                                     \
  namespace XT {                                                                                                       \
  namespace Common {                                                                                                   \
  template <>                                                                                                          \
  struct Typename<NAME>                                                                                                \
  {                                                                                                                    \
    static const char* value()                                                                                         \
    {                                                                                                                  \
      return #NAME;                                                                                                    \
    }                                                                                                                  \
  };                                                                                                                   \
  }                                                                                                                    \
  }                                                                                                                    \
  }

/**
 * \brief Helper macro to be used before DXTC_has_typedef.
 * \see   DXTC_has_typedef
 *
 *        Taken from
 *        http://stackoverflow.com/questions/7834226/detecting-typedef-at-compile-time-template-metaprogramming
 */
#define DXTC_has_typedef_initialize_once(tpdef)                                                                        \
  template <typename T_local>                                                                                          \
  struct DXTC_has_typedef_helper_##tpdef                                                                               \
  {                                                                                                                    \
    template <typename TT_local>                                                                                       \
    struct void_                                                                                                       \
    {                                                                                                                  \
      typedef void type;                                                                                               \
    };                                                                                                                 \
                                                                                                                       \
    template <typename TT_local, typename = void>                                                                      \
    struct helper                                                                                                      \
    {                                                                                                                  \
      static const bool value = false;                                                                                 \
    };                                                                                                                 \
                                                                                                                       \
    template <typename TT_local>                                                                                       \
    struct helper<TT_local, typename void_<typename TT_local::tpdef>::type>                                            \
    {                                                                                                                  \
      static const bool value = true;                                                                                  \
    };                                                                                                                 \
                                                                                                                       \
    static const bool value = helper<T_local>::value;                                                                  \
  }

/**
  * \brief Macro to statically check a type for a typedef.
  *
  *        To check if a type Foo defines a typedef Bar, put this code somewhere, where a generated helper struct
  *        may be defined (obviously only once for each typedef name!, note the trailing `;`!):
\code
DXTC_has_typedef_initialize_once(Bar);
\endcode
  *        You can then check for Bar (this will give you a static const bool):
\code
DXTC_has_typedef(Bar)< Foo >::value
\endcode
  */
#define DXTC_has_typedef(tpdf) DXTC_has_typedef_helper_##tpdf

/**
 * \brief Helper macro to be used before DXTC_has_static_member.
 *
 *        Taken from http://stackoverflow.com/questions/11927032/sfinae-check-for-static-member-using-decltype
 */
#define DXTC_has_static_member_initialize_once(mmbr)                                                                   \
  template <class T_local>                                                                                             \
  struct DXTC_has_static_member_helper_##mmbr                                                                          \
  {                                                                                                                    \
    template <class TT_local,                                                                                          \
              class = typename std::enable_if<!std::is_member_pointer<decltype(&TT_local::mmbr)>::value>::type>        \
    static std::true_type helper(int);                                                                                 \
                                                                                                                       \
    template <class>                                                                                                   \
    static std::false_type helper(...);                                                                                \
                                                                                                                       \
    static const bool value = decltype(helper<T_local>(0))::value;                                                     \
  }

/**
  * \brief Macro to statically check a class or struct for a static member.
  *
  *        To check if a class or struct Foo has a static member bar, put this code somewhere, where a generated helper
  *        struct may be defined (obviously only once for each member name!, note the trailing `;`!):
\code
DXTC_has_static_member_initialize_once(bar);
\endcode
  *        You can then check for bar (this will give you a static const bool):
\code
DXTC_has_static_member(bar)< Foo >::value
\endcode
  */
#define DXTC_has_static_member(mmbr) DXTC_has_static_member_helper_##mmbr


/**
 * \brief Helper macro to be used before DXTC_has_method.
 *
 *        Inspired by https://stackoverflow.com/questions/40329684/sfinae-c-method-check
 */
#define DXTC_has_method_initialize_once(mthd_nm)                                                                       \
  template <class T_local>                                                                                             \
  struct DXTC_has_method_helper_##mthd_nm                                                                              \
  {                                                                                                                    \
    template <class TT_local>                                                                                          \
    static std::true_type helper(decltype(&TT_local::mthd_nm));                                                        \
                                                                                                                       \
    template <class>                                                                                                   \
    static std::false_type helper(...);                                                                                \
                                                                                                                       \
    static const bool value = decltype(helper<T_local>(0))::value;                                                     \
  }


/**
  * \brief Macro to statically check a class or struct for a method.
  *
  *        To check if a class or struct Foo has a method bar, put this code somewhere, where a generated helper
  *        struct may be defined (obviously only once for each method name!, note the trailing `;`!):
\code
DXTC_has_method_initialize_once(bar);
\endcode
  *        You can then check for bar (this will give you a static const bool):
\code
DXTC_has_method(bar)< Foo >::value
\endcode
  */
#define DXTC_has_method(mthd_nm) DXTC_has_method_helper_##mthd_nm


namespace Dune {


// forward

template <class K, int SIZE>
class FieldVector;

template <class K, int ROWS, int COLS>
class FieldMatrix;


namespace XT {
namespace Common {


// forward

template <class K, int SIZE>
class FieldVector;

template <class K, int ROWS, int COLS>
class FieldMatrix;

template <class VecType>
struct VectorAbstraction;

template <class MatType>
struct MatrixAbstraction;


inline std::string demangle_typename(const std::string& mangled_name)
{
#ifdef __GNUC__
  return abi::__cxa_demangle(mangled_name.c_str(), nullptr, nullptr, nullptr);
#else // ifdef __GNUC__
  return mangled_name;
#endif // ifdef __GNUC__
}


//! demangles typeid
template <class T>
std::string demangled_type_id(T& obj)
{
  return demangle_typename(typeid(obj).name());
} // demangled_type_id


//! create output for demangled typeid
template <class T>
void real_type_id(T& obj, std::string name = "", size_t maxlevel = 10000)
{
  std::cout << name << (name == "" ? "" : "'s type is ") << highlight_template(demangled_type_id(obj), maxlevel)
            << std::endl;
}


//! an extensible mechanism to assign "cleartext" names to types
template <typename T>
struct Typename
{
  static std::string value(bool fail_wo_typeid = false)
  {
#if defined(__GNUC__) && defined(__GXX_RTTI)
    DUNE_UNUSED_PARAMETER(fail_wo_typeid);
    return demangle_typename(typeid(T).name());
#else
    if (fail_wo_typeid)
      DUNE_THROW(InvalidStateException, "typeid translation requested to fail with missing rtti");
    return "unknown";
#endif
  }
};

template <class T>
std::string get_typename(const T&)
{
  return Typename<T>::value();
}

template <class T>
std::string get_template_basename(const T&)
{
  const auto str = Typename<T>::value();
  const auto r = str.find("<");
  return str.substr(0, r);
}

template <class T, class Ptr = void>
struct is_smart_ptr
{
  static const bool value = false;
  typedef T type;
};

template <class T>
struct is_smart_ptr<T, typename std::enable_if<std::is_same<std::unique_ptr<typename T::element_type>, T>::value>::type>
{
  static const bool value = true;
  typedef T type;
};

template <class T>
struct is_smart_ptr<T, typename std::enable_if<std::is_same<std::shared_ptr<typename T::element_type>, T>::value>::type>
{
  static const bool value = true;
  typedef T type;
};

template <class T>
struct is_smart_ptr<T, typename std::enable_if<std::is_same<std::weak_ptr<typename T::element_type>, T>::value>::type>
{
  static const bool value = true;
  typedef T type;
};


template <class T, class = void>
struct PtrCaller
{
  static T& call(T& ptr)
  {
    return ptr;
  }
};

template <class T>
struct PtrCaller<T, typename std::enable_if<is_smart_ptr<T>::value || std::is_pointer<T>::value>::type>
{
  static typename T::element_type& call(T& ptr)
  {
    return *ptr;
  }
};


//! workaround for gcc 4.7 missing underlying type, via
//! https://stackoverflow.com/questions/9343329/how-to-know-underlying-type-of-class-enum/10956467#10956467
template <class T>
struct underlying_type
{
#if __GNUC__ == 4 && (__GNUC_MINOR__ < 7)
  typedef typename std::
      conditional<T(-1) < T(0), typename std::make_signed<T>::type, typename std::make_unsigned<T>::type>::type type;
#else
  typedef typename std::underlying_type<T>::type type;
#endif
};


//! gcc < 4.8 fires a static-assert if std::hash< T > () isn't implemented
#if __GNUC__ == 4 && (__GNUC_MINOR__ < 8)
template <typename>
struct is_hashable : std::false_type
{};
#else
//! implementation from https://gcc.gnu.org/ml/libstdc++/2013-03/msg00027.html
template <typename, typename = void>
struct is_hashable : std::false_type
{};

template <typename T>
struct is_hashable<T, typename std::enable_if<!!sizeof(std::declval<std::hash<T>>()(std::declval<T>()))>::type>
  : std::true_type
{};
#endif


template <class T>
struct is_complex : public std::false_type
{};

template <class T>
struct is_complex<std::complex<T>> : public std::true_type
{};


template <class T>
struct is_arithmetic : public std::is_arithmetic<T>
{};

template <int k>
struct is_arithmetic<Dune::bigunsignedint<k>> : public std::true_type
{};


template <class T>
struct is_field_vector : std::false_type
{};

template <class K, int SIZE>
struct is_field_vector<Dune::FieldVector<K, SIZE>> : std::true_type
{};

template <class K, int SIZE>
struct is_field_vector<Dune::XT::Common::FieldVector<K, SIZE>> : std::true_type
{};


template <class T>
struct is_field_matrix : std::false_type
{};

template <class K, int ROWS, int COLS>
struct is_field_matrix<Dune::FieldMatrix<K, ROWS, COLS>> : std::true_type
{};

template <class K, int ROWS, int COLS>
struct is_field_matrix<Dune::XT::Common::FieldMatrix<K, ROWS, COLS>> : std::true_type
{};


namespace internal {


template <class M>
struct is_matrix_helper
{
  DXTC_has_static_member_initialize_once(is_matrix);
  static const bool is_candidate = DXTC_has_static_member(is_matrix)<MatrixAbstraction<M>>::value;
}; // struct is_matrix_helper


template <class V>
struct is_vector_helper
{
  DXTC_has_static_member_initialize_once(is_vector);
  static const bool is_candidate = DXTC_has_static_member(is_vector)<VectorAbstraction<V>>::value;
}; // struct is_vector_helper


} // namespace internal


template <class T, bool candidate = internal::is_matrix_helper<T>::is_candidate>
struct is_matrix
{
  static const bool value = MatrixAbstraction<T>::is_matrix;
};

template <class T>
struct is_matrix<T, false> : public std::false_type
{};


template <class T, bool candidate = internal::is_vector_helper<T>::is_candidate>
struct is_vector
{
  static const bool value = VectorAbstraction<T>::is_vector;
};

template <class T>
struct is_vector<T, false> : public std::false_type
{};


template <class T>
struct field_traits
{
  static_assert(is_arithmetic<T>::value, "If you think otherwise, specialize this class!");
  using field_type = T;
  using real_type = T;
  using complex_type = std::complex<T>;
};

template <class T>
struct field_traits<const T> : public field_traits<T>
{};

template <class T>
struct field_traits<T&> : public field_traits<T>
{};

template <class T>
struct field_traits<const T&> : public field_traits<T>
{};

template <class T>
struct field_traits<std::complex<T>>
{
  static_assert(is_arithmetic<T>::value, "If you think otherwise, specialize this class!");
  using field_type = std::complex<T>;
  using real_type = T;
  using complex_type = std::complex<T>;
};

template <class T>
using field_t = typename field_traits<T>::field_type;

template <class T>
using real_t = typename field_traits<T>::real_type;

template <class T>
using complex_t = typename field_traits<T>::complex_type;

template <class T>
T create_real_or_complex_number(const std::enable_if_t<std::is_arithmetic<T>::value, real_t<T>> real_part,
                                const real_t<T> imag_part)
{
  DUNE_THROW_IF(std::abs(imag_part) > 1e-15,
                Dune::MathError,
                "You are trying to create a real number with non-zero imaginary part!");
  return real_part;
}

template <class T>
T create_real_or_complex_number(const std::enable_if_t<is_complex<T>::value, real_t<T>> real_part,
                                const real_t<T> imag_part)
{
  return T(real_part, imag_part);
}

template <class L, class R>
struct plus_promotion
{
  using type = decltype(std::declval<L>() + std::declval<R>());
};

template <class L, class R>
struct multiplication_promotion
{
  using type = decltype(std::declval<L>() * std::declval<R>());
};

template <class L, class R>
using multiplication_promoted_t = typename multiplication_promotion<L, R>::type;


/**
 * To be used e.g. with AlwaysFalse:
\code
enum class Switch
{
case_1,
case_2,
};

template <Switch important_switch>
struct do_something
{
  static_assert(AlwaysFalse<typename dependent<Switch>::_typename<important_switch>::type>::value,
                "Not Implemented for this case!");
};

template <>
struct do_something<Switch::case_1>
{
  // ...
};
\endcode
 */
template <typename T>
struct dependent
{
  template <T S>
  struct _typename
  {
    typedef void type;
  };
};


template <class K, bool is_number = (is_arithmetic<K>::value || is_complex<K>::value)>
struct suitable_default
{
  static K value()
  {
    return 0;
  }
};

template <class K>
struct suitable_default<K, false>
{
  static K value()
  {
    return K();
  }
};


//! enums for matrix algorithms
enum class StorageLayout
{
  dense_row_major,
  dense_column_major,
  csr,
  csc,
  other
};

enum class MatrixPattern
{
  dense,
  upper_triangular,
  lower_triangular,
  tridiagonal,
  diagonal
};

enum class Transpose
{
  no,
  yes
};


} // namespace Common
} // namespace XT
} // namespace Dune


DUNE_XT_COMMON_TYPENAME(int)
DUNE_XT_COMMON_TYPENAME(double)
DUNE_XT_COMMON_TYPENAME(float)
DUNE_XT_COMMON_TYPENAME(long)
DUNE_XT_COMMON_TYPENAME(unsigned int)
DUNE_XT_COMMON_TYPENAME(unsigned long)
DUNE_XT_COMMON_TYPENAME(char)

#endif // DUNE_XT_COMMON_TYPE_TRAITS_HH
