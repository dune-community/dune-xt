// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)
//   René Fritze     (2019)
//   Tobias Leibner  (2019 - 2020)

#ifndef DUNE_XT_FUNCTIONS_GRID_FUNCTION_HH
#define DUNE_XT_FUNCTIONS_GRID_FUNCTION_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/la/container/eye-matrix.hh>
#include <dune/xt/functions/base/function-as-grid-function.hh>
#include <dune/xt/functions/base/combined-grid-functions.hh>
#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/interfaces/function.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>
#include <dune/xt/functions/generic/function.hh>
#include <dune/xt/functions/generic/grid-function.hh>

namespace Dune::XT::Functions {


/**
 * \brief Wraps a value, a function or a grid function (generic variant, except for scalars and square matrices).
 *
 * Suppose we have a function (this is all pseudo-code below, assuming d is the dimension of the grid element E)
 *
\code
void foo(const GridFunctionInterface<E>& grid_func);
\endcode
 *
 * or a class
 *
\code
struc Bar {
  Bar(const GridFunctionInterface<E>& grid_func)
   : grid_func_(grid_func) {}

  const GridFunctionInterface<E>& grid_func_;
};
\endcode
 *
 * that require a grid function to work, but we want the user to be able to also pass a function or a value, for
 * convenience. In case of the function foo, we would thus need to have two additional ones
 *
\code
void foo(const FunctionInterface<d>& func
{
  foo(FunctionAsGridFunctionWrapper<E>(func));
}

void foo(const double& value)
{
  foo(ConstantFunction<d>(value));
}
\endcode
 *
 * and in case of the class Bar we need additional ctors and a different member, as in
 *
\code
struc Bar {
  Bar(const double& grid_func)
   : grid_func_(new FunctionAsGridFunctionWrapper<E>(new ConstantFunction(value))) {}

  Bar(const FunctionInterface<d>& func)
   : grid_func_(new FunctionAsGridFunctionWrapper<E>(func)) {}

  Bar(const GridFunctionInterface<E>& grid_func)
   : grid_func_(grid_func) {}

  Dune::XT::Common::StorageProvider<GridFunctionInterface<E>> grid_func_;
};
\endcode
 *
 * Similar to the last example, this class simplifies such use cases by providing a grid function, given a grid
 * function, a function or a value. If the range of the grid function is a square matrix, it also accepts a number, a
 * scalar function or scalar grid function. To benefit from this class, one needs to use it as an argument passed by
 * value (its copy and movable), so the above examples simplify to
 *
\code
void foo(GridFunction<E> grid_func);
\endcode
 *
 * and
 *
\code
struc Bar {
  Bar(GridFunction<E> grid_func)
   : grid_func_(grid_func) {}

  GridFunctionInterface<E> grid_func_;
};
\endcode
 *
 * This allows the user to call
 *
\code
foo(1);
foo(some_function);
foo(some_grid_function)
\endcode
 *
 * as well as
 *
\code
Bar bar(1);
Bar bar(some_function);
Bar bar(some_grid_function);
\endcode
 *
 * the magic happens due to the various non explicit ctors.
 */
template <class E, size_t r = 1, size_t rC = 1, class R = double>
class GridFunction : public GridFunctionInterface<E, r, rC, R>
{
  using ThisType = GridFunction;
  using BaseType = GridFunctionInterface<E, r, rC, R>;

public:
  using BaseType::d;
  using typename BaseType::LocalFunctionType;
  using GenericFunctionType = GenericFunction<d, r, rC>;

  GridFunction(const typename RangeTypeSelector<R, r, rC>::type& value)
    : BaseType()
    , storage_(new FunctionAsGridFunctionWrapper<E, r, rC, R>(new ConstantFunction<d, r, rC, R>(value)))
  {}

  GridFunction(const FunctionInterface<d, r, rC, R>& func)
    : BaseType()
    , storage_(new FunctionAsGridFunctionWrapper<E, r, rC, R>(func))
  {}

  GridFunction(FunctionInterface<d, r, rC, R>*&& func_ptr)
    : BaseType()
    , storage_(new FunctionAsGridFunctionWrapper<E, r, rC, R>(std::move(func_ptr)))
  {}

  GridFunction(const GridFunctionInterface<E, r, rC, R>& func)
    : BaseType()
    , storage_(func)
  {}

  GridFunction(GridFunctionInterface<E, r, rC, R>*&& func_ptr)
    : BaseType()
    , storage_(std::move(func_ptr))
  {}

  GridFunction(std::tuple<int, typename GenericFunctionType::GenericEvaluateFunctionType> order_evaluate)
    : BaseType()
    , storage_(new FunctionAsGridFunctionWrapper<E, r, rC, R>(
          new GenericFunctionType(std::get<0>(order_evaluate), std::get<1>(order_evaluate))))
  {}

  GridFunction(std::tuple<int,
                          typename GenericFunctionType::GenericEvaluateFunctionType,
                          typename GenericFunctionType::GenericJacobianFunctionType> order_evaluate_jacobian)
    : BaseType()
    , storage_(
          new FunctionAsGridFunctionWrapper<E, r, rC, R>(new GenericFunctionType(std::get<0>(order_evaluate_jacobian),
                                                                                 std::get<1>(order_evaluate_jacobian),
                                                                                 /*name=*/"",
                                                                                 /*param_type=*/{},
                                                                                 std::get<2>(order_evaluate_jacobian))))
  {}

  GridFunction(const ThisType& other)
    : BaseType(other)
    , storage_(other.storage_)
  {}

  GridFunction(ThisType&& source)
    : BaseType(source)
    , storage_(std::move(source.storage_))
  {}

  std::unique_ptr<LocalFunctionType> local_function() const final
  {
    return storage_.access().local_function();
  }

  std::string name() const final
  {
    return storage_.access().name();
  }

private:
  Common::ConstStorageProvider<GridFunctionInterface<E, r, rC, R>> storage_;
}; // class GridFunction<..., r, rC, ...>


/**
 * \brief Wraps a value, a function or a grid function (variant for square matrices).
 *
 * \sa GridFunction
 */
template <class E, size_t r, class R>
class GridFunction<E, r, r, R> : public GridFunctionInterface<E, r, r, R>
{
  using BaseType = GridFunctionInterface<E, r, r, R>;
  using ThisType = GridFunction;

private:
  static GridFunctionInterface<E, r, r, R>* unit_matrix()
  {
    return new FunctionAsGridFunctionWrapper<E, r, r, R>(new ConstantFunction<d, r, r, R>(
        XT::LA::eye_matrix<typename FunctionInterface<d, r, r, R>::RangeReturnType>(r, r)));
  }

public:
  using BaseType::d;
  using BaseType::rC;
  using typename BaseType::LocalFunctionType;
  using GenericFunctionType = GenericFunction<d, r, rC>;

  GridFunction(const R& value)
    : BaseType()
    , storage_(new ProductGridFunction<GridFunction<E, 1, 1, R>, GridFunctionInterface<E, r, r, R>>(
          new GridFunction<E, 1, 1, R>(value), std::move(unit_matrix()), ""))
  {}

  GridFunction(const FieldMatrix<R, r, r>& value) // <- Must not be XT::Common::FieldMatrix!
    : BaseType()
    , storage_(new FunctionAsGridFunctionWrapper<E, r, r, R>(new ConstantFunction<d, r, r, R>(value)))
  {}

  GridFunction(const FunctionInterface<d, 1, 1, R>& func)
    : BaseType()
    , storage_(new ProductGridFunction<FunctionAsGridFunctionWrapper<E, 1, 1, R>, GridFunctionInterface<E, r, r, R>>(
          new FunctionAsGridFunctionWrapper<E, 1, 1, R>(func), std::move(unit_matrix()), func.name()))
  {}

  GridFunction(FunctionInterface<d, 1, 1, R>*&& func_ptr)
    : BaseType()
    , storage_(new ProductGridFunction<FunctionAsGridFunctionWrapper<E, 1, 1, R>, GridFunctionInterface<E, r, r, R>>(
          new FunctionAsGridFunctionWrapper<E, 1, 1, R>(std::move(func_ptr)),
          std::move(unit_matrix()),
          func_ptr->name()))
  {}

  GridFunction(const FunctionInterface<d, r, r, R>& func)
    : BaseType()
    , storage_(new FunctionAsGridFunctionWrapper<E, r, r, R>(func))
  {}

  GridFunction(FunctionInterface<d, r, r, R>*&& func_ptr)
    : BaseType()
    , storage_(new FunctionAsGridFunctionWrapper<E, r, r, R>(std::move(func_ptr)))
  {}

  GridFunction(const GridFunctionInterface<E, 1, 1, R>& func)
    : BaseType()
    , storage_(new ProductGridFunction<GridFunction<E, 1, 1, R>, GridFunctionInterface<E, r, r, R>>(
          new GridFunction<E, 1, 1, R>(func), std::move(unit_matrix()), func.name()))
  {}

  GridFunction(GridFunctionInterface<E, 1, 1, R>*&& func_ptr)
    : BaseType()
    , storage_(new ProductGridFunction<GridFunctionInterface<E, 1, 1, R>, GridFunctionInterface<E, r, r, R>>(
          std::move(func_ptr), std::move(unit_matrix()), func_ptr->name()))
  {}

  GridFunction(const GridFunctionInterface<E, r, r, R>& func)
    : BaseType()
    , storage_(func)
  {}

  GridFunction(GridFunctionInterface<E, r, r, R>*&& func_ptr)
    : BaseType()
    , storage_(std::move(func_ptr))
  {}

  GridFunction(std::tuple<int, typename GenericFunctionType::GenericEvaluateFunctionType> order_evaluate)
    : BaseType()
    , storage_(new FunctionAsGridFunctionWrapper<E, r, rC, R>(
          new GenericFunctionType(std::get<0>(order_evaluate), std::get<1>(order_evaluate))))
  {}

  GridFunction(std::tuple<int,
                          typename GenericFunctionType::GenericEvaluateFunctionType,
                          typename GenericFunctionType::GenericJacobianFunctionType> order_evaluate_jacobian)
    : BaseType()
    , storage_(
          new FunctionAsGridFunctionWrapper<E, r, rC, R>(new GenericFunctionType(std::get<0>(order_evaluate_jacobian),
                                                                                 std::get<1>(order_evaluate_jacobian),
                                                                                 /*name=*/"",
                                                                                 /*param_type=*/{},
                                                                                 std::get<2>(order_evaluate_jacobian))))
  {}

  GridFunction(const ThisType& other)
    : BaseType(other)
    , storage_(other.storage_)
  {}

  GridFunction(ThisType&& source)
    : BaseType(source)
    , storage_(std::move(source.storage_))
  {}

  std::unique_ptr<LocalFunctionType> local_function() const final
  {
    return storage_.access().local_function();
  }

  std::string name() const final
  {
    return storage_.access().name();
  }

private:
  Common::ConstStorageProvider<GridFunctionInterface<E, r, rC, R>> storage_;
}; // class GridFunction<..., r, r, ...>


/**
 * \brief Wraps a value, a function or a grid function (scalar variant).
 *
 * \sa GridFunction
 */
template <class E, class R>
class GridFunction<E, 1, 1, R> : public GridFunctionInterface<E, 1, 1, R>
{
  using ThisType = GridFunction;
  using BaseType = GridFunctionInterface<E, 1, 1, R>;

public:
  using BaseType::d;
  using BaseType::r;
  using BaseType::rC;
  using typename BaseType::LocalFunctionType;
  using GenericFunctionType = GenericFunction<d, r, rC>;

  GridFunction(const R& value)
    : BaseType()
    , storage_(new FunctionAsGridFunctionWrapper<E, 1, 1, R>(new ConstantFunction<d, 1, 1, R>(value)))
  {}

  GridFunction(const FieldVector<R, 1>& value) // <- Must not be XT::Common::FieldVector!
    : BaseType()
    , storage_(new FunctionAsGridFunctionWrapper<E, 1, 1, R>(new ConstantFunction<d, 1, 1, R>(value)))
  {}

  GridFunction(const FieldMatrix<R, 1, 1>& value) // <- Must not be XT::Common::FieldMatrix!
    : BaseType()
    , storage_(new FunctionAsGridFunctionWrapper<E, 1, 1, R>(new ConstantFunction<d, 1, 1, R>(value[0][0])))
  {}

  GridFunction(const FunctionInterface<d, 1, 1, R>& func)
    : BaseType()
    , storage_(new FunctionAsGridFunctionWrapper<E, 1, 1, R>(func))
  {}

  GridFunction(FunctionInterface<d, 1, 1, R>*&& func_ptr)
    : BaseType()
    , storage_(new FunctionAsGridFunctionWrapper<E, 1, 1, R>(std::move(func_ptr)))
  {}

  GridFunction(const GridFunctionInterface<E, 1, 1, R>& func)
    : BaseType()
    , storage_(func)
  {}

  GridFunction(GridFunctionInterface<E, 1, 1, R>*&& func_ptr)
    : BaseType()
    , storage_(std::move(func_ptr))
  {}

  GridFunction(std::tuple<int, typename GenericFunctionType::GenericEvaluateFunctionType> order_evaluate)
    : BaseType()
    , storage_(new FunctionAsGridFunctionWrapper<E, r, rC, R>(
          new GenericFunctionType(std::get<0>(order_evaluate), std::get<1>(order_evaluate))))
  {}

  GridFunction(std::tuple<int,
                          typename GenericFunctionType::GenericEvaluateFunctionType,
                          typename GenericFunctionType::GenericJacobianFunctionType> order_evaluate_jacobian)
    : BaseType()
    , storage_(
          new FunctionAsGridFunctionWrapper<E, r, rC, R>(new GenericFunctionType(std::get<0>(order_evaluate_jacobian),
                                                                                 std::get<1>(order_evaluate_jacobian),
                                                                                 /*name=*/"",
                                                                                 /*param_type=*/{},
                                                                                 std::get<2>(order_evaluate_jacobian))))
  {}

  GridFunction(const ThisType& other)
    : BaseType(other)
    , storage_(other.storage_)
  {}

  GridFunction(ThisType&& source)
    : BaseType(source)
    , storage_(std::move(source.storage_))
  {}

  std::unique_ptr<LocalFunctionType> local_function() const final
  {
    return storage_.access().local_function();
  }

  std::string name() const final
  {
    return storage_.access().name();
  }

private:
  Common::ConstStorageProvider<GridFunctionInterface<E, 1, 1, R>> storage_;
}; // class GridFunction<..., 1, 1, ...>


} // namespace Dune::XT::Functions

#endif // DUNE_XT_FUNCTIONS_GRID_FUNCTION_HH
