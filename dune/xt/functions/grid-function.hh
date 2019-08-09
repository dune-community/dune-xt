// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_XT_FUNCTIONS_GRID_FUNCTION_HH
#define DUNE_XT_FUNCTIONS_GRID_FUNCTION_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/la/container/eye-matrix.hh>
#include <dune/xt/functions/base/function-as-grid-function.hh>
#include <dune/xt/functions/base/combined-grid-functions.hh>
#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/interfaces/function.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>

namespace Dune {
namespace XT {
namespace Functions {


template <class E, size_t r = 1, size_t rC = 1, class R = double>
class GridFunction : public GridFunctionInterface<E, r, rC, R>
{
  using BaseType = GridFunctionInterface<E, r, rC, R>;
  using ThisType = GridFunction<E, r, rC, R>;

public:
  using BaseType::d;
  using typename BaseType::LocalFunctionType;

  GridFunction(const typename FunctionInterface<d, r, rC, R>::RangeReturnType& value)
    : storage_(new FunctionAsGridFunctionWrapper<E, r, rC, R>(new ConstantFunction<d, r, rC, R>(value)))
  {}

  GridFunction(const FunctionInterface<d, r, rC, R>& func)
    : storage_(new FunctionAsGridFunctionWrapper<E, r, rC, R>(func))
  {}

  GridFunction(FunctionInterface<d, r, rC, R>*&& func_ptr)
    : storage_(new FunctionAsGridFunctionWrapper<E, r, rC, R>(std::move(func_ptr)))
  {}

  GridFunction(const BaseType& func)
    : storage_(func)
  {}

  GridFunction(BaseType*&& func_ptr)
    : storage_(std::move(func_ptr))
  {}

  GridFunction(const ThisType& other)
    : storage_(other.storage_)
  {}

  GridFunction(ThisType&& source)
    : storage_(std::move(source.storage_))
  {}

  std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    return storage_.access().local_function();
  }

  std::string name() const override final
  {
    return storage_.access().name();
  }

private:
  Common::ConstStorageProvider<GridFunctionInterface<E, r, rC, R>> storage_;
}; // class GridFunction<..., r, rC, ...>


template <class E, size_t r, class R>
class GridFunction<E, r, r, R> : public GridFunctionInterface<E, r, r, R>
{
  using BaseType = GridFunctionInterface<E, r, r, R>;
  using ThisType = GridFunction<E, r, r, R>;

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

  GridFunction(const R& value)
    : storage_(make_product(new FunctionAsGridFunctionWrapper<E, 1, 1, R>(new ConstantFunction<d, 1, 1, R>(value)),
                            unit_matrix()))
  {}

  GridFunction(const typename FunctionInterface<d, r, r, R>::RangeReturnType& value)
    : storage_(new FunctionAsGridFunctionWrapper<E, r, r, R>(new ConstantFunction<d, r, r, R>(value)))
  {}

  GridFunction(const FunctionInterface<d, 1, 1, R>& func)
    : storage_(make_product(new FunctionAsGridFunctionWrapper<E, 1, 1, R>(func), std::move(unit_matrix())))
  {}

  GridFunction(FunctionInterface<d, 1, 1, R>*&& func_ptr)
    : storage_(
          make_product(new FunctionAsGridFunctionWrapper<E, 1, 1, R>(std::move(func_ptr)), std::move(unit_matrix())))
  {}

  GridFunction(const FunctionInterface<d, r, r, R>& func)
    : storage_(new FunctionAsGridFunctionWrapper<E, r, r, R>(func))
  {}

  GridFunction(FunctionInterface<d, r, r, R>*&& func_ptr)
    : storage_(new FunctionAsGridFunctionWrapper<E, r, r, R>(std::move(func_ptr)))
  {}

  GridFunction(const GridFunctionInterface<E, 1, 1, R>& func)
    : storage_(make_product(new GridFunction<E, 1, 1, R>(func), std::move(unit_matrix())))
  {}

  GridFunction(GridFunctionInterface<E, 1, 1, R>*&& func_ptr)
    : storage_(make_product(std::move(func_ptr), std::move(unit_matrix())))
  {}

  GridFunction(const GridFunctionInterface<E, r, r, R>& func)
    : storage_(func)
  {}

  GridFunction(GridFunctionInterface<E, r, r, R>*&& func_ptr)
    : storage_(std::move(func_ptr))
  {}

  GridFunction(const ThisType& other)
    : storage_(other.storage_)
  {}

  GridFunction(ThisType&& source)
    : storage_(std::move(source.storage_))
  {}

  std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    return storage_.access().local_function();
  }

  std::string name() const override final
  {
    return storage_.access().name();
  }

private:
  Common::ConstStorageProvider<GridFunctionInterface<E, r, rC, R>> storage_;
}; // class GridFunction<..., r, r, ...>


template <class E, class R>
class GridFunction<E, 1, 1, R> : public GridFunctionInterface<E, 1, 1, R>
{
  using BaseType = GridFunctionInterface<E, 1, 1, R>;
  using ThisType = GridFunction<E, 1, 1, R>;

public:
  using BaseType::d;
  using BaseType::r;
  using BaseType::rC;
  using typename BaseType::LocalFunctionType;

  GridFunction(const typename FunctionInterface<d, r, rC, R>::RangeReturnType& value)
    : storage_(new FunctionAsGridFunctionWrapper<E, r, rC, R>(new ConstantFunction<d, r, rC, R>(value)))
  {}

  GridFunction(const FunctionInterface<d, r, rC, R>& func)
    : storage_(new FunctionAsGridFunctionWrapper<E, r, rC, R>(func))
  {}

  GridFunction(FunctionInterface<d, r, rC, R>*&& func_ptr)
    : storage_(new FunctionAsGridFunctionWrapper<E, r, rC, R>(std::move(func_ptr)))
  {}

  GridFunction(const GridFunctionInterface<E, r, rC, R>& func)
    : storage_(func)
  {}

  GridFunction(GridFunctionInterface<E, r, rC, R>*&& func_ptr)
    : storage_(std::move(func_ptr))
  {}

  GridFunction(const ThisType& other)
    : storage_(other.storage_)
  {}

  GridFunction(ThisType&& source)
    : storage_(std::move(source.storage_))
  {}

  std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    return storage_.access().local_function();
  }

  std::string name() const override final
  {
    return storage_.access().name();
  }

private:
  Common::ConstStorageProvider<GridFunctionInterface<E, r, rC, R>> storage_;
}; // class GridFunction<..., 1, 1, ...>


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_GRID_FUNCTION_HH
