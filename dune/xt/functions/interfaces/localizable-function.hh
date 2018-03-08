// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)
//   Tobias Leibner  (2017)

#ifndef DUNE_XT_FUNCTIONS_INTERFACES_LOCALIZABLE_FUNCTION_HH
#define DUNE_XT_FUNCTIONS_INTERFACES_LOCALIZABLE_FUNCTION_HH

#include <memory>
#include <ostream>
#include <string>
#include <type_traits>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/io/file/vtk.hh>

#include <dune/xt/common/filesystem.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/grid/dd/subdomains/grid.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/grid/view/from-part.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/xt/functions/exceptions.hh>
#include <dune/xt/functions/type_traits.hh>

#include "local-functions.hh"

namespace Dune {
namespace XT {
namespace Functions {


// forward, required in LocalizableFunctionInterface::visualize
template <class GridViewType, size_t dimRange, size_t dimRangeCols, class RangeFieldImp>
class VisualizationAdapter;

template <class MinuendType, class SubtrahendType>
class DifferenceFunction;

template <class LeftSummandType, class RightSummandType>
class SumFunction;

template <class LeftSummandType, class RightSummandType>
class ProductFunction;

// template <class FunctionImp>
// class DivergenceFunction;


/**
 * \brief Interface for functions which can be localized to an element.
 *
 *        We follow one strong assumption in dune-xt-functions: after being localized to an element, a local function
 *        (or a set of local functions) has to provide the order which is required to integrate it _exactly_ on this
 *        element and _should_ be polynomial on an element. While there exists non-polynomial data functions, we only
 * have
 *        quadratures of finite order and no way to detect discontinuities of functions (thus the order of the local
 *        polynomial approximation must be provided). In addition, we demand that discontinuities may only lie on
 *        intersections (which is not uncommen). This forces the user to really think about the data functions when
 *        implementing/using them, but avoids situations which could not be handled generically later on.
 */
template <class ElementImp, size_t rangeDim = 1, size_t rangeDimCols = 1, class RangeFieldImp = double>
class LocalizableFunctionInterface : public Common::ParametricInterface
{
  using ThisType = LocalizableFunctionInterface<ElementImp, rangeDim, rangeDimCols, RangeFieldImp>;

public:
  using LocalFunctionType = LocalFunctionInterface<ElementImp, rangeDim, rangeDimCols, RangeFieldImp>;

  using ElementType = typename LocalFunctionType::ElementType;
  using DomainFieldType = typename LocalFunctionType::DomainFieldType;
  static const constexpr size_t dimDomain = LocalFunctionType::dimDomain;
  using RangeFieldType = typename LocalFunctionType::RangeFieldType;
  static const constexpr size_t dimRange = LocalFunctionType::dimRange;
  static const constexpr size_t dimRangeCols = LocalFunctionType::dimRangeCols;

  using E = typename LocalFunctionType::E;
  using D = typename LocalFunctionType::D;
  static const constexpr size_t d = LocalFunctionType::d;
  using R = typename LocalFunctionType::R;
  static const constexpr size_t r = LocalFunctionType::r;
  static const constexpr size_t rC = LocalFunctionType::rC;

  static const constexpr bool available = false;

  typedef Functions::DifferenceFunction<ThisType, ThisType> DifferenceType;
  typedef Functions::SumFunction<ThisType, ThisType> SumType;

  virtual ~LocalizableFunctionInterface() = default;

  static std::string static_id()
  {
    return "localizable_function";
  }

  /**
   * \name ´´These methods have to be implemented.''
   * \{
   **/

  virtual std::unique_ptr<LocalFunctionType> local_function() const = 0;

  virtual std::unique_ptr<LocalFunctionType> local_function(const ElementType& /*element*/) const = 0;

  /**
   * \}
   * \name ´´These methods should be implemented in order to identify the function.''
   * \{
   */

  virtual std::string type() const
  {
    return "localizable_function";
  }

  virtual std::string name() const
  {
    return "localizable_function";
  }

  /// \}

  DifferenceType operator-(const ThisType& other) const
  {
    return DifferenceType(*this, other);
  }

  SumType operator+(const ThisType& other) const
  {
    return SumType(*this, other);
  }

  template <class OtherType>
  typename std::enable_if<is_localizable_function<OtherType>::value,
                          Functions::ProductFunction<ThisType, OtherType>>::type
  operator*(const OtherType& other) const
  {
    return Functions::ProductFunction<ThisType, OtherType>(*this, other);
  }

  /**
   * \note  We use the SubsamplingVTKWriter (which is better for higher orders) by default: the grid you see in the
   *        visualization may thus be a refinement of the actual grid!
   */
  template <class GridLayerType>
  typename std::enable_if<Grid::is_layer<GridLayerType>::value, void>::type
  visualize(const GridLayerType& grid_layer,
            const std::string path,
            const bool subsampling = true,
            const VTK::OutputType vtk_output_type = VTK::appendedraw,
            const XT::Common::Parameter& param = {}) const
  {
    if (path.empty())
      DUNE_THROW(Exceptions::wrong_input_given, "path must not be empty!");
    const auto directory = Common::directory_only(path);
    Common::test_create_directory(directory);
    const auto tmp_grid_view = Grid::make_tmp_view(grid_layer);
    const auto& grid_view = tmp_grid_view.access();
    using GridViewType = std::decay_t<decltype(grid_view)>;
    const auto adapter =
        std::make_shared<VisualizationAdapter<GridViewType, dimRange, dimRangeCols, RangeFieldType>>(*this, "", param);
    std::unique_ptr<VTKWriter<GridViewType>> vtk_writer =
        subsampling ? Common::make_unique<SubsamplingVTKWriter<GridViewType>>(grid_view, /*subsampling_level=*/2)
                    : Common::make_unique<VTKWriter<GridViewType>>(grid_view, VTK::nonconforming);
    vtk_writer->addVertexData(adapter);
    if (MPIHelper::getCollectiveCommunication().size() == 1)
      vtk_writer->write(path, vtk_output_type);
    else
      vtk_writer->pwrite(Common::filename_only(path), directory, "", vtk_output_type);
  } // ... visualize(...)
}; // class LocalizableFunctionInterface


} // namespace Functions
} // namespace XT
} // namespace Dune

#include "../combined.hh"
#include <dune/xt/functions/base/visualization.hh>
//#include "../derived.hh"

#endif // DUNE_XT_FUNCTIONS_INTERFACES_LOCALIZABLE_FUNCTION_HH
