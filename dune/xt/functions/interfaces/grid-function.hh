// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2019)
//   René Fritze     (2018 - 2019)
//   Tim Keil        (2018)
//   Tobias Leibner  (2017, 2019 - 2020)

#ifndef DUNE_XT_FUNCTIONS_INTERFACES_GRID_FUNCTION_HH
#define DUNE_XT_FUNCTIONS_INTERFACES_GRID_FUNCTION_HH

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
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/grid/view/from-part.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/xt/functions/exceptions.hh>
#include <dune/xt/functions/type_traits.hh>
#include <dune/xt/functions/base/visualization.hh>

#include "element-functions.hh"

namespace Dune {
namespace XT {
namespace Functions {


template <class MinuendType, class SubtrahendType>
class DifferenceGridFunction;

template <class LeftSummandType, class RightSummandType>
class SumGridFunction;

template <class LeftSummandType, class RightSummandType>
class ProductGridFunction;

// template <class Function>
// class DivergenceFunction;


/**
 * \brief Interface for functions which can be localized to an element.
 *
 *        We follow one strong assumption in dune-xt-functions: after being localized to an element, a local function
 *        (or a set of local functions) has to provide the order which is required to integrate it _exactly_ on this
 *        element, and it _should_ thus coincide with a polynomial of respective order on an element. While there exists
 *        non-polynomial data functions, we only have quadratures of finite order and no way to detect
 *        discontinuities of functions (thus the order of the local polynomial approximation must be provided). This
 *        encourages the user to really think about the data functions when implementing/using them, but avoids
 *        situations which could not be handled generically later on.
 */
template <class Element, size_t rangeDim = 1, size_t rangeDimCols = 1, class RangeField = double>
class GridFunctionInterface : public Common::ParametricInterface
{
  using ThisType = GridFunctionInterface<Element, rangeDim, rangeDimCols, RangeField>;

public:
  using LocalFunctionType = ElementFunctionInterface<Element, rangeDim, rangeDimCols, RangeField>;

  using ElementType = typename LocalFunctionType::ElementType;
  using DomainFieldType = typename LocalFunctionType::DomainFieldType;
  static constexpr size_t domain_dim = LocalFunctionType::domain_dim;
  using RangeFieldType = typename LocalFunctionType::RangeFieldType;
  static constexpr size_t range_dim = LocalFunctionType::range_dim;
  static constexpr size_t range_dim_cols = LocalFunctionType::range_dim_cols;

  using E = typename LocalFunctionType::E;
  using D = typename LocalFunctionType::D;
  static constexpr size_t d = LocalFunctionType::d;
  using R = typename LocalFunctionType::R;
  static constexpr size_t r = LocalFunctionType::r;
  static constexpr size_t rC = LocalFunctionType::rC;

  static constexpr bool available = false;

  using DifferenceType = Functions::DifferenceGridFunction<ThisType, ThisType>;
  using SumType = Functions::SumGridFunction<ThisType, ThisType>;

  GridFunctionInterface(const Common::ParameterType& param_type = {})
    : Common::ParametricInterface(param_type)
  {}

  virtual ~GridFunctionInterface() = default;

  static std::string static_id()
  {
    return "dune.xt.functions.gridfunction";
  }

  /**
   * \name ´´This method has to be implemented.''
   * \{
   **/

  virtual std::unique_ptr<LocalFunctionType> local_function() const = 0;

  /**
   * \}
   * \name ´´These methods should be implemented in order to identify the function.''
   * \{
   */

  virtual std::string name() const
  {
    return "dune.xt.functions.gridfunction";
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
  typename std::enable_if<is_grid_function<OtherType>::value, Functions::ProductGridFunction<ThisType, OtherType>>::type
  operator*(const OtherType& other) const
  {
    return Functions::ProductGridFunction<ThisType, OtherType>(*this, other);
  }

  /**
   * \note  We use the SubsamplingVTKWriter (which is better for higher orders) by default: the grid you see in the
   *        visualization may thus be a refinement of the actual grid!
   */
  template <class GridViewType>
  void visualize(const GridViewType& grid_view,
                 const std::string path,
                 const bool subsampling = true,
                 const VTK::OutputType vtk_output_type = VTK::appendedraw,
                 const XT::Common::Parameter& param = {},
                 const VisualizerInterface<r, rC, R>& visualizer = default_visualizer<r, rC, R>()) const
  {
    static_assert(Grid::is_view<GridViewType>::value);
    auto vtk_writer = create_vtkwriter(grid_view, subsampling);
    add_to_vtkwriter(*vtk_writer, param, visualizer);
    write_visualization(*vtk_writer, path, vtk_output_type);
  } // ... visualize(...)

  /**
   * \note  We use the SubsamplingVTKWriter (which is better for higher orders) by default: the grid you see in the
   *        visualization may thus be a refinement of the actual grid!
   * \note  Not yet implemented for vector-valued functions.
   */
  template <class GridViewType>
  void visualize_gradient(const GridViewType& grid_view,
                          const std::string path,
                          const bool subsampling = true,
                          const VTK::OutputType vtk_output_type = VTK::appendedraw,
                          const XT::Common::Parameter& param = {},
                          const VisualizerInterface<d, 1, R>& visualizer = default_visualizer<d, 1, R>()) const
  {
    static_assert(Grid::is_view<GridViewType>::value);
    auto vtk_writer = create_vtkwriter(grid_view, subsampling);
    add_gradient_to_vtkwriter(*vtk_writer, param, visualizer);
    write_visualization(*vtk_writer, path, vtk_output_type);
  } // ... visualize_gradient(...)

  template <class GridViewType>
  std::unique_ptr<VTKWriter<GridViewType>> create_vtkwriter(const GridViewType& grid_view,
                                                            const bool subsampling = true) const
  {
    static_assert(Grid::is_view<GridViewType>::value);
    return subsampling ? std::make_unique<SubsamplingVTKWriter<GridViewType>>(grid_view, /*subsampling_level=*/2)
                       : std::make_unique<VTKWriter<GridViewType>>(grid_view, VTK::nonconforming);
  }

  template <class GridViewType>
  void add_to_vtkwriter(VTKWriter<GridViewType>& vtk_writer,
                        const XT::Common::Parameter& param = {},
                        const VisualizerInterface<r, rC, R>& visualizer = default_visualizer<r, rC, R>()) const
  {
    static_assert(Grid::is_view<GridViewType>::value);
    const auto adapter =
        std::make_shared<VisualizationAdapter<GridViewType, range_dim, range_dim_cols, RangeFieldType>>(
            *this, visualizer, "", param);
    vtk_writer.addVertexData(adapter);
  }

  template <class GridViewType>
  void add_gradient_to_vtkwriter(VTKWriter<GridViewType>& vtk_writer,
                                 const XT::Common::Parameter& param = {},
                                 const VisualizerInterface<d, 1, R>& visualizer = default_visualizer<d, 1, R>()) const
  {
    static_assert(Grid::is_view<GridViewType>::value);
    const auto adapter =
        std::make_shared<GradientVisualizationAdapter<GridViewType, range_dim, range_dim_cols, RangeFieldType>>(
            *this, visualizer, "", param);
    vtk_writer.addCellData(adapter);
  }

  template <class GridViewType>
  auto write_visualization(VTKWriter<GridViewType>& vtk_writer,
                           const std::string path,
                           const VTK::OutputType vtk_output_type = VTK::appendedraw) const
  {
    static_assert(Grid::is_view<GridViewType>::value);
    if (path.empty())
      DUNE_THROW(Exceptions::wrong_input_given, "path must not be empty!");
    const auto directory = Common::directory_only(path);
    Common::test_create_directory(directory);
    if (MPIHelper::getCollectiveCommunication().size() == 1)
      vtk_writer.write(path, vtk_output_type);
    else
      vtk_writer.pwrite(Common::filename_only(path), directory, "", vtk_output_type);
  }
}; // class GridFunctionInterface


} // namespace Functions
} // namespace XT
} // namespace Dune

#include <dune/xt/functions/base/combined-grid-functions.hh>
#include <dune/xt/functions/base/visualization.hh>

#endif // DUNE_XT_FUNCTIONS_INTERFACES_GRID_FUNCTION_HH
