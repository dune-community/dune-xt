// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)
//   René Fritze     (2020)
//   Tobias Leibner  (2020)

#ifndef DUNE_XT_FUNCTIONS_VISUALIZATION_HH
#define DUNE_XT_FUNCTIONS_VISUALIZATION_HH

#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/base/visualization.hh>
#include <dune/xt/functions/interfaces/function.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>
#include <dune/xt/functions/grid-function.hh>

namespace Dune::XT::Functions {
namespace internal {


template <class GridViewType>
std::unique_ptr<VTKWriter<GridViewType>>
create_vtkwriter(const GridViewType& grid_view, const bool subsampling = true, const int subsampling_level = 2)
{
  static_assert(Grid::is_view<GridViewType>::value);
  return subsampling ? std::make_unique<SubsamplingVTKWriter<GridViewType>>(grid_view, subsampling_level)
                     : std::make_unique<VTKWriter<GridViewType>>(grid_view, VTK::nonconforming);
}


template <class GridViewType, class E, size_t r, size_t rC, class R>
void add_to_vtkwriter(VTKWriter<GridViewType>& vtk_writer,
                      const GridFunctionInterface<E, r, rC, R>& grid_function,
                      const XT::Common::Parameter& param = {},
                      const VisualizerInterface<r, rC, R>& visualizer = default_visualizer<r, rC, R>())
{
  static_assert(Grid::is_view<GridViewType>::value);
  static_assert(std::is_same<E, XT::Grid::extract_entity_t<GridViewType>>::value);
  vtk_writer.addVertexData(
      std::make_shared<VisualizationAdapter<GridViewType, r, rC, R>>(grid_function, visualizer, "", param));
}

template <class GridViewType, size_t d, size_t r, size_t rC, class R>
void add_to_vtkwriter(VTKWriter<GridViewType>& vtk_writer,
                      const FunctionInterface<d, r, rC, R>& function,
                      const XT::Common::Parameter& param = {},
                      const VisualizerInterface<r, rC, R>& visualizer = default_visualizer<r, rC, R>())
{
  static_assert(Grid::is_view<GridViewType>::value);
  static_assert(GridViewType::dimension == d);
  add_to_vtkwriter(vtk_writer, make_grid_function<GridViewType>(function), param, visualizer);
}


template <class GridViewType, class E, size_t r, size_t rC, class R>
void add_gradient_to_vtkwriter(VTKWriter<GridViewType>& vtk_writer,
                               const GridFunctionInterface<E, r, rC, R>& grid_function,
                               const XT::Common::Parameter& param = {},
                               const VisualizerInterface<GridViewType::dimension, 1, R>& visualizer =
                                   default_visualizer<GridViewType::dimension, 1, R>())
{
  static_assert(Grid::is_view<GridViewType>::value);
  static_assert(std::is_same<E, XT::Grid::extract_entity_t<GridViewType>>::value);
  const auto adapter =
      std::make_shared<GradientVisualizationAdapter<GridViewType, r, rC, R>>(grid_function, visualizer, "", param);
  vtk_writer.addCellData(adapter);
}

template <class GridViewType, size_t d, size_t r, size_t rC, class R>
void add_gradient_to_vtkwriter(VTKWriter<GridViewType>& vtk_writer,
                               const FunctionInterface<d, r, rC, R>& function,
                               const XT::Common::Parameter& param = {},
                               const VisualizerInterface<GridViewType::dimension, 1, R>& visualizer =
                                   default_visualizer<GridViewType::dimension, 1, R>())
{
  static_assert(Grid::is_view<GridViewType>::value);
  static_assert(GridViewType::dimension == d);
  add_gradient_to_vtkwriter(vtk_writer, make_grid_function<GridViewType>(function), param, visualizer);
}


template <class GridViewType>
auto write_visualization(VTKWriter<GridViewType>& vtk_writer,
                         const std::string path,
                         const VTK::OutputType vtk_output_type = VTK::appendedraw)
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


} // namespace internal


/**
 * \note  We use the SubsamplingVTKWriter (which is better for higher orders) by default: the grid you see in the
 *        visualization may thus be a refinement of the actual grid!
 */
template <class GridViewType, size_t d, size_t r, size_t rC, class R>
void visualize(const FunctionInterface<d, r, rC, R>& function,
               const GridViewType& grid_view,
               const std::string path,
               const bool subsampling = true,
               const VTK::OutputType vtk_output_type = VTK::appendedraw,
               const XT::Common::Parameter& param = {},
               const VisualizerInterface<r, rC, R>& visualizer = default_visualizer<r, rC, R>())
{
  static_assert(Grid::is_view<GridViewType>::value);
  static_assert(GridViewType::dimension == d);
  auto vtk_writer = internal::create_vtkwriter(grid_view, subsampling);
  internal::add_to_vtkwriter(*vtk_writer, function, param, visualizer);
  internal::write_visualization(*vtk_writer, path, vtk_output_type);
} // ... visualize(...)


/**
 * \note  We use the SubsamplingVTKWriter (which is better for higher orders) by default: the grid you see in the
 *        visualization may thus be a refinement of the actual grid!
 */
template <class GridViewType, class E, size_t r, size_t rC, class R>
void visualize(const GridFunctionInterface<E, r, rC, R>& grid_function,
               const GridViewType& grid_view,
               const std::string path,
               const bool subsampling = true,
               const VTK::OutputType vtk_output_type = VTK::appendedraw,
               const XT::Common::Parameter& param = {},
               const VisualizerInterface<r, rC, R>& visualizer = default_visualizer<r, rC, R>())
{
  static_assert(Grid::is_view<GridViewType>::value);
  static_assert(std::is_same<E, XT::Grid::extract_entity_t<GridViewType>>::value);
  auto vtk_writer = internal::create_vtkwriter(grid_view, subsampling);
  internal::add_to_vtkwriter(*vtk_writer, grid_function, param, visualizer);
  internal::write_visualization(*vtk_writer, path, vtk_output_type);
} // ... visualize(...)


/**
 * \note  We use the SubsamplingVTKWriter (which is better for higher orders) by default: the grid you see in the
 *        visualization may thus be a refinement of the actual grid!
 * \note  Not yet implemented for vector-valued functions.
 */
template <class GridViewType, size_t d, size_t r, size_t rC, class R>
void visualize_gradient(const FunctionInterface<d, r, rC, R>& function,
                        const GridViewType& grid_view,
                        const std::string path,
                        const bool subsampling = true,
                        const VTK::OutputType vtk_output_type = VTK::appendedraw,
                        const XT::Common::Parameter& param = {},
                        const VisualizerInterface<d, 1, R>& visualizer = default_visualizer<d, 1, R>())
{
  static_assert(Grid::is_view<GridViewType>::value);
  auto vtk_writer = internal::create_vtkwriter(grid_view, subsampling);
  add_gradient_to_vtkwriter(*vtk_writer, function, param, visualizer);
  internal::write_visualization(*vtk_writer, path, vtk_output_type);
} // ... visualize_gradient(...)


/**
 * \note  We use the SubsamplingVTKWriter (which is better for higher orders) by default: the grid you see in the
 *        visualization may thus be a refinement of the actual grid!
 * \note  Not yet implemented for vector-valued functions.
 */
template <class GridViewType, class E, size_t r, size_t rC, class R>
void visualize_gradient(const GridFunctionInterface<E, r, rC, R>& grid_function,
                        const GridViewType& grid_view,
                        const std::string path,
                        const bool subsampling = true,
                        const VTK::OutputType vtk_output_type = VTK::appendedraw,
                        const XT::Common::Parameter& param = {},
                        const VisualizerInterface<GridViewType::dimension, 1, R>& visualizer =
                            default_visualizer<GridViewType::dimension, 1, R>())
{
  static_assert(Grid::is_view<GridViewType>::value);
  auto vtk_writer = internal::create_vtkwriter(grid_view, subsampling);
  internal::add_gradient_to_vtkwriter(*vtk_writer, grid_function, param, visualizer);
  internal::write_visualization(*vtk_writer, path, vtk_output_type);
} // ... visualize_gradient(...)


} // namespace Dune::XT::Functions

#endif // DUNE_XT_FUNCTIONS_VISUALIZATION_HH
