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
#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/grid/view/from-part.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/xt/functions/base/visualization.hh>
#include <dune/xt/functions/exceptions.hh>
#include <dune/xt/functions/type_traits.hh>

#include "element-functions.hh"

namespace Dune::XT::Functions {


// forwards, includes are below
template <class, class, class CombinationType>
struct CombinedHelper;

template <class, class>
class DifferenceGridFunction;

template <class, class>
class SumGridFunction;

template <class, class>
class ProductGridFunction;

template <class, class>
class FractionGridFunction;


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
class GridFunctionInterface
  : public Common::ParametricInterface
  , public Common::WithLogger<GridFunctionInterface<Element, rangeDim, rangeDimCols, RangeField>>
{
  using ThisType = GridFunctionInterface;
  using Logger = Common::WithLogger<GridFunctionInterface<Element, rangeDim, rangeDimCols, RangeField>>;

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

private:
  std::string logging_id() const
  {
    return "GridFunctionInterface<" + Common::to_string(size_t(r)) + "," + Common::to_string(size_t(rC)) + ">";
  }

public:
  GridFunctionInterface(const Common::ParameterType& param_type = {},
                        const std::string& logging_prefix = "",
                        const std::array<bool, 3>& logging_state = {false, false, true})
    : Common::ParametricInterface(param_type)
    , Logger(logging_prefix.empty() ? "GridFunctionInterface" : logging_prefix, logging_state)
  {
    LOG_(debug) << logging_id() << "(param_type=" << param_type << ")" << std::endl;
  }

  GridFunctionInterface(const ThisType& other)
    : Common::ParametricInterface(other)
    , Logger(other)
  {}

  GridFunctionInterface(ThisType&) = default;

  virtual ~GridFunctionInterface() = default;

  ThisType& operator=(const ThisType&) = delete;

  ThisType& operator=(ThisType&&) = delete;

  static std::string static_id()
  {
    return "dune.xt.functions.gridfunction";
  }

  /**
   * \name ´´These methods have to be implemented.''
   * \{
   **/

  /**
   * \brief Returns a (shallow) copy of the function.
   * actual implementation work is delegated to the private `copy_as_function_impl`
   * combined with hiding `copy_as_function` in dervived classes, this allows us the a
   * unique_ptr with correct type at all levels of the polymorphic hierarchy
   *
   * \note This is intended to be cheap, so make sure to share resources (but in a thread-safe way)!
   */
  std::unique_ptr<ThisType> copy_as_grid_function() const
  {
    return std::unique_ptr<ThisType>(this->copy_as_grid_function_impl());
  }

private:
  virtual ThisType* copy_as_grid_function_impl() const = 0;

public:
  /**
   * \brief Returns the local function which can be bound to grid elements.
   *
   * \note If possible, the returned function should be able to live on its own, e.g. by copying the grid function.
   */
  virtual std::unique_ptr<LocalFunctionType> local_function() const = 0;

  /**
   * \}
   * \name ´´These methods should be implemented in order to identify the function.''
   * \{
   */

  virtual std::string name() const
  {
    LOG_(debug) << logging_id() << "::name()\n   returning \"GridFunction\"" << std::endl;
    return "GridFunction";
  }

  /// \}

  /// \name Operators emulating numeric types.
  /// \{

  Functions::DifferenceGridFunction<ThisType, ThisType> operator-(const ThisType& other) const
  {
    std::string derived_logging_prefix = "";
    if (this->logger.debug_enabled || other.logger.debug_enabled) {
      derived_logging_prefix = "(" + this->logger.prefix + " - " + other.logger.prefix + ")";
      this->logger.debug() << logging_id() << "::operator-(other=" << &other << ")" << std::endl;
    }
    return Functions::DifferenceGridFunction<ThisType, ThisType>(
        *this, other, "(" + this->name() + " - " + other.name() + ")", derived_logging_prefix);
  }

  Functions::SumGridFunction<ThisType, ThisType> operator+(const ThisType& other) const
  {
    std::string derived_logging_prefix = "";
    if (this->logger.debug_enabled || other.logger.debug_enabled) {
      derived_logging_prefix = "(" + this->logger.prefix + " - " + other.logger.prefix + ")";
      this->logger.debug() << logging_id() << "::operator+(other=" << &other << ")" << std::endl;
    }
    return Functions::SumGridFunction<ThisType, ThisType>(
        *this, other, "(" + this->name() + " + " + other.name() + ")", derived_logging_prefix);
  }

  template <class OtherType>
  std::enable_if_t<is_grid_function<OtherType>::value
                       && internal::CombinedHelper<ThisType, OtherType, CombinationType::product>::available,
                   Functions::ProductGridFunction<ThisType, as_grid_function_interface_t<OtherType>>>
  operator*(const OtherType& other) const
  {
    std::string derived_logging_prefix = "";
    if (this->logger.debug_enabled || other.logger.debug_enabled) {
      derived_logging_prefix = "(" + this->logger.prefix + "*" + other.logger.prefix + ")";
      this->logger.debug() << logging_id() << "::operator*(other=" << &other << ")" << std::endl;
    }
    return Functions::ProductGridFunction<ThisType, as_grid_function_interface_t<OtherType>>(
        *this, other, "(" + this->name() + "*" + other.name() + ")", derived_logging_prefix);
  }

  template <class OtherType>
  std::enable_if_t<is_grid_function<OtherType>::value
                       && internal::CombinedHelper<ThisType, OtherType, CombinationType::fraction>::available,
                   Functions::FractionGridFunction<ThisType, as_grid_function_interface_t<OtherType>>>
  operator/(const OtherType& other) const
  {
    std::string derived_logging_prefix = "";
    if (this->logger.debug_enabled || other.logger.debug_enabled) {
      derived_logging_prefix = "(" + this->logger.prefix + "/" + other.logger.prefix + ")";
      this->logger.debug() << logging_id() << "::operator/(other=" << &other << ")" << std::endl;
    }
    return Functions::FractionGridFunction<ThisType, as_grid_function_interface_t<OtherType>>(
        *this, other, "(" + this->name() + "/" + other.name() + ")", derived_logging_prefix);
  }

  /// \}

  /**
   * \note  We use the SubsamplingVTKWriter (which is better for higher orders) by default: the grid you see in the
   *        visualization may thus be a refinement of the actual grid!
   */
  template <class GridViewType>
  [[deprecated("Use visualize(*this, ...) from <dune/xt/functions/visualization.hh> instead (24.09.2020)!")]] void
  visualize(const GridViewType& grid_view,
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
  [[deprecated(
      "Use visualize_gradient(*this, ...) from <dune/xt/functions/visualization.hh> instead (24.09.2020)!")]] void
  visualize_gradient(const GridViewType& grid_view,
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
  [[deprecated(
      "Use internal::create_vtkwriter(...) from <dune/xt/functions/visualization.hh> instead (24.09.2020)!")]] std::
      unique_ptr<VTKWriter<GridViewType>>
      create_vtkwriter(const GridViewType& grid_view, const bool subsampling = true) const
  {
    static_assert(Grid::is_view<GridViewType>::value);
    return subsampling ? std::make_unique<SubsamplingVTKWriter<GridViewType>>(grid_view, /*subsampling_level=*/2)
                       : std::make_unique<VTKWriter<GridViewType>>(grid_view, VTK::nonconforming);
  }

  template <class GridViewType>
  [[deprecated(
      "Use internal::add_to_vtkwriter(...) from <dune/xt/functions/visualization.hh> instead (24.09.2020)!")]] void
  add_to_vtkwriter(VTKWriter<GridViewType>& vtk_writer,
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
  [[deprecated("Use internal::add_gradient_to_vtkwriter(...) from <dune/xt/functions/visualization.hh> instead "
               "(24.09.2020)!")]] void
  add_gradient_to_vtkwriter(VTKWriter<GridViewType>& vtk_writer,
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
  [[deprecated("Use internal::write_visualization(...) from <dune/xt/functions/visualization.hh> instead "
               "(24.09.2020)!")]] void
  write_visualization(VTKWriter<GridViewType>& vtk_writer,
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


} // namespace Dune::XT::Functions

#include <dune/xt/functions/base/combined-grid-functions.hh>

#endif // DUNE_XT_FUNCTIONS_INTERFACES_GRID_FUNCTION_HH
