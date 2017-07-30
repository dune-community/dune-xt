// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
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

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/filesystem.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/grid/dd/subdomains/grid.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/grid/view/from-part.hh>
#include <dune/xt/grid/type_traits.hh>

#include "../type_traits.hh"

#include "local-functions.hh"

namespace Dune {
namespace XT {
namespace Functions {


template <class GridViewType, size_t dimRange, size_t dimRangeCols = 1>
class VisualizationAdapterFunction;


template <class MinuendType, class SubtrahendType>
class DifferenceFunction;

template <class LeftSummandType, class RightSummandType>
class SumFunction;

template <class LeftSummandType, class RightSummandType>
class ProductFunction;

template <class FunctionImp>
class DivergenceFunction;


/**
 * \brief Interface for functions which provide a LocalfunctionInterface for an entity.
 */
template <class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols = 1>
class LocalizableFunctionInterface : public Common::ParametricInterface
{
  typedef LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
      ThisType;

public:
  typedef EntityImp EntityType;

  typedef DomainFieldImp DomainFieldType;
  static const constexpr size_t dimDomain = domainDim;

  typedef RangeFieldImp RangeFieldType;
  static const constexpr size_t dimRange = rangeDim;
  static const constexpr size_t dimRangeCols = rangeDimCols;

  typedef EntityType E;
  typedef DomainFieldType D;
  static const constexpr size_t d = dimDomain;
  typedef RangeFieldType R;
  static const constexpr size_t r = dimRange;
  static const constexpr size_t rC = dimRangeCols;

  typedef LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols>
      LocalfunctionType;

  typedef typename LocalfunctionType::DomainType DomainType;
  typedef typename LocalfunctionType::RangeType RangeType;
  typedef typename LocalfunctionType::JacobianRangeType JacobianRangeType;

  static const bool available = false;

  typedef Functions::DifferenceFunction<ThisType, ThisType> DifferenceType;
  typedef Functions::SumFunction<ThisType, ThisType> SumType;

  typedef void DivergenceType;

  virtual ~LocalizableFunctionInterface() = default;

  static std::string static_id()
  {
    return "function";
  }

  /**
   * \name ´´These methods have to be implemented.''
   * \{
   **/
  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& /*entity*/) const = 0;
  /* \} */

  /**
   * \name ´´These methods should be implemented in order to identify the function.''
   * \{
   */
  virtual std::string type() const
  {
    return "function";
  }

  virtual std::string name() const
  {
    return "function";
  }
  /* \} */

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

  void DUNE_DEPRECATED_MSG("use make_divergence(func) instead of func.divergence() (09.05.2017)!") divergence() const
  {
  }

  /**
   * \note  We use the SubsamplingVTKWriter (which is better for higher orders) by default. The grid you see in the
   *  visualization may thus be a refinement of the actual grid!
   */
  template <class GridLayerType>
  typename std::enable_if<Grid::is_layer<GridLayerType>::value, void>::type
  visualize(const GridLayerType& grid_layer,
            const std::string path,
            const bool subsampling = true,
            const VTK::OutputType vtk_output_type = VTK::appendedraw) const
  {
    if (path.empty())
      DUNE_THROW(RangeError, "Empty path given!");
    const auto directory = Common::directory_only(path);
    const auto tmp_grid_view = Grid::make_tmp_view(grid_layer);
    const auto& grid_view = tmp_grid_view.access();
    using GridViewType = std::decay_t<decltype(grid_view)>;
    const auto adapter = std::make_shared<VisualizationAdapterFunction<GridViewType, dimRange, dimRangeCols>>(*this);
    std::unique_ptr<VTKWriter<GridViewType>> vtk_writer =
        subsampling ? Common::make_unique<SubsamplingVTKWriter<GridViewType>>(grid_view, /*subsampling_level=*/2)
                    : Common::make_unique<VTKWriter<GridViewType>>(grid_view, VTK::nonconforming);
    vtk_writer->addVertexData(adapter);
    Common::test_create_directory(directory);
    if (MPIHelper::getCollectiveCommunication().size() == 1)
      vtk_writer->write(path, vtk_output_type);
    else
      vtk_writer->pwrite(Common::filename_only(path), directory, "", vtk_output_type);
  } // ... visualize(...)

  virtual void report(std::ostream& out, const std::string prefix = "") const
  {
    out << prefix << "function '" << name() << "' (of type " << type() << ")";
  }

private:
  template <class T>
  friend std::ostream& operator<<(std::ostream& /*out*/, const ThisType& /*function*/);
}; // class LocalizableFunctionInterface


template <class E, class D, size_t d, class R, size_t r, size_t rC>
std::ostream& operator<<(std::ostream& out, const LocalizableFunctionInterface<E, D, d, R, r, rC>& function)
{
  function.report(out);
  return out;
} // ... operator<<(...)


} // namespace Functions
} // namespace XT
} // namespace Dune

#include "../combined.hh"
#include "../default.hh"
#include "../derived.hh"

#endif // DUNE_XT_FUNCTIONS_INTERFACES_LOCALIZABLE_FUNCTION_HH
