// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_FUNCTIONS_VISUALIZATION_HH
#define DUNE_STUFF_FUNCTIONS_VISUALIZATION_HH

#if HAVE_DUNE_GRID
#include <dune/grid/io/file/vtk/function.hh>
#endif

#include <dune/stuff/common/float_cmp.hh>
//#include <dune/stuff/grid/entity.hh>

#include "interfaces.hh"

namespace Dune {
namespace Stuff {
namespace Function {


#if HAVE_DUNE_GRID
template <class GridViewType, int dimRange>
class VisualizationAdapter : public VTKFunction<GridViewType>
{
public:
  typedef typename GridViewType::template Codim<0>::Entity EntityType;

  typedef typename GridViewType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridViewType::dimension;
  typedef FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, double, dimRange> FunctionType;

  VisualizationAdapter(const FunctionType& function, const std::string nm = "")
    : function_(function)
    , tmp_value_(0)
    , name_(nm)
  {
  }

  virtual int ncomps() const DS_OVERRIDE
  {
    return dimRange;
  }

  virtual std::string name() const DS_OVERRIDE
  {
    if (name_.empty())
      return function_.name();
    else
      return name_;
  }

  virtual double evaluate(int comp, const EntityType& en, const DomainType& xx) const DS_OVERRIDE
  {
    assert(comp >= 0);
    assert(comp < dimRange);
    const auto local_func = function_.local_function(en);
    local_func->evaluate(xx, tmp_value_);
    return tmp_value_[comp];
  }

private:
  const FunctionType& function_;
  mutable FieldVector<double, dimRange> tmp_value_;
  const std::string name_;
}; // class VisualizationAdapter
#endif // HAVE_DUNE_GRID


// template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int
// rangeDimCols = 1 >
// class LocalDifferentiableFunctionDefault;


// template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp >
// class LocalDifferentiableFunctionDefault< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1 >
//  : public LocalfunctionInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1 >
//{
//  typedef LocalfunctionInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1 > BaseType;
// public:
//  typedef typename BaseType::EntityType EntityType;

//  typedef typename BaseType::DomainFieldType  DomainFieldType;
//  static const unsigned int                   dimDomain = BaseType::dimDomain;
//  typedef typename BaseType::DomainType       DomainType;

//  typedef typename BaseType::RangeFieldType RangeFieldType;
//  static const unsigned int                 dimRange = BaseType::dimRange;
//  static const unsigned int                 dimRangeRows = BaseType::dimRangeRows;
//  static const unsigned int                 dimRangeCols = BaseType::dimRangeCols;
//  typedef typename BaseType::RangeType      RangeType;

//  typedef typename BaseType::JacobianRangeType JacobianRangeType;

//  LocalDifferentiableFunctionDefault(const EntityType& ent)
//    : entity_(ent)
////    , max_variation_(DomainFieldType(0))
//    , disturbed_point_(DomainFieldType(0))
//    , point_value_(RangeFieldType(0))
//    , disturbed_value_(RangeFieldType(0))
//    , hh_(std::max(100.0 * std::numeric_limits< DomainFieldType >::min(), 1e-2))
//  {
////    // get corners
////    const size_t num_corners = entity_.geometry().corners();
////    std::vector< DomainType > corners(num_corners, DomainType(DomainFieldType(0)));
////    for (size_t cc = 0; cc < num_corners; ++cc)
////      corners[cc] = entity_.geometry().corner(cc);
////    // compute max distance per dimension
////    DomainType difference(DomainFieldType(0));
////    for (size_t cc = 0; cc < num_corners; ++cc) {
////      const auto& corner = corners[cc];
////      for (size_t oo = cc + 1; oo < num_corners; ++oo) {
////        const auto& other_corner = corners[oo];
////        difference = corner;
////        difference -= other_corner;
////        max_variation_[cc] = std::max(max_variation_[cc], std::abs(difference[cc]));
////      }
////    }
//    assert(Common::FloatCmp::ne(hh_, DomainFieldType(0)));
//  } // LocalDifferentiableFunctionDefault(...)

//  virtual ~LocalDifferentiableFunctionDefault() {}

//  virtual const EntityType& entity() const DS_OVERRIDE
//  {
//    return entity_;
//  }

//  virtual void evaluate(const DomainType& /*xx*/, RangeType& /*ret*/) const = 0;

//  virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const DS_OVERRIDE
//  {
//    assert(this->is_a_valid_point(xx));
//    // clear
//    ret *= RangeFieldType(0);
//    // evaluate
//    evaluate(xx, point_value_);
//    // loop over all dimensions
//    for (size_t dd = 0; dd < dimDomain; ++dd) {
//      // disturbe the point
//      disturbed_point_ = xx;
//      disturbed_point_[dd] += hh_;
//      assert(this->is_a_valid_point(disturbed_point_));
////      // find a closer point if that one is not contained in the entity
////      size_t num_tries = 0;
////      while(!this->is_a_valid_point(disturbed_point_) && num_tries < max_tries_) {
////        disturbed_point_[dd] = 0.5*(xx[dd] + disturbed_point_[dd]);
////        ++num_tries;
////      }
////      if (num_tries == max_tries_)
////        DUNE_THROW(InvalidStateException, "Could not find a disturbed point inside the entity!");
//      // compute gradient
//      evaluate(disturbed_point_, disturbed_value_);
////      const DomainFieldType hh = std::abs(disturbed_point_[dd] - xx[dd]);
//      ret[0][dd] = (point_value_ - disturbed_value_) / hh_;
//    } // loop over all dimensions
//  } // ... jacobian(...)

// private:
//  const EntityType& entity_;
////  const size_t max_tries_;
//  DomainFieldType hh_;
////  DomainType max_variation_;
//  mutable DomainType disturbed_point_;
//  mutable RangeType point_value_;
//  mutable RangeType disturbed_value_;
//}; // class LocalDifferentiableFunctionDefault< ..., 1 >


} // namespace Function
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTIONS_VISUALIZATION_HH
