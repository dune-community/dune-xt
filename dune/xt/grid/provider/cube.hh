// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Andreas Buhr     (2014)
//   Barbara Verfürth (2015)
//   Felix Schindler  (2012 - 2015)
//   Kirsten Weber    (2012)
//   Rene Milk        (2012 - 2015)
//   Tobias Leibner   (2014)

#ifndef DUNE_XT_GRID_PROVIDER_CUBE_HH
#define DUNE_XT_GRID_PROVIDER_CUBE_HH

#include <array>
#include <memory>
#include <sstream>
#include <type_traits>
#include <vector>

#include <boost/numeric/conversion/cast.hpp>

#if HAVE_DUNE_GRID // clang-format off
# define DUNE_AVOID_SGRID_DEPRE_WARNING_BECAUSE_I_KNOW_WHAT_IM_DOING 1
# include <dune/grid/sgrid.hh>
# undef DUNE_AVOID_SGRID_DEPRE_WARNING_BECAUSE_I_KNOW_WHAT_IM_DOING
# include <dune/grid/yaspgrid.hh>
# if HAVE_ALBERTA
#   include <dune/xt/common/disable_warnings.hh>
#     include <dune/grid/albertagrid.hh>
#   include <dune/xt/common/reenable_warnings.hh>
# endif
# if HAVE_DUNE_ALUGRID
#   include <dune/alugrid/grid.hh>
# endif
# if HAVE_ALUGRID
#   include <dune/grid/alugrid.hh>
# endif
# if HAVE_DUNE_SPGRID
#   include <dune/grid/spgrid.hh>
# endif
# include <dune/xt/grid/structuredgridfactory.hh>
#endif // HAVE_DUNE_GRID // clang-format on

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/memory.hh>

#include "default.hh"

namespace Dune {
namespace XT {
namespace Grid {
namespace Providers {

#if HAVE_DUNE_GRID

namespace Configs {

static inline Common::Configuration Cube_default(const std::string sub_name = "")
{
  Common::Configuration config;
  config["lower_left"]      = "[0 0 0 0]";
  config["upper_right"]     = "[1 1 1 1]";
  config["num_elements"]    = "[8 8 8 8]";
  config["num_refinements"] = "0";
  config["overlap"] = "1";
  if (sub_name.empty())
    return config;
  else {
    Common::Configuration tmp;
    tmp.add(config, sub_name);
    return tmp;
  }
} // ... Cube_default(...)

} // namespace Configs
namespace internal {

template <typename GridType>
struct ElementVariant;

template <typename GridType>
struct ElementVariant
{
  static const int id = 2;
};

template <int dd, int dw>
struct ElementVariant<Dune::SGrid<dd, dw>>
{
  static const int id = 1;
};

template <int dim, class Coords>
struct ElementVariant<Dune::YaspGrid<dim, Coords>>
{
  static const int id = 1;
};

#if HAVE_DUNE_SPGRID
template <class ct, int dim, template <int> class Refinement, class Comm>
struct ElementVariant<Dune::SPGrid<ct, dim, Refinement, Comm>>
{
  static const int id = 1;
};
#endif

#if HAVE_ALUGRID || HAVE_DUNE_ALUGRID

template <int dimGrid, int dimWorld, class MpiCommImp>
struct ElementVariant<Dune::ALUGrid<dimGrid, dimWorld, Dune::cube, Dune::conforming, MpiCommImp>>
{
  static const int id = 1;
};

template <int dimGrid, int dimWorld, class MpiCommImp>
struct ElementVariant<Dune::ALUGrid<dimGrid, dimWorld, Dune::cube, Dune::nonconforming, MpiCommImp>>
{
  static const int id = 1;
};

#endif // HAVE_ALUGRID || HAVE_DUNE_ALUGRID

} // namespace internal

/**
 *  \brief  Creates a grid of a cube in various dimensions.
 *
 *          Default implementation using the Dune::StructuredGridFactory to create a grid of a cube in 1, 2 or 3
 *          dimensions. Tested with
 *          <ul><li> \c YASPGRID, \c variant 1, dim = 1, 2, 3,
 *          <li> \c SGRID, \c variant 1, dim = 1, 2, 3,
 *          <li> \c ALUGRID_SIMPLEX, \c variant 2, dim = 2, 3,
 *          <li> \c ALUGRID_CONFORM, \c variant 2, dim = 2, 2 and
 *          <li> \c ALUGRID_CUBE, \c variant 1, dim = 2, 3.</ul>
 *  \tparam GridImp
 *          Type of the underlying grid.
 *  \tparam variant
 *          Type of the codim 0 elements:
 *          <ul><li>\c 1: cubes
 *          <li>2: simplices</ul>
 **/
template <typename GridImp, int variant = internal::ElementVariant<GridImp>::id>
class Cube : public ProviderInterface<GridImp>
{
  typedef ProviderInterface<GridImp> BaseType;
  typedef Cube<GridImp, variant> ThisType;

public:
  using typename BaseType::GridType;
  static const size_t dimDomain = BaseType::dimDomain;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;

  static const std::string static_id() { return BaseType::static_id() + ".cube"; }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    return Configs::Cube_default(sub_name);
  }

  static std::unique_ptr<ThisType> create(const Common::Configuration config = default_config(),
                                          const std::string sub_name = static_id())
  {
    // get correct config
    const Common::Configuration cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const auto overlap              = cfg.get("overlap", default_config().get<unsigned int>("overlap"));
    const auto overlap_array = Common::make_array<unsigned int, dimDomain>(overlap);
    return Common::make_unique<ThisType>(
        cfg.get("lower_left", default_config().get<DomainType>("lower_left")),
        cfg.get("upper_right", default_config().get<DomainType>("upper_right")),
        cfg.get("num_elements", default_config().get<std::vector<unsigned int>>("num_elements"), dimDomain),
        cfg.get("num_refinements", default_config().get<size_t>("num_refinements")),
        overlap_array);
  } // ... create(...)

  /**
   *  \brief      Creates a cube.
   *  \param[in]  lower_left
   *              Used as a lower left corner (in each dimension, if scalar).
   *  \param[in]  upper_right
   *              Used as an upper right corner (in each dimension, if scalar).
   *  \param[in]  num_elements (optional)
   *              Number of elements.
   **/
  explicit Cube(const DomainFieldType lower_left = default_config().get<DomainType>("lower_left")[0],
                const DomainFieldType upper_right = default_config().get<DomainType>("upper_right")[0],
                const unsigned int num_elements   = default_config().get<std::vector<unsigned int>>("num_elements")[0],
                const size_t num_refinements = default_config().get<size_t>("num_refinements"),
                const std::array<unsigned int, dimDomain> overlap =
                    Common::make_array<unsigned int, dimDomain>(default_config().get<unsigned int>("overlap")))
    : grid_ptr_(create_grid(DomainType(lower_left), DomainType(upper_right), parse_array(num_elements), num_refinements,
                            overlap))
  {
  }

  Cube(const Common::FieldVector<DomainFieldType, dimDomain>& lower_left,
       const Common::FieldVector<DomainFieldType, dimDomain>& upper_right,
       const unsigned int num_elements = default_config().get<std::vector<unsigned int>>("num_elements")[0],
       const size_t num_refinements = default_config().get<size_t>("num_refinements"),
       const std::array<unsigned int, dimDomain> overlap =
           Common::make_array<unsigned int, dimDomain>(default_config().get<unsigned int>("overlap")))
    : grid_ptr_(create_grid(lower_left, upper_right, parse_array(num_elements), num_refinements, overlap))
  {
  }

  Cube(const Common::FieldVector<DomainFieldType, dimDomain>& lower_left,
       const Common::FieldVector<DomainFieldType, dimDomain>& upper_right, const std::vector<unsigned int> num_elements,
       const size_t num_refinements = default_config().get<size_t>("num_refinements"),
       const std::array<unsigned int, dimDomain> overlap =
           Common::make_array<unsigned int, dimDomain>(default_config().get<unsigned int>("overlap")))
    : grid_ptr_(create_grid(lower_left, upper_right, parse_array(num_elements), num_refinements, overlap))
  {
  }

  virtual GridType& grid() override { return *grid_ptr_; }

  virtual const GridType& grid() const override { return *grid_ptr_; }

  std::shared_ptr<GridType> grid_ptr() { return grid_ptr_; }

  const std::shared_ptr<const GridType>& grid_ptr() const { return grid_ptr_; }

private:
  static std::array<unsigned int, dimDomain> parse_array(const unsigned int in)
  {
    std::array<unsigned int, dimDomain> ret;
    std::fill(ret.begin(), ret.end(), in);
    return ret;
  } // ... parse_array(...)

  static std::array<unsigned int, dimDomain> parse_array(const std::vector<unsigned int>& in)
  {
    if (in.size() < dimDomain)
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "Given vector is too short: should be " << dimDomain << ", is " << in.size() << "!");
    std::array<unsigned int, dimDomain> ret;
    for (size_t ii = 0; ii < dimDomain; ++ii)
      ret[ii] = in[ii];
    return ret;
  } // ... parse_array(...)

  /// TODO simplex grid overlap
  static std::shared_ptr<GridType>
      create_grid(DomainType lower_left, DomainType upper_right,
                  const std::array<unsigned int, dimDomain>& num_elements, const size_t num_refinements,
                  const std::array<unsigned int, dimDomain> overlap =
                      Common::make_array<unsigned int, dimDomain>(default_config().get<size_t>("overlap")))
  {
    static_assert(variant == 1 || variant == 2, "variant has to be 1 or 2!");
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      if (!(lower_left[dd] < upper_right[dd]))
        DUNE_THROW(Common::Exceptions::wrong_input_given,
                   "lower_left has to be elementwise smaller than upper_right!\n\n" << lower_left[dd] << " vs. "
                                                                                    << upper_right[dd]);
    }
    std::shared_ptr<GridType> grd_ptr(nullptr);
    switch (variant) {
      case 1:
        grd_ptr = StructuredGridFactory<GridType>::createCubeGrid(lower_left, upper_right, num_elements, overlap);
        break;
      case 2:
      default:
        grd_ptr = StructuredGridFactory<GridType>::createSimplexGrid(lower_left, upper_right, num_elements);
        break;
    }
    grd_ptr->loadBalance();
#if HAVE_ALBERTA
    if (!std::is_same<GridType, AlbertaGrid<dimDomain, dimDomain>>::value)
#endif
      grd_ptr->preAdapt();
    grd_ptr->globalRefine(boost::numeric_cast<int>(num_refinements));
#if HAVE_ALBERTA
    if (!std::is_same<GridType, AlbertaGrid<dimDomain, dimDomain>>::value)
#endif
      grd_ptr->postAdapt();
    grd_ptr->loadBalance();
    return grd_ptr;
  } // ... create_grid(...)

  std::shared_ptr<GridType> grid_ptr_;
}; // class Cube

#else // HAVE_DUNE_GRID

template <typename GridImp, int variant = 1>
class Cube
{
  static_assert(AlwaysFalse<GridImp>::value, "You are missing dune-grid!");
};

#endif // HAVE_DUNE_GRID

} // namespace Providers
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_PROVIDER_CUBE_HH
