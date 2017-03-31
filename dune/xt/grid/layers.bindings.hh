#include "layers.hh"

namespace Dune {
namespace XT {
namespace Grid {
namespace bindings {


template <Layers layer>
struct layer_name
{
  static_assert(AlwaysFalse<typename XT::Grid::internal::layer_dependent_typename<layer>::type>::value,
                "Please add a specialization for this layer!");

  static std::string value()
  {
    return "";
  }
};

template <>
struct layer_name<Layers::adaptive_leaf>
{
  static std::string value()
  {
    return "adaptive_leaf";
  }
};

template <>
struct layer_name<Layers::dd_subdomain>
{
  static std::string value()
  {
    return "dd_subdomain";
  }
};

template <>
struct layer_name<Layers::leaf>
{
  static std::string value()
  {
    return "leaf";
  }
};

template <>
struct layer_name<Layers::level>
{
  static std::string value()
  {
    return "level";
  }
};


} // namespace bindings
} // namespace Grid
} // namespace XT
} // namespace Dune
