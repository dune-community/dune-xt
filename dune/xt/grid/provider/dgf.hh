// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Andreas Buhr    (2014)
//   Felix Schindler (2012 - 2015)
//   Kirsten Weber   (2012)
//   Rene Milk       (2012 - 2013, 2015)
//   Tobias Leibner  (2014)

#ifndef DUNE_XT_GRID_PROVIDER_DGF_HH
#define DUNE_XT_GRID_PROVIDER_DGF_HH

#include <memory>

#if HAVE_DUNE_GRID
#include <dune/grid/io/file/dgfparser/gridptr.hh> // How convenient that GridPtr requires DGFGridFactory but
#include <dune/grid/io/file/dgfparser/dgfgridfactory.hh> // does not include it!
#include <dune/grid/io/file/dgfparser/dgfoned.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#if HAVE_ALUGRID
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#endif
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/dgf.hh>
#endif
#endif // HAVE_DUNE_GRID

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/grid/provider/interface.hh>

namespace Dune {
namespace Stuff {
namespace Grid {
namespace Providers {

#if HAVE_DUNE_GRID

template <class GridImp>
class DGF : public Grid::ProviderInterface<GridImp>
{
  typedef Grid::ProviderInterface<GridImp> BaseType;
  typedef DGF<GridImp> ThisType;

public:
  using typename BaseType::GridType;
  using BaseType::dimDomain;

  static const std::string static_id()
  {
    return BaseType::static_id() + ".dgf";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  { // size_t(...) required, else linker error with clang
    Common::Configuration config("filename", "dgf_" + Common::toString(size_t(dimDomain)) + "d_interval.dgf");
    if (sub_name.empty())
      return config;
    else {
      Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr<ThisType> create(const Common::Configuration config = default_config(),
                                          const std::string sub_name = static_id())
  {
    const Common::Configuration cfg         = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration default_cfg = default_config();
    return Common::make_unique<ThisType>(cfg.get("filename", default_cfg.get<std::string>("filename")));
  }

  DGF(const std::string filename)
  {
    grid_ = std::shared_ptr<GridType>(GridPtr<GridType>(filename).release());
  }

  DGF(ThisType&& source) = default;
  DGF(const ThisType& other) = default;

  virtual ~DGF() = default;

  ThisType& operator=(const ThisType& other) = default;
  ThisType& operator=(ThisType&& source) = default;

  virtual const GridType& grid() const override final
  {
    return *grid_;
  }

  virtual GridType& grid() override final
  {
    return *grid_;
  }

  const std::shared_ptr<const GridType> grid_ptr() const
  {
    return grid_;
  }

  std::shared_ptr<GridType> grid_ptr()
  {
    return grid_;
  }

private:
  std::shared_ptr<GridType> grid_;
}; // class DGF

#else // HAVE_DUNE_GRID

template <class GridImp>
class DGF
{
  static_assert(AlwaysFalse<GridImp>::value, "You are missing dune-grid!");
};

#endif // HAVE_DUNE_GRID

} // namespace Providers
} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // DUNE_XT_GRID_PROVIDER_DGF_HH
