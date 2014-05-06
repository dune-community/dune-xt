// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Kirsten Weber

#ifndef DUNE_STUFF_FUNCTIONS_CONSTANT_HH
#define DUNE_STUFF_FUNCTIONS_CONSTANT_HH

#include <memory>

#include <dune/stuff/common/configtree.hh>

#include "interfaces.hh"

namespace Dune {
namespace Stuff {
namespace Functions {


template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class Constant
    : public GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
{

public:
  typedef GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> BaseType;
  typedef Constant<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> ThisType;
  typedef typename BaseType::RangeType RangeType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::JacobianRangeType JacobianRangeType;

  using typename BaseType::LocalfunctionType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".constant";
  }

  static Common::ConfigTree default_config(const std::string sub_name = "")
  {
    Common::ConfigTree config;
    config["value"] = "1.0";
    config["name"] = static_id();
    if (sub_name.empty())
      return config;
    else {
      Common::ConfigTree tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr<ThisType> create(const Common::ConfigTree config = default_config(),
                                          const std::string sub_name = static_id())
  {
    // get correct config
    const Common::ConfigTree cfg         = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::ConfigTree default_cfg = default_config();
    // create
    return Common::make_unique<ThisType>(cfg.get<RangeType>("value", default_cfg.get<RangeType>("value")),
                                         cfg.get("name", default_cfg.get<std::string>("name")));
  } // ... create(...)

  explicit Constant(const RangeType& constant, const std::string name = static_id())
    : constant_(constant)
    , name_(name)
  {
  }

  explicit Constant(const RangeFieldImp& constant, const std::string name = static_id())
    : constant_(constant)
    , name_(name)
  {
  }

  Constant(const ThisType& other)
    : constant_(other.constant_)
    , name_(other.name_)
  {
  }

  virtual size_t order() const DS_OVERRIDE DS_FINAL
  {
    return 0;
  }

  virtual void evaluate(const DomainType& /*x*/, RangeType& ret) const DS_OVERRIDE DS_FINAL
  {
    ret = constant_;
  }

  virtual void jacobian(const DomainType& /*x*/, JacobianRangeType& ret) const DS_OVERRIDE DS_FINAL
  {
    ret *= 0.0;
  }

  virtual std::string name() const DS_OVERRIDE DS_FINAL
  {
    return name_;
  }

private:
  const RangeType constant_;
  const std::string name_;
};


} // namespace Functions
} // namespace Stuff
} // namespace Dune

#ifdef DUNE_STUFF_FUNCTIONS_TO_LIB
#define DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGE(etype, ddim)                                                       \
  DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGECOLS(etype, ddim, 1)                                                      \
  DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGECOLS(etype, ddim, 2)                                                      \
  DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGECOLS(etype, ddim, 3)

#define DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGECOLS(etype, ddim, rdim)                                             \
  DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DOMAINFIELDTYPES(etype, ddim, rdim, 1)                                            \
  DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DOMAINFIELDTYPES(etype, ddim, rdim, 2)                                            \
  DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DOMAINFIELDTYPES(etype, ddim, rdim, 3)

#define DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DOMAINFIELDTYPES(etype, ddim, rdim, rcdim)                                  \
  DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_RANGEFIELDTYPES(etype, double, ddim, rdim, rcdim)

#define DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_RANGEFIELDTYPES(etype, dftype, ddim, rdim, rcdim)                           \
  DUNE_STUFF_FUNCTIONS_CONSTANT_LAST_EXPANSION(etype, dftype, ddim, double, rdim, rcdim)                               \
  DUNE_STUFF_FUNCTIONS_CONSTANT_LAST_EXPANSION(etype, dftype, ddim, long double, rdim, rcdim)

#define DUNE_STUFF_FUNCTIONS_CONSTANT_LAST_EXPANSION(etype, dftype, ddim, rftype, rdim, rcdim)                         \
  extern template class Dune::Stuff::Functions::Constant<etype, dftype, ddim, rftype, rdim, rcdim>;

DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGE(DuneStuffFunctionsInterfacesFake1dEntityType, 1)
DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGE(DuneStuffFunctionsInterfacesFake2dEntityType, 2)
DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGE(DuneStuffFunctionsInterfacesFake3dEntityType, 3)

#if HAVE_DUNE_GRID

DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGE(DuneStuffFunctionsInterfacesSGrid1dEntityType, 1)
DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGE(DuneStuffFunctionsInterfacesSGrid2dEntityType, 2)
DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGE(DuneStuffFunctionsInterfacesSGrid3dEntityType, 3)

DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGE(DuneStuffFunctionsInterfacesYaspGrid1dEntityType, 1)
DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGE(DuneStuffFunctionsInterfacesYaspGrid2dEntityType, 2)
DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGE(DuneStuffFunctionsInterfacesYaspGrid3dEntityType, 3)

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H

DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGE(DuneStuffFunctionsInterfacesAluSimplexGrid2dEntityType, 2)
DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGE(DuneStuffFunctionsInterfacesAluSimplexGrid3dEntityType, 3)
DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGE(DuneStuffFunctionsInterfacesAluCubeGrid3dEntityType, 3)

#endif // HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
#endif // HAVE_DUNE_GRID

#undef DUNE_STUFF_FUNCTIONS_CONSTANT_LAST_EXPANSION
#undef DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_RANGEFIELDTYPES
#undef DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DOMAINFIELDTYPES
#undef DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGECOLS
#undef DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGE
#endif // DUNE_STUFF_FUNCTIONS_TO_LIB

#endif // DUNE_STUFF_FUNCTIONS_CONSTANT_HH
