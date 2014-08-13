// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Kirsten Weber

#ifndef DUNE_STUFF_FUNCTIONS_CONSTANT_HH
#define DUNE_STUFF_FUNCTIONS_CONSTANT_HH

#include <memory>

#include <dune/stuff/common/configuration.hh>

#include "interfaces.hh"

namespace Dune {
namespace Stuff {
namespace Functions {


template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class Constant
    : public GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
{
  typedef GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> BaseType;
  typedef Constant<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> ThisType;

  template <class R, int r, int rC>
  struct Get
  {
    static std::string value_str()
    {
      std::string str = "[";
      for (size_t rr = 0; rr < r; ++rr) {
        if (rr > 0)
          str += "; ";
        for (size_t cc = 0; cc < rC; ++cc) {
          if (cc > 0)
            str += " ";
          if (cc == rr)
            str += "1";
          else
            str += "0";
        }
      }
      str += "]";
      return str;
    }
  };

  template <class R, int rC>
  struct Get<R, 1, rC>
  {
    static std::string value_str()
    {
      std::string str = "[";
      for (size_t cc = 0; cc < rC; ++cc) {
        if (cc > 0)
          str += " ";
        str += "1";
      }
      str += "]";
      return str;
    }
  };

  template <class R, int r>
  struct Get<R, r, 1>
  {
    static std::string value_str()
    {
      return Get<R, 1, r>::value_str();
    }
  };

  template <class R>
  struct Get<R, 1, 1>
  {
    static std::string value_str()
    {
      return "1";
    }
  };

public:
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeFieldType RangeFieldType;
  typedef typename BaseType::RangeType RangeType;
  typedef typename BaseType::JacobianRangeType JacobianRangeType;

  using typename BaseType::LocalfunctionType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".constant";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["value"] = Get<RangeFieldImp, rangeDim, rangeDimCols>::value_str();
    config["name"] = static_id();
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
    // get correct config
    const Common::Configuration cfg         = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration default_cfg = default_config();
    return Common::make_unique<ThisType>(cfg.get("value", default_cfg.get<RangeType>("value")),
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

  virtual std::string type() const DS_OVERRIDE DS_FINAL
  {
    return BaseType::static_id() + ".constant";
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

#if 0 // def DUNE_STUFF_FUNCTIONS_TO_LIB
the whole  DUNE_STUFF_FUNCTIONS_TO_LIB stuff threw loads of errors, so I disabled it everywhere (Rene)

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
