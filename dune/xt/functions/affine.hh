// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2013 - 2015)
//   Kirsten Weber   (2013)
//   Rene Milk       (2013 - 2015)
//   Tobias Leibner  (2014 - 2015)

#ifndef DUNE_STUFF_FUNCTIONS_AFFINE_HH
#define DUNE_STUFF_FUNCTIONS_AFFINE_HH

#include <memory>

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/functions/constant.hh>

#include "interfaces.hh"

namespace Dune {
namespace Stuff {
namespace Functions {

/**
 * \brief Simple affine function of the form f(x) = A*x + b
 */
template <class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim,
          size_t rangeDimCols = 1>
class Affine
{
  Affine()
  {
    static_assert(AlwaysFalse<EntityImp>::value, "Not available for rangeDimCols > 1!");
  }
};

template <class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim>
class Affine<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
    : public GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
{
  typedef GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> BaseType;
  typedef Affine<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> ThisType;

public:
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;
  typedef typename BaseType::JacobianRangeType JacobianRangeType;
  typedef typename Dune::FieldMatrix<RangeFieldImp, rangeDim, domainDim> MatrixType;

  using typename BaseType::LocalfunctionType;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::static_id() + ".affine";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["A"]    = internal::Get<RangeFieldImp, rangeDim, domainDim>::value_str();
    config["b"]    = internal::Get<RangeFieldImp, rangeDim, 1>::value_str();
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
    return Common::make_unique<ThisType>(cfg.get("A", default_cfg.get<MatrixType>("A")),
                                         cfg.get("b", default_cfg.get<RangeType>("b")),
                                         cfg.get("name", default_cfg.get<std::string>("name")));
  } // ... create(...)

  explicit Affine(const MatrixType& matrix, const RangeType& vector = RangeType(0),
                  const std::string name_in = static_id())
    : A_(matrix)
    , b_(vector)
    , name_(name_in)
  {
  }

  Affine(const ThisType& other) = default;

  virtual std::string type() const override final
  {
    return BaseType::static_id() + ".affine";
  }

  virtual size_t order() const override final
  {
    return 1;
  }

  virtual void evaluate(const DomainType& x, RangeType& ret) const override final
  {
    A_.mv(x, ret);
    ret += b_;
  }

  virtual void jacobian(const DomainType& /*x*/, JacobianRangeType& ret) const override final
  {
    ret = A_;
  }

  virtual std::string name() const override final
  {
    return name_;
  }

private:
  const MatrixType A_;
  const RangeType b_;
  const std::string name_;
};

} // namespace Functions
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTIONS_AFFINE_HH
