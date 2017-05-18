#ifndef DUNE_XT_FUNCTIONS_INTERFACES_TIMEDEPENDENT_FUNCTION_HH
#define DUNE_XT_FUNCTIONS_INTERFACES_TIMEDEPENDENT_FUNCTION_HH

#include <memory>
#include <string>

#include <dune/common/deprecated.hh>

namespace Dune {
namespace XT {
namespace Functions {


template <class TimeIndependentFunctionImp, class TimeFieldImp = double>
class DUNE_DEPRECATED_MSG("All functions are parametric by now, use '_t' as the parameter for time, see "
                          "ParametricExpressionFunction! (18.05.2017)") TimeDependentFunctionInterface
{
public:
  typedef TimeIndependentFunctionImp TimeIndependentFunctionType;
  typedef TimeFieldImp TimeFieldType;

  static const bool available = false;

  virtual ~TimeDependentFunctionInterface()
  {
  }

  static std::string static_id()
  {
    return "xt.functions.timedependentfunction";
  }

  /**
   * \name ´´These methods have to be implemented.''
   * \{
   **/
  virtual std::unique_ptr<TimeIndependentFunctionType> evaluate_at_time(const TimeFieldType t) const = 0;
  /* \} */

  /**
   * \name ´´These methods should be implemented in order to identify the function.''
   * \{
   */
  virtual std::string type() const
  {
    return "xt.functions.timedependentfunction";
  }

  virtual std::string name() const
  {
    return "xt.functions.timedependentfunction";
  }
  /* \} */
};


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_INTERFACES_TIMEDEPENDENT_FUNCTION_HH
