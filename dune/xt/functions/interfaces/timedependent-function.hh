#ifndef DUNE_XT_FUNCTIONS_INTERFACES_TIMEDEPENDENT_FUNCTION_HH
#define DUNE_XT_FUNCTIONS_INTERFACES_TIMEDEPENDENT_FUNCTION_HH

#include <memory>
#include <string>

namespace Dune {
namespace XT {
namespace Functions {


template <class TimeIndependentFunctionImp, class TimeFieldImp = double>
class TimeDependentFunctionInterface
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
