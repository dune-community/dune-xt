// This file is part of the dune-xt-common project:
//   https://github.com/dune-community/dune-xt-common
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014, 2016)
//   Rene Milk       (2013, 2015)

#include "config.h"

#include "signals.hh"

#include <dune/xt/common/logging.hh>
#include <dune/xt/common/string.hh>

namespace Dune {
namespace XT {
namespace Common {

//! reset given signal to default handler
void reset_signal(int signal)
{
  struct sigaction new_action;

  new_action.sa_handler = SIG_DFL;
  sigemptyset(&new_action.sa_mask);
  new_action.sa_flags = 0;
  sigaction(signal, &new_action, NULL);
} // reset_signal

//! example signal handler
void handle_interrupt(int signal)
{
  DXTC_LOG_INFO << "forcefully terminated at " << stringFromTime() << std::endl;
  // reset signal handler and commit suicide
  reset_signal(signal);
  kill(getpid(), signal);
} // handle_interrupt

//! type of handler functions
typedef void handler_type(int);

//! calling this from your main() will install handler as callback when signal is received
void install_signal_handler(int signal, handler_type handler)
{
  struct sigaction new_action;

  /* Set up the structure to specify the new action. */
  new_action.sa_handler = handler;
  sigemptyset(&new_action.sa_mask);
  new_action.sa_flags = 0;

  sigaction(signal, &new_action, NULL);
} // install_signal_handler

} // namepsace Common
} // namespace XT
} // namepsace Dune
