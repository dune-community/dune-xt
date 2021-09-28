// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012, 2014, 2016 - 2017)
//   René Fritze     (2010, 2012 - 2013, 2015 - 2016, 2018 - 2020)
//   Tobias Leibner  (2020)

#ifndef DUNE_XT_COMMON_SIGNALS
#define DUNE_XT_COMMON_SIGNALS

#include <signal.h>

namespace Dune::XT::Common {

//! POSIX signal management utils

//! reset given signal to default handler
void reset_signal(int signal);

//! example signal handler
void handle_interrupt(int signal);

//! type of handler functions
typedef void handler_type(int);

//! calling this from your main() will install handler as callback when signal is received
void install_signal_handler(int signal = SIGINT, handler_type handler = handle_interrupt);

} // namespace Dune::XT::Common

#endif // DUNE_XT_COMMON_SIGNALS
