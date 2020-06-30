#!/usr/bin/env python3
#
# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017)
#   René Fritze     (2017 - 2019)
#   Tobias Leibner  (2019 - 2020)
# ~~~
"""Run a speedup series benchmark of a given executable on a node with libnuma

Usage: numa_speedup.py MIN_CORE MAX_CORE ARG_FMT...

Arguments:
    ARG_FMT         A python format string to allow core count injection in args.
                    Example: 'mpirun -np {CORECOUNT} my_executable -global.datadir data_{CORECOUNT}'
    MIN_CORE        Start run at this many cores
    MAX_CORE        Stop run at this many cores
"""

import subprocess
import docopt
import itertools
from collections import deque
from copy import deepcopy


def coremap():
    try:
        import numa
    except ImportError:
        print('This script requires the libnuma python bindings')
        raise RuntimeError("Numa not available")
    if not numa.available():
        raise RuntimeError("Numa not available")
    node_to_core = {int(i): deque([int(k) for k in numa.node_to_cpus(i)]) for i in range(numa.get_max_node() + 1)}
    total_core = max(itertools.chain(*node_to_core.values())) + 1
    return node_to_core, total_core


def _single_run(cores, arg_fmt):
    args = [a.format(CORES=cores, CORECOUNT=len(cores)) for a in arg_fmt]
    corestring = ','.join([str(i) for i in cores])
    subprocess.call(['numactl', '-C', corestring] + args)


def _assign_core(node_to_core, core_count):
    cores = []
    node_to_core = deepcopy(node_to_core)
    nodes = deque([k for k in node_to_core.keys()])
    left = True
    while core_count > 0:
        if left:
            cores.append(node_to_core[nodes[0]].popleft())
        else:
            cores.append(node_to_core[nodes[0]].pop())
        left = not left
        nodes.rotate(1)
        core_count -= 1
    return cores


def run(min_core, max_core, arg_fmt):
    node_to_core, total_cores = coremap()
    core_count = min_core
    max_core = min(max_core, total_cores)
    while core_count <= max_core:
        cores = _assign_core(node_to_core, core_count)
        _single_run(cores, arg_fmt)
        core_count *= 2


if __name__ == '__main__':
    arguments = docopt.docopt(__doc__)
    min_core = int(arguments['MIN_CORE'])
    max_core = int(arguments['MAX_CORE'])
    arg_fmt = arguments['ARG_FMT']
    run(min_core, max_core, arg_fmt)
