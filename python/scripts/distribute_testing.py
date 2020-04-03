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
#   René Fritze     (2016 - 2019)
#   Tobias Leibner  (2019 - 2020)
# ~~~

import logging
import os
import pickle
import json
import sys
from pprint import pprint, pformat
import subprocess
import time
from contextlib import contextmanager
import binpacking
from multiprocessing import Pool, cpu_count
from collections import OrderedDict
from statistics import mean, stdev

MAXTIME = 23 * 60
pickle_file = 'totals.pickle'
PROTOCOL = 0


@contextmanager
def elapsed_timer():
    clock = time.time
    start = clock()
    elapser = lambda: clock() - start
    yield lambda: elapser()
    end = clock()
    elapser = lambda: end - start


def _dump(obj, fn):
    json.dump(obj, open(fn, 'wt'), sort_keys=True)


def _load(fn):
    try:
        return json.load(open(fn, 'rt'))
    except FileNotFoundError:
        pass
    try:
        return pickle.load(open(fn, 'rb'))
    except FileNotFoundError:
        return None


def _compile(binary):
    with elapsed_timer() as timer:
        try:
            _ = subprocess.check_output(['cmake', '--build', '.', '--', '-j1', binary],
                                        universal_newlines=True,
                                        stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as cpe:
            if 'Timeout' not in cpe.output:
                raise cpe
            print('Timeout in compile {}'.format(binary))
        return timer()


def _run_tests(tpl):
    binary, teststrings = tpl
    testtimes = 0
    for test in teststrings.split(';'):
        with elapsed_timer() as timer:
            try:
                _ = subprocess.check_output(['ctest', '-j1', '-N', '-R', test],
                                            universal_newlines=True,
                                            stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as cpe:
                if 'Timeout' not in cpe.output:
                    raise cpe
                # be pessimistic and double the timeout value as time for this run
                testtimes += timer()
                print('Timeout in {} from {}'.format(test, binary))
            testtimes += timer()
    return testtimes


def _redo(processes, keys, *args):
    try:
        with Pool(processes=processes) as pool:
            result = pool.map(*args)
        return {k: v for k, v in zip(keys, result)}
    except subprocess.CalledProcessError as cpe:
        logging.error('*' * 79)
        try:
            logging.error(cpe.stdout)
            logging.error(cpe.stderr)
        except AttributeError:
            logging.error(cpe.output)
        logging.error('*' * 79)
        raise cpe


def do_timings(builddir, pickledir, binaries, testnames, processes, headerlibs):
    os.chdir(builddir)
    testlimit = -1

    binaries = binaries[:testlimit]
    headerlibs = headerlibs[:testlimit]
    targets = binaries + headerlibs
    compiles_fn = os.path.join(pickledir, 'compiles_' + pickle_file)
    try:
        compiles = _load(compiles_fn) or {}
        if set(compiles.keys()) != set(targets):
            removed = set(compiles.keys()) - set(targets)
            new = set(targets) - set(compiles.keys())
            logging.error('redoing compiles due to mismatched binaries')
            logging.error('Removed: {}'.format(pformat(removed)))
            logging.error('New: {}'.format(pformat(new)))
            compiles = {k: v for k, v in compiles.items() if k not in removed}
            compiles.update(_redo(processes, new, _compile, new))
    except FileNotFoundError:
        logging.error('redoing compiles due to missing pickle')
        compiles = _redo(processes, targets, _compile, targets)
    _dump(compiles, compiles_fn)

    testnames = testnames[:testlimit]
    testruns_fn = os.path.join(pickledir, 'testruns_' + pickle_file)
    try:
        loaded_testnames, testruns = _load(testruns_fn) or ([], dict())
        if set(compiles.keys()) != set(targets) or loaded_testnames != testnames:
            removed = [f for f in loaded_testnames if f not in set(testnames)]
            new = [n for n in testnames if n not in set(loaded_testnames)]
            logging.error('redoing tests due to mismatched binaries/testnames')
            logging.error('Removed: {}'.format(pformat(removed)))
            logging.error('New: {}'.format(pformat(new)))
            testruns = {k: v for k, v in testruns.items() if k not in removed}
            new_map = [(b, t) for b, t in zip(binaries, testnames) if t not in removed]
            testruns.update(_redo(processes, new, _run_tests, new_map))
    except FileNotFoundError:
        logging.error('redoing tests due to missing pickle')
        testruns = _redo(processes, binaries, _run_tests, zip(binaries, testnames))
    _dump((testnames, testruns), testruns_fn)

    totals = {n: compiles[n] for n in binaries}

    testname_map = {b: t.strip().split(';') for b, t in zip(binaries, testnames)}
    for binary, testnames in testname_map.items():
        totals[binary] += testruns[';'.join(testnames)]

    # add totals for headerlib compiles that do not have associated testruns
    totals.update({n: compiles[n] for n in headerlibs})
    _dump(totals, os.path.join(pickledir, pickle_file))
    # print('totals')
    # pprint.pprint(totals)
    return totals


builddir = sys.argv[1]
testdir = sys.argv[2]
cmake_outfile = os.path.join(builddir, 'builder_definitions.cmake')
binaries = open(sys.argv[3], 'rt').read().split(';')
# list comes with a leading empty entry
all_testnames = open(sys.argv[4], 'rt').read().split('/')[1:]
headerlibs = open(sys.argv[5], 'rt').read().split(';')
try:
    bincount = int(sys.argv[6])
except:
    bincount = 13
logging.basicConfig(level=logging.DEBUG)
testname_map = {b: t.strip().split(';') for b, t in zip(binaries, all_testnames)}
processes = 1     #cpu_count()

totals = do_timings(builddir, builddir, binaries, all_testnames, processes, headerlibs)

#bins = binpacking.to_constant_volume(totals, MAXTIME)
bins = binpacking.to_constant_bin_number(totals, bincount)
vols = [sum(bi.values()) for bi in bins]
norm = 100 / MAXTIME
print('Generated {} bins.\nRelative volumes:\n\t\tMin {:.2f}%\n\t\tMax {:.2f}%\n\t\tAvg {:.2f}%\n'.format(
    len(bins),
    min(vols) * norm,
    max(vols) * norm,
    mean(vols) * norm, stdev(vols)))

set_tpl = '''set_tests_properties({} PROPERTIES LABELS "builder_{}")\n'''
with open(cmake_outfile, 'wt') as out:
    out.write('set(DXT_BIN_COUNT "{}" CACHE STRING "number of bins for test targets" )\n'.format(len(bins)))
    for idx, bin in enumerate(bins):
        out.write('add_custom_target(test_binaries_builder_{} DEPENDS {})\n'.format(idx, ' '.join(sorted(bin.keys()))))
        for binary in sorted(bin.keys()):
            if binary not in headerlibs:
                out.write(set_tpl.format(' '.join(testname_map[binary]), idx))
