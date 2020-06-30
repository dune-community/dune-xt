# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017, 2019)
#   René Fritze     (2017 - 2019)
#   Tobias Leibner  (2018 - 2020)
# ~~~

import os
import pprint


def parse_cache(filepath):
    import pyparsing as p
    EQ = p.Literal('=').suppress()
    DP = p.Literal(':').suppress()
    token = p.Word(p.alphanums + '/_- .')
    line = p.Group(token('key') + DP + token('type') + EQ + p.restOfLine('value'))
    line.ignore(p.pythonStyleComment)
    line.ignore(p.dblSlashComment)
    grammar = p.OneOrMore(line)

    kv = {}
    types = {}
    for key, type, value in grammar.parseFile(filepath, parseAll=True):
        kv[key] = value.strip()
        types[key] = type
        if key.endswith('_DIR'):
            kv[key[:-4]] = os.path.isdir(value)
    return kv, types


if __name__ == '__main__':
    kv, types = parse_cache('/tmp/CMakeCache.txt')
    pprint.pprint(kv)
