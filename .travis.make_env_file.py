#!/usr/bin/env python

import os
from os.path import expanduser
home = expanduser("~")

prefixes = os.environ.get('ENV_PREFIXES', 'TRAVIS CI encrypt TOKEN TESTS').split(' ')
env_file = os.environ.get('ENV_FILE', os.path.join(home, 'env'))
with open(env_file, 'wt') as env:
    for k,v in os.environ.items():
        for pref in prefixes:
            if k.startswith(pref):
                env.write('{}="{}"\n'.format(k,v))