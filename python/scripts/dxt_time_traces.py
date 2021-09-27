#!/usr/bin/env python3
import os
import sys
import json
import pandas as pd

categories = ['Source', 'InstantiateClass', 'InstantiateFunction', 'ParseClass', 'ParseTemplate']
categories = [
    'Source',
    'InstantiateClass',
    'InstantiateFunction',
]


def read_files(filenames):
    df = None
    for fn in filenames:
        assert os.path.isfile(fn), fn

        data = pd.DataFrame(json.load(open(fn))['traceEvents'])
        df = df.append(data) if df is not None else data

    # from micro to full seconds
    df[['dur']] /= 1e6
    df['reldur'] = df[['dur']] / df[['dur']].max()
    return df


def _process_sources(df, strip_dir, detail_name):
    sources = df[df['name'] == detail_name]
    sources['detail'] = [f['detail'].replace(strip_dir, 'PROJ_ROOT/') for f in sources['args']]
    sources = sources.drop('args', axis=1)
    sources = sources.sort_values(by=['dur'], ascending=False)
    total_time = sources['dur'].sum()

    grouped = sources.groupby(['detail'])
    counts = grouped.count()['name']
    grouped = grouped[['dur']].sum().sort_values(by=['dur'], ascending=False).reset_index()
    sl = grouped.head(20)
    sl['count'] = [counts[f] for f in sl['detail']]
    sl = sl.rename(columns={'detail': detail_name, 'dur': 'cum. time (s)', 'count': 'event count'})
    return sl


def report(trace_files, out_fn_base='time_trace_results', strip_dir=''):
    df = read_files(trace_files)
    df.to_pickle(f'{out_fn_base}.dataframe.pickle.xz')
    totals = df[df['name'].str.contains('Total')]
    for cat in categories:
        with open(f'{out_fn_base}.{cat.lower()}.md', 'wt') as out:
            data = _process_sources(df, strip_dir, detail_name=cat)
            out.write(data.to_markdown(index=False))
    with open(f'{out_fn_base}.totals.md', 'wt') as out:
        out.write('Totals:\n')
        out.write(totals.groupby(['name'])[['dur']].sum().rename(columns={'dur': 'cum. time (s)'}).to_markdown())


if __name__ == '__main__':
    strip_dir = '/home/r_milk01/master_xt/'
    report(sys.argv[1:], strip_dir=strip_dir)
