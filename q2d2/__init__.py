#!/usr/bin/env python

__version__ = "0.0.0-dev"

import random
import io
from collections import defaultdict
import hashlib
import os
import shutil
import glob
import math

import marisa_trie
import numpy as np
import pandas as pd
from IPython.html import widgets
from IPython.html.widgets import interactive, fixed, IntSlider
from IPython.display import display
from scipy.optimize import minimize_scalar

import skbio
from skbio.diversity.beta import pw_distances
from skbio.stats.ordination import PCoA
from skbio.stats import subsample_counts

def exact(trie, seq):
    return [(trie[str(seq)], 1.)]

def split(trie, seq):
    oids = trie.items(str(seq))
    oid_count = 1./len(oids)
    return [(oid[1], oid_count) for oid in oids]

def rand(trie, seq):
    oids = trie.items(str(seq))
    return [(random.choice(oids)[1], 1.)]

def last(trie, seq):
    oids = trie.items(str(seq))
    return [(oids[-1][1], 1.)]

def all(trie, seq):
    oids = trie.items(str(seq))
    return [(oid[1], 1.) for oid in oids]

count_fs = {'exact': exact, 'split': split, 'rand': rand, 'last': last,
            'all': all}

def build_trie(seq_generator):
    sequence_summary = {'sample_ids': {},
                        'count': 0,
                        'lengths': []}

    def unpacker(seq):
        sequence_summary['count'] += 1
        sequence_summary['lengths'].append(len(seq))
        sample_id = seq.metadata['id'].split('_', maxsplit=1)[0]
        sample_ids = sequence_summary['sample_ids']
        if not sample_id in sample_ids:
            sample_ids[sample_id] = len(sample_ids)
        return str(seq)

    return marisa_trie.Trie(map(unpacker, seq_generator)), sequence_summary

def biom_from_trie(trie, sids, seq_iter, count_f=exact):
    """
    """
    # Build numpy array to store data. This will ultimately need to be a
    # sparse table, though the sparse DataFrames seem to be pretty slow to
    # populate.
    data = np.zeros((len(sids), len(trie)))
    # Create a result DataFrame. This will be our biom table.
    result = pd.DataFrame(data, columns=range(len(trie)), index=sids)

    for seq in seq_iter:
        sid = seq.metadata['id'].split('_', maxsplit=1)[0]
        for oid, count in count_f(trie, seq):
            result[oid][sid] += count
    # this requires two passes, there must be a better way to drop columns with zero count in place
    result.drop([i for i,e in enumerate(result.sum()) if e == 0], axis=1, inplace=True)
    return result

def pcoa_from_biom(biom, sample_md, color_by):

    title = "Samples colored by %s." % color_by

    dm = pw_distances(biom, biom.index)
    pcoa_results = PCoA(dm).scores()
    _ = pcoa_results.plot(df=sample_md,
                          column=color_by,
                          axis_labels=['PC 1', 'PC 2', 'PC 3'],
                          title=title,
                          s=35)


def table_summary(df):
    print("Samples: ", len(df.index))
    print("Observations: ", len(df.columns))
    print("Sequence/sample count detail:")
    print(df.T.sum().describe())

def load_table():
    return pd.DataFrame.from_csv('.table.biom')

def store_table(table):
    table_path = os.path.abspath('.table.biom')
    return table.to_csv(table_path)

def load_rarified_table():
    return pd.DataFrame.from_csv('.rarified-table.biom')

def store_rarified_table(table):
    table_path = os.path.abspath('.rarified-table.biom')
    return table.to_csv(table_path)

def get_markdown_template(fn):
    base_dir = os.path.abspath(os.path.split(__file__)[0])
    return open(os.path.join(base_dir, "markdown", fn)).read()

def get_seqs_to_biom_markdown(seqs_fp, count_f, command, output_fp):
    shutil.copy(seqs_fp, os.path.join(output_fp, '.seqs'))
    seqs_to_biom_md_template = get_markdown_template('seqs-to-biom.md')
    result = seqs_to_biom_md_template.format('.seqs', count_f, __version__, "dummy-md5", command, seqs_fp)
    return result

def get_biom_to_pcoa_markdown(map_fp, color_by, command, output_fp):
    shutil.copy(map_fp, os.path.join(output_fp, '.sample-md'))
    biom_to_pcoa_md_template = get_markdown_template('biom-to-pcoa.md')
    result = biom_to_pcoa_md_template.format('.sample-md', color_by, __version__, "dummy-md5", command, map_fp)
    return result

def get_biom_to_adiv_markdown(map_fp, collated_alpha_fp, command, output_fp):
    shutil.copy(map_fp, os.path.join(output_fp, '.sample-md'))
    #shutil.copy(collated_alpha_fp, os.path.join(output_fp, '.collated-alpha'))
    md_template = get_markdown_template('biom-to-adiv.md')
    print(collated_alpha_fp)
    result = md_template.format('.sample-md', collated_alpha_fp, __version__, "dummy-md5", command, map_fp)
    return result

def get_index_markdown(analysis_root):
    index_md_template = get_markdown_template('index.md')
    md_fps = glob.glob(os.path.join(analysis_root, '*.md'))
    md_fps.sort()
    toc = []
    for md_fp in md_fps:
        md_fn = os.path.split(md_fp)[1]
        toc.append(' * [%s](%s)' % (md_fn.split('.')[1].replace('-', ' ').title(), md_fn))
    toc = '\n'.join(toc)
    result = index_md_template.format(toc)
    return result

def _summarize_even_sampling_depth(even_sampling_depth, counts):
    samples_retained = (counts >= even_sampling_depth)
    num_samples_retained = samples_retained.sum()
    num_sequences_retained = num_samples_retained * even_sampling_depth
    return samples_retained, num_samples_retained, num_sequences_retained

def _get_depth_for_max_sequence_count(counts):
    """Find the even sampling depth that retains the most sequences."""
    count_summary = counts.describe()
    def f(d):
        return -1 * _summarize_even_sampling_depth(d, counts)[2]

    res = minimize_scalar(f,
                          bounds=(count_summary['min'], count_summary['max']),
                          method='bounded')
    return int(np.floor(res.x))

def get_default_even_sampling_depth(biom):
    counts = biom.T.sum()
    return _get_depth_for_max_sequence_count(counts)

def explore_sampling_depth(biom):
    import seaborn as sns
    counts = biom.T.sum()
    count_summary = counts.describe()
    total_num_samples = len(counts)
    total_num_sequences = counts.sum()
    depth_for_max_sequence_count = _get_depth_for_max_sequence_count(counts)
    sampling_depth_slider = IntSlider(min=count_summary['min'],
                                      max=count_summary['max'],
                                      step=10 ** (math.log(count_summary['max'], 10) - 2),
                                      value=depth_for_max_sequence_count)
    default_samples_retained, default_num_samples_retained, default_num_sequences_retained = \
            _summarize_even_sampling_depth(depth_for_max_sequence_count, counts)

    default_percent_samples_retained = default_num_samples_retained * 100 / total_num_samples
    default_percent_sequences_retained = default_num_sequences_retained * 100 / total_num_sequences

    label_s = "Depth {0}: {1:.2f}% of sequences and {2:.2f}% of samples retained."

    def f(even_sampling_depth):
        samples_retained, num_samples_retained, num_sequences_retained = \
            _summarize_even_sampling_depth(even_sampling_depth, counts)
        percent_samples_retained = num_samples_retained * 100 / total_num_samples
        percent_sequences_retained = num_sequences_retained * 100 / total_num_sequences
        ax = sns.distplot(counts)
        ax.set_xlabel("Number of sequences per sample")
        ax.set_ylabel("Frequency")
        line_label = label_s.format(depth_for_max_sequence_count,
                                    default_percent_sequences_retained,
                                    default_percent_samples_retained)
        ax.plot([depth_for_max_sequence_count, depth_for_max_sequence_count], ax.get_ylim(),
                'k--', label=line_label)

        line_label = label_s.format(even_sampling_depth,
                                    percent_sequences_retained,
                                    percent_samples_retained)
        ax.plot([even_sampling_depth, even_sampling_depth], ax.get_ylim(),
                'k-', label=line_label)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    def reset_depth(_):
        sampling_depth_slider.value = depth_for_max_sequence_count

    reset = widgets.Button(icon='fa-refresh')
    reset.on_click(reset_depth)

    w = interactive(f, even_sampling_depth=sampling_depth_slider)
    display(widgets.HBox(children=[w, reset]))

def rarify(biom, even_sampling_depth):
    data = []
    sample_ids = []
    for e in biom.index:
        count_vector = biom.loc[e]
        if count_vector.sum() < even_sampling_depth:
            continue
        else:
            sample_ids.append(e)
            data.append(subsample_counts(count_vector.astype(int), even_sampling_depth))
    return pd.DataFrame(data, index=sample_ids, columns=biom.columns)


markdown_templates = {'seqs-to-biom': get_seqs_to_biom_markdown,
                      'biom-to-pcoa': get_biom_to_pcoa_markdown,
                      'biom-to-adiv': get_biom_to_adiv_markdown,
                      'index': get_index_markdown}
