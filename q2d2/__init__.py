#!/usr/bin/env python

__version__ = "0.0.0-dev"

import random
import io
import itertools
from collections import defaultdict, namedtuple
import hashlib
import os
import shutil
import glob
import math
from functools import partial

import marisa_trie
import numpy as np
import pandas as pd

from IPython.html import widgets
from IPython.html.widgets import interactive, fixed, IntSlider
from IPython.display import display
from scipy.optimize import minimize_scalar

import skbio
from skbio.diversity.beta import pw_distances
import skbio.diversity.alpha
from skbio.stats.ordination import PCoA
from skbio.stats import subsample_counts

from q2d2.wui import metadata_controls

WorkflowCategory = namedtuple('WorkflowCategory', ['name', 'title', 'workflows'])
Workflow = namedtuple('Workflow', ['name', 'title', 'inputs'])

workflows = [
    WorkflowCategory('no-biom', 'No BIOM table', []),
    WorkflowCategory('raw-biom', 'Raw (unnormalized) BIOM table', [
        Workflow('rarefy-biom', 'Rarefy BIOM table', ['unrarefied_biom']),
        Workflow('biom-to-taxa-plots', 'Taxonomy plots',
                 ['unrarefied_biom', 'sample_metadata', 'otu_metadata']),
    ]),
    WorkflowCategory('normalized-biom', 'Normalized BIOM table', [
        Workflow('biom-to-adiv', 'Alpha diversity', ['rarefied_biom', 'sample_metadata']),
        Workflow('biom-to-bdiv', 'Beta diversity', ['rarefied_biom', 'sample_metadata']),
    ])
]

def create_index(study_name, command):
    markdown_s = get_index_markdown(study_name, command)
    output_filepath = os.path.join(study_name, 'index.md')
    open(output_filepath, 'w').write(markdown_s)

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


type_to_study_filepath = {'sample_metadata': '.sample-md',
                          'otu_metadata': '.otu-md',
                          'unrarefied_biom': '.biom',
                          'rarefied_biom': '.rarefied-biom',
                          'tree': '.tree'}

def create_input_files(study_name, **kwargs):
    for input_type, input_filepath in kwargs.items():
        study_filepath = type_to_study_filepath[input_type]
        study_filepath = os.path.join(study_name, study_filepath)
        shutil.copy(input_filepath, study_filepath)

def load_table(rarefied=False):
    if rarefied:
        table_path = type_to_study_filepath['rarefied_biom']
    else:
        table_path = type_to_study_filepath['unrarefied_biom']
    result = pd.read_csv(table_path, sep='\t', skiprows=1, index_col=0)
    result.index = result.index.astype(str)
    if 'taxonomy' in result:
        result.drop('taxonomy', axis=1, inplace=True)
    return result

def store_table(table, rarefied=False):
    if rarefied:
        table_path = type_to_study_filepath['rarefied_biom']
    else:
        table_path = type_to_study_filepath['unrarefied_biom']
    with open(table_path, 'w') as table_file:
        table_file.write('# Constructed by [q2d2](github.com/gregcaporaso/q2d2)\n')
        table.to_csv(table_file, index_label="#OTU ID", sep='\t')

load_rarefied_table = partial(load_table, rarefied=True)
store_rarefied_table = partial(store_table, rarefied=True)

def load_sample_metadata():
    return pd.read_csv(type_to_study_filepath['sample_metadata'], sep='\t', index_col=0)

def load_otu_metadata():
    return pd.read_csv(type_to_study_filepath['otu_metadata'], sep='\t', names=['OTU ID', 'taxonomy'],
                       index_col=0, usecols=[0, 1], dtype=object)

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
    return result.T



def biom_to_adiv(metric, biom):
    metric_f = getattr(skbio.diversity.alpha, metric)
    results = []
    for e in biom.columns:
        results.append(metric_f(biom[e]))
    return pd.Series(results, index=biom.columns)

def biom_to_dm(metric, biom):
    return pw_distances(metric=metric, counts=biom.T, ids=biom.columns)

def dm_to_pcoa(dm, sample_md, category):
    title = "Samples colored by %s." % category
    pcoa_results = PCoA(dm).scores()
    _ = pcoa_results.plot(df=sample_md,
                          column=category,
                          axis_labels=['PC 1', 'PC 2', 'PC 3'],
                          title=title,
                          s=35)


def table_summary(df):
    print("Samples: ", len(df.columns))
    print("Observations: ", len(df.index))
    print("Sequence/sample count detail:")
    print(df.sum().describe())

def get_workflow_template_filepath(workflow_id):
    base_dir = os.path.abspath(os.path.split(__file__)[0])
    return os.path.join(base_dir, "markdown", "%s.md" % workflow_id)

def get_seqs_to_biom_markdown(seqs_fp, count_f, command, output_fp):
    shutil.copy(seqs_fp, os.path.join(output_fp, '.seqs'))
    seqs_to_biom_md_template = get_markdown_template('seqs-to-biom.md')
    result = seqs_to_biom_md_template.format('.seqs', count_f, __version__, "dummy-md5", command, seqs_fp)
    return result

def create_workflow(workflow_id, study_id):
    workflow_template_filepath = get_workflow_template_filepath(workflow_id)
    output_fn = os.path.split(workflow_template_filepath)[1]
    workflow_filepath = os.path.join(study_id, output_fn)
    if not os.path.exists(workflow_filepath):
        shutil.copy(workflow_template_filepath, workflow_filepath)
    return workflow_filepath

def get_index_markdown(study_name, command):
    index_md_template = open(get_workflow_template_filepath('index')).read()
    md_fps = glob.glob(os.path.join(study_name, '*.md'))
    md_fps.sort()
    toc = []
    for md_fp in md_fps:
        md_fn = os.path.split(md_fp)[1]
        title = os.path.splitext(md_fn)[0].replace('-', ' ').title()
        toc.append(' * [%s](%s)' % (title, md_fn))
    toc = '\n'.join(toc)
    result = index_md_template.format(toc, study_name, __version__, command)
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
    counts = biom.sum()
    return _get_depth_for_max_sequence_count(counts)

def explore_sampling_depth(biom):
    import seaborn as sns
    counts = biom.sum()
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
    for e in biom.columns:
        count_vector = biom[e]
        if count_vector.sum() < even_sampling_depth:
            continue
        else:
            sample_ids.append(e)
            data.append(subsample_counts(count_vector.astype(int), even_sampling_depth))
    return pd.DataFrame(np.asarray(data).T, index=biom.index, columns=sample_ids)


def filter_dm_and_map(dm, map_df):
    ids_to_exclude = set(dm.ids) - set(map_df.index.values)
    ids_to_keep = set(dm.ids) - ids_to_exclude
    filtered_dm = dm.filter(ids_to_keep)
    filtered_map = map_df.loc[ids_to_keep]

    return filtered_dm, filtered_map

def get_within_between_distances(map_df, dm, col):
    filtered_dm, filtered_map = filter_dm_and_map(dm, map_df)
    groups = []
    distances = []
    map_dict = filtered_map[col].to_dict()
    for id_1, id_2 in itertools.combinations(filtered_map.index.tolist(), 2):
        row = []
        if map_dict[id_1] == map_dict[id_2]:
            groups.append('Within')
        else:
            groups.append('Between')
        distances.append(filtered_dm[(id_1, id_2)])
    groups = zip(groups, distances)
    distances_df = pd.DataFrame(data=list(groups), columns=['Groups', 'Distance'])

    return distances_df

def distance_histogram(dm, category, metadata, metric='Distance', order=['Within', 'Between']):
    import seaborn as sns
    within_bw_distances = get_within_between_distances(metadata, dm, category)
    ax = sns.violinplot(x='Groups', y='Distance', data=within_bw_distances, order=order, orient='v')
    ax.set_xlabel(category)
    ax.set_ylabel(metric)

def interactive_distance_histograms(dm, sample_metadata):
    def on_update(category, metadata, check_within, check_between):
        order = []
        if check_within:
            order.append('Within')
        if check_between:
            order.append('Between')
        distance_histogram(dm, category, metadata, order=order)
    check_within = widgets.Checkbox(description='Show within category', value=True)
    check_between = widgets.Checkbox(description='Show between category', value=True)
    extras = widgets.VBox(children=[check_within, check_between])
    return metadata_controls(sample_metadata, on_update, extras)
