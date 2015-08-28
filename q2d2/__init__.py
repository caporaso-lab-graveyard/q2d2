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

import numpy as np
import pandas as pd

from IPython.html import widgets
from IPython.html.widgets import interactive, fixed, IntSlider
from IPython.display import display
from scipy.optimize import minimize_scalar

import skbio
from skbio.diversity.beta import pw_distances
import skbio.diversity.alpha
from skbio.stats.ordination import pcoa
from skbio.stats import subsample_counts
from skbio.util import safe_md5

from q2d2.wui import metadata_controls

data_type_to_study_filename = {'sample_metadata': '.sample-md',
                               'otu_metadata': '.otu-md',
                               'unrarefied_biom': '.biom',
                               'rarefied_biom': '.rarefied-biom',
                               'tree': '.tree'}

# this may make sense as a database schema. can we use an existing schema, e.g. Qiita?
WorkflowCategory = namedtuple('WorkflowCategory', ['title'])
Workflow = namedtuple('Workflow', ['title', 'inputs', 'outputs', 'category_id'])

workflow_categories = {
    'no-biom': WorkflowCategory('No BIOM table'),
    'raw-biom': WorkflowCategory('Raw (unnormalized) BIOM table'),
    'normalized-biom': WorkflowCategory('Normalized BIOM table')
}

workflows = {
    'rarefy-biom': Workflow(
        'Rarefy BIOM table', {'unrarefied_biom'}, {'rarefied_biom'}, 'raw-biom'),
    'biom-to-taxa-plots': Workflow(
        'Taxonomy plots', {'unrarefied_biom', 'sample_metadata', 'otu_metadata'}, {}, 'raw-biom'),
    'biom-to-adiv': Workflow(
        'Alpha diversity', {'rarefied_biom', 'sample_metadata'}, {}, 'normalized-biom'),
    'biom-to-bdiv': Workflow(
        'Beta diversity', {'rarefied_biom', 'sample_metadata'}, {}, 'normalized-biom')
}

def get_data_info(study_id):
    existing_data_types = get_existing_data_types(study_id)
    data_info = []
    for data_type in data_type_to_study_filename:
        filename = data_type_to_study_filename[data_type]
        exists = data_type in existing_data_types
        data_info.append((data_type, filename, exists))
    return data_info

def get_workflow_info(workflow_id):
    workflow = workflows[workflow_id]
    return {
        'workflow-id': workflow_id,
        'title': workflow.title,
        'inputs': list(workflow.inputs),
        'outputs': list(workflow.outputs),
        'category-id': workflow.category_id
    }

def get_workflow_category_info(category_id):
    return {
        'category-id': category_id,
        'title': workflow_categories[category_id].title
    }

def get_study_state(study_id):
    existing_data_types = get_existing_data_types(study_id)

    state = {
        'study-id': study_id,
        'workflow': {'exe': [], 'nexe': []},
        'data': {}
    }

    for workflow_id in workflows:
        workflow = workflows[workflow_id]
        if workflow.inputs.issubset(existing_data_types):
            state['workflow']['exe'].append(workflow_id)
        else:
            state['workflow']['nexe'].append(workflow_id)

    for data_type in existing_data_types:
        data_filepath = get_data_filepath(data_type, study_id)
        with open(data_filepath, 'rb') as data_file:
            # should we be using sha256 instead?
            md5 = safe_md5(data_file).hexdigest()
        state['data'][data_filepath] = md5

    return state

def get_system_info():
    # what other info goes here? dependencies?
    return {'version': __version__}

def get_existing_data_types(study_id):
    data_types = set()
    for data_type in data_type_to_study_filename:
        try:
            get_data_filepath(data_type, study_id)
        except FileNotFoundError:
            pass
        else:
            data_types.add(data_type)
    return data_types

def create_index(study_id, command):
    markdown_s = get_index_markdown(study_id, command)
    output_filepath = os.path.join(study_id, 'index.md')
    open(output_filepath, 'w').write(markdown_s)

def get_data_filepath(data_type, study_id):
    data_filepath = os.path.join(study_id, data_type_to_study_filename[data_type])
    if not os.path.exists(data_filepath):
        raise FileNotFoundError(data_filepath)
    return data_filepath

def create_input_files(study_id, **kwargs):
    for input_type, input_filepath in kwargs.items():
        study_filepath = data_type_to_study_filename[input_type]
        study_filepath = os.path.join(study_id, study_filepath)
        shutil.copy(input_filepath, study_filepath)

def load_table(rarefied=False):
    if rarefied:
        table_path = data_type_to_study_filename['rarefied_biom']
    else:
        table_path = data_type_to_study_filename['unrarefied_biom']
    result = pd.read_csv(table_path, sep='\t', skiprows=1, index_col=0)
    result.index = result.index.astype(str)
    if 'taxonomy' in result:
        result.drop('taxonomy', axis=1, inplace=True)
    return result

def store_table(table, rarefied=False):
    if rarefied:
        table_path = data_type_to_study_filename['rarefied_biom']
    else:
        table_path = data_type_to_study_filename['unrarefied_biom']
    with open(table_path, 'w') as table_file:
        table_file.write('# Constructed by [q2d2](github.com/gregcaporaso/q2d2)\n')
        table.to_csv(table_file, index_label="#OTU ID", sep='\t')

load_rarefied_table = partial(load_table, rarefied=True)
store_rarefied_table = partial(store_table, rarefied=True)

def load_tree():
    return skbio.TreeNode.read(data_type_to_study_filename['tree'], format='newick')

def load_sample_metadata():
    return pd.read_csv(data_type_to_study_filename['sample_metadata'], sep='\t', index_col=0)

def load_otu_metadata():
    return pd.read_csv(data_type_to_study_filename['otu_metadata'], sep='\t', names=['OTU ID', 'taxonomy'],
                       index_col=0, usecols=[0, 1], dtype=object)

def biom_to_adiv(metric, biom, tree=None):
    metric_f = getattr(skbio.diversity.alpha, metric)
    results = []
    for e in biom.columns:
        if metric == 'faith_pd':
            results.append(metric_f(biom[e], biom.index, tree))
        else:
            results.append(metric_f(biom[e]))
    return pd.Series(results, index=biom.columns)

def compute_alphas(otu_table, tree=None,
                   metrics=['chao1',
                            'faith_pd',
                            'observed_otus']):
    alphas = {}
    for metric in metrics:
        alpha = biom_to_adiv(metric, otu_table, tree)
        alphas[metric] = alpha 
    return alphas

def biom_to_dm(metric, biom, tree=None):
    return pw_distances(metric=metric, counts=biom.T, ids=biom.columns)

def dm_to_pcoa(dm, sample_md, category):
    title = "Samples colored by %s." % category
    pcoa_results = pcoa(dm)
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

def get_workflow_filepath(workflow_id, study_id):
    return os.path.join(study_id, "%s.md" % workflow_id)

def create_workflow(workflow_id, study_id):
    workflow_template_filepath = get_workflow_template_filepath(workflow_id)
    workflow_filepath = get_workflow_filepath(workflow_id, study_id)
    if not os.path.exists(workflow_filepath):
        shutil.copy(workflow_template_filepath, workflow_filepath)
    return workflow_filepath

def delete_workflow(workflow_id, study_id):
    workflow_filepath = get_workflow_filepath(workflow_id, study_id)
    os.remove(workflow_filepath)

def get_index_markdown(study_id, command):
    index_md_template = open(get_workflow_template_filepath('index')).read()
    md_fps = glob.glob(os.path.join(study_id, '*.md'))
    md_fps.sort()
    toc = []
    for md_fp in md_fps:
        md_fn = os.path.split(md_fp)[1]
        title = os.path.splitext(md_fn)[0].replace('-', ' ').title()
        toc.append(' * [%s](%s)' % (title, md_fn))
    toc = '\n'.join(toc)
    result = index_md_template.format(toc, study_id, __version__, command)
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

def distance_violinplots(dm, category, metadata, metric=None, order=['Within', 'Between']):
    import seaborn as sns
    within_bw_distances = get_within_between_distances(metadata, dm, category)
    ax = sns.violinplot(x='Groups', y='Distance', data=within_bw_distances, order=order, orient='v')
    ax.set_xlabel(category)
    ax.set_ylabel(metric)
    return ax

def interactive_distance_violinplots(dms, sample_metadata):
    
    def on_update(category, metadata, metric, check_within, check_between):
        order = []
        if check_within:
            order.append('Within')
        if check_between:
            order.append('Between')
        
        dm = dms[metric]
        distance_violinplots(dm, category, metadata, metric, order=order)
        
    check_within = widgets.Checkbox(description='Show within category', value=True)
    check_between = widgets.Checkbox(description='Show between category', value=True)
    metric_but = widgets.Dropdown(options=list(dms.keys()), description='Metrics')

    
    extras = widgets.VBox(children=[metric_but, check_within, check_between])
    return metadata_controls(sample_metadata, on_update, extras)

def compute_distance_matrices(
               otu_table,
               tree=None,
               metrics=['weighted_unifrac', 'unweighted_unifrac', 'braycurtis', 'jaccard']):
    dms = {}
    for metric in metrics:
        dm = pw_distances(metric, otu_table.T.values, otu_table.columns.tolist(), 
                             tree=tree, otu_ids=otu_table.index.tolist())
        dms[metric] = dm
    return dms

def interactive_plot_pcoa(metadata, dms):

    def on_update(category, metadata, metric):
        dm = dms[metric]
        filtered_dm, _ = filter_dm_and_map(dm, metadata)
        pc = pcoa(filtered_dm)
        pc.plot(df=metadata,
        column=category,
        axis_labels=['PC 1', 'PC 2', 'PC 3'],
        s=35).set_size_inches(12, 9)
        
    metric_but = widgets.Dropdown(options=list(dms.keys()), description='Metrics')
    extras = widgets.VBox(children=[metric_but])
    
    return metadata_controls(metadata, on_update, extras)