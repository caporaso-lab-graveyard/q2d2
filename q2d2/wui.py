from IPython.html import widgets
from IPython.display import clear_output
import pandas as pd

import matplotlib.pyplot as plt
import numpy as np


#This is the function that creates the interface
def metadata_controls(df, callback, extras=None):

    def filter_category(filters):
        newdf = df
        for filter_ in filters.children:
            key = filter_.children[0].value
            value = filter_.children[1].value
            newdf = newdf[newdf[key].isin(value)]
        return newdf

    def create_updater(group):
        def group_updater(_, value):
            group.options = [str(e) for e in df[value].dropna().unique().tolist()]
        return group_updater

    def create_category(filters):
        def category_widget():
            cat = widgets.Dropdown(options=df.columns.tolist())
            group = widgets.SelectMultiple(width=130)
            group_updater = create_updater(group)

            cat.on_trait_change(group_updater, 'value')
            group.on_displayed(lambda: group_updater(cat.value))

            category = widgets.VBox(children=[cat, group])
            filters.children += (category,)

        return category_widget

    def remove_category(filters):
        def remove_cat_widget():

            remove = filters.children[-1]
            filters.children = filters.children[:-1]
            remove.close()


        return remove_cat_widget


    def onclick(category, filters):
        extra_values = []
        if extras is not None:
            extra_values = [w.value for w in extras.children]

        clear_output()
        callback(category, filter_category(filters), *extra_values)

    filter_add = widgets.Button(description="Add Filter")
    filter_rem = widgets.Button(description="Remove Filter")
    go = widgets.Button(description="Go!")

    plt_cat = widgets.Dropdown(options=df.columns.tolist(), description='Plot')

    menu = widgets.HBox(children=[plt_cat, go, filter_add, filter_rem])
    filters = widgets.HBox(children=[])

    filter_add.on_click(lambda _: create_category(filters)())
    filter_rem.on_click(lambda _: remove_category(filters)())
    go.on_click(lambda _: onclick(plt_cat.value, filters))

    children = [menu, filters]
    if extras is not None:
        children.append(extras)
    return widgets.VBox(children=children)

#These are functions necessary for the alpha diversity plots
def get_df_intersection(df1, df2):
    intersect_ids = set.intersection(set(df1.index.tolist()), set(df2.index.tolist()))
    df1 = df1.loc[intersect_ids, ]
    df2 = df2.loc[intersect_ids, ]
    return df1, df2


def merge_metadata_alpha_div(metadata, alpha_div):
    metadata, alpha_div = get_df_intersection(metadata, alpha_div)
    #return pd.concat([metadata, alpha_div], axis=1)
    metadata['Alpha diversity'] = alpha_div
    return metadata


def plot_alpha(metadata, category, hue, metric):
    import seaborn as sns
    with plt.rc_context(dict(sns.axes_style("darkgrid"),
                             **sns.plotting_context("notebook", font_scale=2))):
        width = len(metadata[category].unique())
        plt.figure(figsize=(width*4, 8))
        fig = sns.boxplot(x=category, y='Alpha diversity',
                    data=metadata.sort(category), hue=hue, palette='cubehelix')
        fig.set_title(metric)

        
def plot_alpha_diversity(metadata, alpha_div, category, hue=None, metric=None):
    metadata_alpha_div = merge_metadata_alpha_div(metadata, alpha_div)
    plot_alpha(metadata_alpha_div, category, hue, metric)

    
def interactive_plot_alpha_diversity(metadata, alpha_divs):
    def on_update(category, metadata, Hue, check, metric):
        alpha_diversity = alpha_divs[metric]
        if not check:
            Hue = None
        plot_alpha_diversity(metadata, alpha_diversity, category, Hue, metric)

    check = widgets.Checkbox(Description='Plot Hue', Value=True)
    plt_hue = widgets.Dropdown(options=metadata.columns.tolist(), description='Hue')
    
    metric_but = widgets.Dropdown(options=list(alpha_divs.keys()), description='Metrics')
    extras = widgets.HBox(children=[plt_hue, check, metric_but])

    return metadata_controls(metadata, on_update, extras)

###########Functions for the taxa sumamry plots#####################

def get_taxa_counts(metadata_df, otu_df, category):
    cols = metadata_df[category].dropna().unique()
    indices = otu_df.index.tolist()
    taxa_counts_df = pd.DataFrame(index=indices, columns=cols)
    for col in cols:
        id_s = metadata_df[metadata_df[category] == col].index
        otu_sums = otu_df[id_s].sum(axis=1)
        taxa_counts_df.loc[otu_sums.index, col] = otu_sums
    return taxa_counts_df

def normalize(df):
    for col in df.columns:
        normalized_col = df[col]/df[col].sum()
        df[col] = normalized_col
    return df

#plotting functions

def plot_stacked_bar(df):
    df = df.sort(axis=1)
    import seaborn as sns
    with plt.rc_context(dict(sns.axes_style("darkgrid"),
                         **sns.plotting_context("notebook", font_scale=1.8))):
        f, ax = plt.subplots(1, figsize=(10, 10))
        x = list(range(len(df.columns)))
        bottom = np.array([0] * len(df.columns))
        cat_percents = []
        for id_ in df.index:
            color = '#' + ''.join(np.random.choice(list('ABCDEF123456789'), 6))
            ax.bar(x, df.loc[id_], color=color, bottom=bottom, align='center')
            bottom = df.loc[id_] + bottom
            cat_percents.append(''.join(["[{0:.2f}] ".format(x) for x in df.loc[id_].tolist()]))

        legend_labels = [' '.join(e) for e in zip(cat_percents, df.index.tolist())]

        ax.set_xticks(x)
        ax.set_xticklabels(df.columns.tolist())
        ax.set_ylim([0, 1])
        ax.legend(legend_labels, loc='center left', bbox_to_anchor=(1, 0.5))

def filter_by_dic(df, dic):
    for key, value in dic.items():
        df = df[df[key].isin(value)]
    return df


def plot_taxa_summary(otu_table, metadata, taxonomy, category, level='Phylum', min_percent=1):
        otu_counts = get_taxa_counts(metadata, otu_table, category)
        normalized_otus = normalize(otu_counts)
        normalized_otus = normalized_otus[normalized_otus.sum(axis=1) >= min_percent/100]
        normalized_otus = normalize(normalized_otus)

        normalized_taxa = pd.concat([normalized_otus, taxonomy], axis=1)
        normalized_taxa = normalized_taxa.groupby(level).sum()
        normalized_taxa = normalized_taxa.dropna()

        plot_stacked_bar(normalized_taxa)


def interactive_plot_taxa_summary(metadata, otu_df, taxa_df):
    shared_sample_ids = set.intersection(set(metadata.index.tolist()), set(otu_df.columns.tolist()))
    shared_otu_ids = set.intersection(set(taxa_df.index.tolist()), set(otu_df.index.tolist()))
    if len(shared_otu_ids) == 0:
        raise ValueError("There are no OTU ids in common between your OTU table and your "
                         "OTU metadata file.")
    metadata = metadata.loc[shared_sample_ids, ]
    otu_df = otu_df[list(shared_sample_ids)]
    for index, level in enumerate(['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']):
        taxa_df[level] = taxa_df['taxonomy'].apply(lambda x: ' '.join(x.split(' ')[:index + 1]))
    def on_update(category, metadata, level, min_percent):

        plot_taxa_summary(otu_df, metadata, taxa_df, category, level, min_percent=min_percent)

    plt_level = widgets.Dropdown(options=['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'],
                                 description='Level')
    
    min_percent = widgets.BoundedFloatText(width=40,
    value=1,
    min=0.0,
    max=10.0,
    description='Min Percent:')
    
    extras = widgets.HBox(children=[plt_level, min_percent])

    return metadata_controls(metadata, on_update, extras)