First, we'll load the table and the sample metadata.

```python
>>> import pandas as pd
>>> sample_md = pd.read_csv('{0}', sep='\t', index_col=0)
>>> from q2d2 import load_rarified_table
>>> metric = 'braycurtis'
>>> category = '{1}'
>>> table = load_rarified_table()
```

Next, we compute Bray-Curtis distances, and then principal coordinates. This process
is described in detail in [IAB's *Studying Biological Diversity*](http://readiab.org/book/0.1.1/3/1#4).

```python
>>> %matplotlib inline
>>> from q2d2 import pcoa_from_biom
>>> pcoa_from_biom(table, sample_md, category, metric)
```

Now we'll visually explore the distance distributions themselves.

```python
>>> from skbio.diversity.beta import pw_distances
>>> from q2d2 import interactive_distance_histograms
>>> dm = pw_distances(metric=metric, counts=table, ids=table.index)
>>> interactive_distance_histograms(dm, sample_md)
```

Summary of what you did to prepare this analysis.
```
q2d2 version: {2}
map filepath: {5}
command: {4}
```
