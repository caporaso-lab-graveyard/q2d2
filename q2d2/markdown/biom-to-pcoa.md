First, we'll load the table and the sample metadata.

```python
>>> import pandas as pd
>>> sample_md = pd.read_csv('{0}', sep='\t', index_col=0)
>>> from q2d2 import load_rarefied_table
>>> metric = 'braycurtis'
>>> category = '{1}'
>>> table = load_rarefied_table()
```

Next, we compute Bray-Curtis distances. This process
is described in detail in [IAB's *Studying Biological Diversity*](http://readiab.org/book/0.1.1/3/1#4). We can visualize the distance matrix directly as a heatmap...

```python
>>> %matplotlib inline
>>> from q2d2 import biom_to_dm
>>> dm = biom_to_dm(metric, table)
>>> _ = dm.plot(cmap='Greens')
```

... or more usefully, we can apply principal coordinates analysis, and view the first three principal coordinates as a 3D scatter plot.

```python
>>> from q2d2 import dm_to_pcoa
>>> dm_to_pcoa(dm, sample_md, category)
```

Now we'll visually explore the distance distributions themselves.

```python
>>> from q2d2 import interactive_distance_histograms
>>> interactive_distance_histograms(dm, sample_md)
```

Summary of what you did to prepare this analysis.
```
q2d2 version: {2}
map filepath: {5}
command: {4}
```
