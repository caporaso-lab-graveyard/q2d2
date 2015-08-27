First, we'll load the table and the sample metadata.

```python
>>> import pandas as pd
>>> from q2d2 import load_sample_metadata
>>> sample_md = load_sample_metadata()
>>> from q2d2 import load_rarefied_table
>>> table = load_rarefied_table()
>>> from q2d2 import load_tree
>>> tree = load_tree()
```

Next, we compute beta diversity distance matrices. This process
is described in detail in [IAB's *Studying Biological Diversity*](http://readiab.org/book/0.1.1/3/1#4). We can visualize the distance matrix directly as a heatmap...

```python
>>> %matplotlib inline
>>> from q2d2 import compute_distance_matrices
>>> dms = compute_distance_matrices(table, tree)
>>> _ = dms['unweighted_unifrac'].plot(cmap='Greens')
```

... or more usefully, we can apply principal coordinates analysis, and view the first three principal coordinates as a 3D scatter plot.

```python
>>> from q2d2 import interactive_plot_pcoa
>>> interactive_plot_pcoa(sample_md, dms)
```

Now we'll visually explore the distance distributions themselves.

```python
>>> from q2d2 import interactive_distance_violinplots
>>> interactive_distance_violinplots(dm, sample_md)
```
