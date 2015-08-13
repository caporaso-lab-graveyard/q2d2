First, we'll load the table and the sample metadata.

```python
>>> %matplotlib inline
>>> import pandas as pd
>>> sample_md = pd.read_csv('{0}', sep='\t', index_col=0)
>>> alpha_diversity = pd.read_csv('{1}', sep='\t', index_col=1)
>>> # from q2d2 import load_rarified_table
>>> # table = load_rarified_table()
```

Some text...

```python
>>> from q2d2.wui import interactive_plot_alpha_diversity
>>> interactive_plot_alpha_diversity(sample_md, alpha_diversity)
```

Summary of what you did to prepare this analysis.
```
q2d2 version: {2}
collated alpha diversity data: {1}
map filepath: {5}
command: {4}
```
