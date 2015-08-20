First, we'll load the table and the sample metadata.

```python
>>> %matplotlib inline
>>> import pandas as pd
>>> sample_md = pd.read_csv('.sample-md', sep='\t', index_col=0)
>>> from q2d2 import load_rarefied_table
>>> table = load_rarefied_table()
>>> metric = 'observed_otus'
```

Some text...

```python
>>> from q2d2 import biom_to_adiv
>>> alpha_diversities = biom_to_adiv(metric, table)
```

Some text...

```python
>>> from q2d2.wui import interactive_plot_alpha_diversity
>>> interactive_plot_alpha_diversity(sample_md, alpha_diversities)
```

Summary of what you did to prepare this analysis.
```
q2d2 version: {2}
collated alpha diversity data: {1}
map filepath: {5}
command: {4}
```
