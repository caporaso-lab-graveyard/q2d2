First, we'll load the table and the sample metadata.

```python
>>> import pandas as pd
>>> sample_md = pd.read_csv('{0}', sep='\t', index_col=0)
>>> from q2d2 import load_table
>>> table = load_table()
>>> taxa = pd.read_csv('{2}', sep='\t', names=['TaxID', 'Species'], index_col=0, usecols=[0, 1], dtype=object)
```

Text describing workflow...

```python
>>> %matplotlib inline
>>> from q2d2.wui import interactive_plot_taxa_summary
>>> interactive_plot_taxa_summary(sample_md, table, taxa, min_percent=1)
```

Summary of what you did to prepare this analysis.
```
q2d2 version: {3}
map filepath: {6}
tax filepath: {2}
biom filepath: {1}
command: {5}
```
