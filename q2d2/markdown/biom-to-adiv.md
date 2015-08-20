First, we'll load the table and the sample metadata.

```python
>>> %matplotlib inline
>>> import pandas as pd
>>> from q2d2 import load_sample_metadata
>>> sample_md = load_sample_metadata()
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
