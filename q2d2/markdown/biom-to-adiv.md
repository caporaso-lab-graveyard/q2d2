First, we'll load the table and the sample metadata.

```python
>>> %matplotlib inline
>>> import pandas as pd
>>> from q2d2 import load_sample_metadata
>>> sample_md = load_sample_metadata()
>>> from q2d2 import load_rarefied_table
>>> table = load_rarefied_table()
>>> from q2d2 import load_tree
>>> tree = load_tree()
```

Some text...

```python
>>> from q2d2 import compute_alphas
>>> alpha_diversities = compute_alphas(table, tree)
```

Some text...

```python
>>> from q2d2.wui import interactive_plot_alpha_diversity
>>> interactive_plot_alpha_diversity(sample_md, alpha_diversities)
```
