First, we'll load the table and the sample metadata.

```python
>>> import pandas as pd
>>> from q2d2 import load_sample_metadata
>>> sample_md = load_sample_metadata()
>>> from q2d2 import load_table
>>> table = load_table()
>>> from q2d2 import load_otu_metadata
>>> taxa = load_otu_metadata()
```

Text describing workflow...

```python
>>> %matplotlib inline
>>> from q2d2.wui import interactive_plot_taxa_summary
>>> interactive_plot_taxa_summary(sample_md, table, taxa)
```
