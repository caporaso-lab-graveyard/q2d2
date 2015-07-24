First, we'll load the table and the sample metadata.

```python
>>> import pandas as pd
>>> sample_md = pd.read_csv('{0}', sep='\t', index_col=0)
>>> from q2d2 import load_rarified_table
>>> table = load_rarified_table()
```

Next, we compute Bray-Curtis distances, and then principal coordinates. This process
is described in detail in [IAB's *Studying Biological Diversity*](http://readiab.org/book/0.1.1/3/1#4).

```python
>>> %matplotlib inline
>>> from q2d2 import pcoa_from_biom
>>> pcoa_from_biom(table, sample_md, '{1}')
```
