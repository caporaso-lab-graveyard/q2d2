# Construction of a BIOM table.

First, we'll get summary statistics on our table.

```python
>>> from q2d2 import table_summary, load_table
>>> biom = load_table()
>>> table_summary(biom)
```

We can next explore how choosing an even sampling depth (important for analyses that require rarefaction, such as computing community richness or between-sample distances) will affect the number of sequences as samples that are included in those analyses. In the Jupyter Notebook, this analysis will be interactive. The dashed line in the plot illustrates the even sampling depth at which the largest number of sequences are retained, but it is very likely that you want to explore the space around this value in the context of which of your samples will be retained and which will be discarded.

```python
>>> %matplotlib inline
>>> from q2d2 import explore_sampling_depth
>>> explore_sampling_depth(biom)
```

Next, we'll rarify the table to an even number of sequences per sample. If you want to overwrite the default even sampling depth, do that by setting the value of ``user_supplied_even_sampling_depth`` here.

```python
>>> from q2d2 import get_default_even_sampling_depth
>>> user_supplied_even_sampling_depth = None
>>> even_sampling_depth = user_supplied_even_sampling_depth or get_default_even_sampling_depth(biom)
```

```python
>>> from q2d2 import rarify
>>> rarefied_biom = rarify(biom, even_sampling_depth=even_sampling_depth)
>>> print("Rarefied table to", even_sampling_depth, "sequences per sample.")
```

Finally, we'll store the table so we can use it again another time.

```python
>>> from q2d2 import store_rarefied_table
>>> store_rarefied_table(rarefied_biom)
```

Summary of what you did to prepare this analysis.
```
q2d2 version: {0}
input filepath md5: {1}
command: {2}
```
