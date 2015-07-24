# Construction of a BIOM table.

Summary of what you did to prepare this analysis.
```
q2d2 version: {2}
input filepath: {0}
input filepath md5: {3}
command: {4}
```

In this first step, we'll construct the [Trie](https://en.wikipedia.org/wiki/Trie)
data structure, which is used for grouping sequences into OTUs.

```python
>>> from q2d2 import build_trie
>>> import skbio
>>> seq_gen = skbio.io.read('{0}', format='fasta',)
>>> trie, sequence_summary = build_trie(seq_gen)
>>> sample_count = len(sequence_summary['sample_ids'])
>>> sequence_count = sequence_summary['count']
>>> print("%d sequences from %d samples were used to build the Trie." % (sequence_count, sample_count))
```

Next, we'll group our sequences, computing counts using the ``{1}`` method that you passed in on the command line.

```python
>>> from q2d2 import biom_from_trie, count_fs
>>> seq_gen = skbio.io.read('{0}', format='fasta',)
>>> biom = biom_from_trie(trie, sequence_summary['sample_ids'], seq_gen, count_f=count_fs['{1}'])
```

Then, we'll get summary statistics on our table.

```python
>>> from q2d2 import table_summary
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
>>> rarified_biom = rarify(biom, even_sampling_depth=even_sampling_depth)
>>> print("Rarified table to", even_sampling_depth, "sequences per sample.")
```

Finally, we'll store the table so we can use it again another time.

```python
>>> from q2d2 import store_table, store_rarified_table
>>> store_table(biom)
>>> store_rarified_table(rarified_biom)
```
