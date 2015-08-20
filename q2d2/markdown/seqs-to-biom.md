# Construction of a BIOM table.

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

Finally, we'll store the table so we can use it again another time.

```python
>>> from q2d2 import store_table
>>> store_table(biom)
```
