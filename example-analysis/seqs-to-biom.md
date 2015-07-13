

# Construction of a BIOM table.

Summary of what you did to prepare this analysis.
```
q2d2 version: 0.0.0-dev
input filepath: /Users/caporaso/Google Drive/Teaching/NAU/2015-spring/bio290/assignments/qiime_assignment_materials/forensic-seqs.fna
input filepath md5: dummy-md5
command: ./scripts/q2d2 seqs-to-biom example-analysis /Users/caporaso/Google Drive/Teaching/NAU/2015-spring/bio290/assignments/qiime_assignment_materials/forensic-seqs.fna exact
```

In this first step, we'll construct the [Trie](https://en.wikipedia.org/wiki/Trie)
data structure, which is used for grouping sequences into OTUs.

```python
>>> from q2d2 import build_trie
>>> import skbio
>>> seq_gen = skbio.io.read('/Users/caporaso/Google Drive/Teaching/NAU/2015-spring/bio290/assignments/qiime_assignment_materials/forensic-seqs.fna', format='fasta',)
>>> trie, sequence_summary = build_trie(seq_gen)
>>> sample_count = len(sequence_summary['sample_ids'])
>>> sequence_count = sequence_summary['count']
>>> print("%d sequences from %d samples were used to build the Trie." % (sequence_count, sample_count))
34214 sequences from 104 samples were used to build the Trie.
```

Next, we'll group our sequences, computing counts using the ``exact`` method that you passed in on the command line.

```python
>>> from q2d2 import biom_from_trie, count_fs
>>> seq_gen = skbio.io.read('/Users/caporaso/Google Drive/Teaching/NAU/2015-spring/bio290/assignments/qiime_assignment_materials/forensic-seqs.fna', format='fasta',)
>>> biom = biom_from_trie(trie, sequence_summary['sample_ids'], seq_gen, count_f=count_fs['exact'])
```

Finally, we'll get summary statistics on our table.

```python
>>> from q2d2 import table_summary
>>> table_summary(biom)
Samples:  104
Observations:  8998
Sequence/sample count detail:
count    104.000000
mean     328.980769
std      100.871348
min       77.000000
25%      287.000000
50%      324.500000
75%      388.250000
max      638.000000
dtype: float64
```

```python

```
