#!/usr/bin/env python

__version__ = "0.0.0-dev"

import random
import io
from collections import defaultdict
import hashlib

import skbio
import marisa_trie
import numpy as np
import pandas as pd


def exact(trie, seq):
    return [(trie[str(seq)], 1.)]

def split(trie, seq):
    oids = trie.items(str(seq))
    oid_count = 1./len(oids)
    return [(oid[1], oid_count) for oid in oids]

def random(trie, seq):
    oids = trie.items(str(seq))
    return [(random.choice(oids)[1], 1.)]

def last(trie, seq):
    oids = trie.items(str(seq))
    return [(oids[-1][1], 1.)]

def all(trie, seq):
    oids = trie.items(str(seq))
    return [(oid[1], 1.) for oid in oids]

count_fs = {'exact': exact, 'split': split, 'random': random, 'last': last,
            'all': all}

def build_trie(seq_generator):
    sequence_summary = {'sample_ids': {},
                        'count': 0,
                        'lengths': []}

    def unpacker(seq):
        sequence_summary['count'] += 1
        sequence_summary['lengths'].append(len(seq))
        sample_id = seq.metadata['id'].split('_', maxsplit=1)[0]
        sample_ids = sequence_summary['sample_ids']
        if not sample_id in sample_ids:
            sample_ids[sample_id] = len(sample_ids)
        return str(seq)

    return marisa_trie.Trie(map(unpacker, seq_generator)), sequence_summary

def biom_from_trie(trie, sids, seq_iter, count_f=exact):
    """
    """
    # Build numpy array to store data. This will ultimately need to be a
    # sparse table.
    data = np.zeros((len(sids), len(trie)))
    # Create a result DataFrame. This will be our biom table.
    result = pd.DataFrame(data, columns=range(len(trie)), index=sids)

    for seq in seq_iter:
        sid = seq.metadata['id'].split('_', maxsplit=1)[0]
        for oid, count in count_f(trie, seq):
            result[oid][sid] += count
    # this requires two passes, there must be a better way to drop columns with zero count in place
    result.drop([i for i,e in enumerate(result.sum()) if e == 0], axis=1, inplace=True)
    return result

def table_summary(df):
    print("Samples: ", len(df.index))
    print("Observations: ", len(df.columns))
    print("Sequence/sample count detail:")
    print(df.T.sum().describe())

seqs_to_biom_md_template = """

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

Finally, we'll get summary statistics on our table.

```python
>>> from q2d2 import table_summary
>>> table_summary(biom)
```



"""

def get_seqs_to_biom_markdown(seqs_fp, count_f, command):

    return seqs_to_biom_md_template.format(seqs_fp, count_f, __version__, "dummy-md5", command)

markdown_templates = {'seqs-to-biom': get_seqs_to_biom_markdown}
