The data in this directory is derived from the following study:

Forensic identification using skin bacterial communities.
Fierer N, Lauber CL, Zhou N, McDonald D, Costello EK, Knight R.
PNAS. 2010 Apr 6;107(14):6477-81.
doi: 10.1073/pnas.1000162107.

The data in the q191 directory was generated with [QIIME](http://www.qiime.org) 1.9.1. The following commands were run:

```bash
pick_open_reference_otus.py -i forensic-seqs.fna -o or-otus-q191
core_diversity_analyses.py -i or-otus-q191/otu_table_mc2_w_tax_no_pynast_failures.biom -e 200 -o cd200_q191 -t or-otus-q191/rep_set.tre -m forensic-map.txt
```

Files were then shuffled around and renamed for simplicity. 
