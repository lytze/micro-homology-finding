### Find Micro-homology Pairs Across Give DNA Sequence

#### Definition: Micro-homology Pair

Micro-homology pairs (MH pairs) are short identical subsequences flanking a interval of sequnce.

For example, in sequence `**ATCG**AGCGCACTTAG**ATCG**`, MH pair with identical sequence `ATCG` are flanking with a interval.

The indel size means the length from the first base of the first MH to the first base of the second MH. This indel size is the size of the indel event when the sequence "duplicated" or "collapsed" induced by the MH.


#### Usage

For the `C` version:

	find_mh output_prefix < input.fa

An old `R` version is also provided.
