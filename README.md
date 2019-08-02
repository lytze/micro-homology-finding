## 1. Find MH pairs across given sequence(s)

#### Description

Micro-homology pairs (MH pairs) are short identical subsequences flanking a interval of sequnce. For example, in sequence `**ATCG**AGCGCACTTAG**ATCG**`, MH pair with identical sequence `ATCG` are flanking with a interval.

The indel size means the length from the first base of the first MH to the first base of the second MH. This indel size is the size of the indel event when the sequence "duplicated" or "collapsed" induced by the MH.

The first part is to survey a given sequence for all MH pair occurence within limititions. Here we want to find MH pairs with homologous k-mers' size not smaller than 4 bps, interspace greater than 3 pbs, and indel size not over 100 bps.

#### Source

	find_mh.c

#### Usage

	find_mh [output prefix] < [input fasta file]
	e.g:
	find_mh test_data/test < test_data/test.fa

changing the macros for parameter definition will change the limitition settings.

#### Output Format

The program outputs file(s) with name `[output prefix].number.sequence_name.out` contains three fields seperated by tabs. These fields are:

1. starting position of the first MH in a pair
2. starting position of the second MH in a pair
3. the size of the MH (k-mer's k)

## 2. Generate signatures for MTD reads detection

#### Description

The second part is to generate a list of signature sequence features for sequencing reads that carry a MH induced tandem duplication/collapse. Here we used the GNU's `gawk` for scripting.

#### Script

	generate_signatures.awk

#### Usage

Change the `#!` interpreter indicator on the first line in the script to your `gwak` path, and run `chmod u+x` on the script, then you can

	generate_signatures.awk <SZ=[feature size]> <TH=[feature diff threshold]> reference.fa mh.1.out <mh.2.out ...>

Or directly call `gawk` for

	gawk -f generate_signatures.awk <SZ=[feature size]> <TH=[feature diff threshold]> reference.fa mh.1.out <mh.2.out ...>


