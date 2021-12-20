[![DOI](https://zenodo.org/badge/440225959.svg)](https://zenodo.org/badge/latestdoi/440225959)

# Description

This Perl 5 script was written to quickly scan GenBank files created by
antiSMASH v5 for non-ribosomal peptide synthase (NRPS) type biosynthetic gene
clusters (BGCs) containing at least one "condensation starter" class C-domain
(C-start). Such NRPS BGCs are predicted to produce lipopeptide products (NRLP).
Predicted NRLPs are screened for predicted cationic building blocks (cationic
NRLPs, CNRLPs).

Output is a tab-separate values (TSV) file with each line describing a
potential NRLP with the number of cationic building blocks and what they are
predicted to be. A gzip-compressed FASTA file is also written out on STDERR
containing sequencs for the C

# Usage

Screen one GenBank antiSMASH result file:

```bash
script/find_cnrp.pl bgc.gbk > out.tsv 2> out.fasta.gz
```

Logging is automatically written to a file in the current working directory
called `cnrp.log`.

Parallelize using the `find` command and GNU `parallel`:

```bash
# Using GNU parallel and find for a large number of files
find input/ -name '*gbk' -print0 |
    parallel -0 -j 64 \
    'find_cnrp.pl {} > output/{/.}.tsv 2 > output/{/.}.fasta.gz'
```

The output can then be combined using cat after parallel execution
is complete.

# Requirements

This script requires at least perl 5.26. It also requires the following perl modules:

- Bio::SeqIO (from [BioPerl](https://bioperl.org/))
- List::MoreUtils
- Log::Any

# Author

Yozen Hernandez <yzhernand at gmail dot com>

