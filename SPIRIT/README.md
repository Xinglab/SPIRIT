# **SPIRIT** (**SP**lint **I**mproved **R**epeat **I**dentifier for **T**ranscripts)

SPIRIT is a computational pipeline for obtaining consensus sequences from long concatemeric cDNA reads with known splint sequences.

Please contact [Yuan Gao](mailto:gy.james@outlook.com) for questions, feedback, or bugs. 

## Table of Contents
* [Dependencies](#dependencies)
* [Usage](#usage)
* [Example](#example)

## Dependencies

SPIRIT requires the following to be installed and available on `$PATH`:

* [Perl](https://www.perl.org/) >= 5.8 built with threading enabled
  + Check for thread support with `perl -e 'use threads; print("ok\n")'`
* [hmmer](http://www.hmmer.org/) >= 3.3.1
  + T. J. Wheeler, S. R. Eddy, nhmmer: DNA homology search with profile HMMs. *Bioinformatics* **29**, 2487-2489 (2013).
* [abPOA](https://github.com/Xinglab/abPOA/) >= 1.0.6
  + Y. Gao, Y. Liu, Y. Ma, B. Liu, Y. Wang, Y. Xing, abPOA: an SIMD-based C library for fast partial order alignment using adaptive band. *Bioinformatics* **37**, 2209-2211 (2021).

## Usage

Prior to running SPIRIT, we need to first build a HMMER binary database file from a file containing our splint sequence. This can be done using the following command:

```
makehmmerdb [/path/to/splint/sequence/file] [/path/to/HMMER/binary/database/prefix]
```

After building this database file, we can run SPIRIT. The basic usage of SPIRIT is:

```
bash spirit.sh -I [/path/to/directory/of/input/read/fastq] \
    -O [/path/to/output/directory] \
    -D [/path/to/HMMER/binary/database/prefix] \
    -L [UMI length] \
    -S [/path/to/SPIRIT/directory] \
    -T [number of threads]
```

## Example

A small test dataset is provided in the `toy_data` folder. This folder contains the following files:

* `fastq_runid_head400_*_0.fastq.gz` - A gzip compressed FASTQ file containing long concatemeric cDNA reads (**note:** decompress these files prior to running SPIRIT)
* database/
  + `splint_50nt.fa` - A FASTA file containing the splint sequence
  + `splint_50nt.fa.*` - A collection of HMMER binary database files prepared from `splint_50nt.fa`

SPIRIT can be run on files from the `toy_data` folder as follows, assuming that 5 worker threads are available and we have a UMI length of 60:

```
bash spirit.sh -I [/path/to/toy_data] \
    -O "out" \
    -D [/path/to/toy_data/splint_50nt.fa] \
    -L "60" -S ./ -T "5"
```

Within a few minutes, SPIRIT should finish processing the files from the `toy_data` folder and generate the following output files in the `out` folder:

* `extract_unit_splints_consensus_versatile_out.fa` - A FASTA file containing the consensus sequences for each long concatemeric cDNA read
* `extract_unit_splints_consensus_versatile_out.fa.cDNA` - A FASTA file containing both UMI sequences and long concatemeric cDNA sequences
