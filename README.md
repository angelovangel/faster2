# faster2
Ultra-fast summary for fastx files

## Description
This is a stripped-down version of my other Rust program, [faster](), with less functionality but focused on speed. It is arguably the fastest tool that gives summary information about fastq files. The program uses the [`kseq`]() library for reading fastx files and outputs various summary metrics.

## Usage and benchmarks

```bash
# simple
faster2 -t file.fastq

# the query file can be fasta, fastq, fastq.gz...
```

The benchmarks were performed with sampled Nanopore ([Zymo mock community dataset](https://github.com/LomanLab/mockcommunity)) data using [hyperfine](https://github.com/sharkdp/hyperfine). `faster2` was compared to `faster`, `seqkit stats` and `seqstats`.

![img](img/zymo.png)
