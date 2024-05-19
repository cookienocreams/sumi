# sumi: A simple small RNA umi analysis
[![Rust](https://github.com/cookienocreams/sumi/actions/workflows/rust.yml/badge.svg)](https://github.com/cookienocreams/sumi/actions/workflows/rust.yml)
[![Continuous integration](https://github.com/cookienocreams/sumi/actions/workflows/CI.yaml/badge.svg)](https://github.com/cookienocreams/sumi/actions/workflows/CI.yaml)
[![pages-build-deployment](https://github.com/cookienocreams/sumi/actions/workflows/pages/pages-build-deployment/badge.svg)](https://github.com/cookienocreams/sumi/actions/workflows/pages/pages-build-deployment)

A simple analysis for small RNA libraries with UMIs

Can be used to check for the presence of isomiRs in addition to canonical miRNAs.

It performs UMI error correction and deduplication using a directional graph algorithm. This script 
implements a slightly modified directional graph algorithm that allows for a Hamming Distance of 1 
between UMIs, see Fu, Y., et al, (2018). Elimination of PCR duplicates in RNA-seq and 
small RNA-seq using unique molecular identifiers. https://doi.org/10.1186/s12864-018-4933-1. 
It uses a five fold threshold, while the original algorithm uses a two fold count threshold.

The UMI deduplication process is multithreaded to significantly speed up the process.

## Script dependencies
- bowtie2
- samtools
- cutadapt

## Installation instructions

Create executable to run on local machine using the compiled binary.


You will need to have Rust installed on your computer before starting. Rust can be installed from here: [Install Rust](https://www.rust-lang.org/tools/install)

Download the git repository using cargo.
```bash
cargo install sumi
```

The compiled binary will be located in `~/.cargo/bin/`. Make it executable with the following command:
```bash
chmod +x ./target/release/sumi 
```

There are numerous options that can be changed if desired. Use `-h` or `--help` flags to see options.

```bash
./target/release/sumi --help
```

## Basic usage

The app can be run using the `sumi` executable. To this command to analyze files with a miRNA bowtie2 reference located
in `/home/user/data/`, a 12 bp UMI on the 5' end with the structure "NNNNNCCANNTCANNNNN", and with 8 threads.

```bash
cd fastqs
./sumi --reference /home/user/data/miRNA --umi-regex "^(.{5})CCA(.{2})TCA(.{5})" --threads 8
```

To check for isomiRs and generate counts and the alignment information, run the following command. 
Note that the isomiR and canonical miRNA percentages are separated.

```bash
./sumi --reference /home/user/data/miRNA --isomir --write-metrics
```

This can be used to analyze libraries with a 12 bp 3' UMI with a 1 bp mismatch allowed during alignment.
Make sure the regex pattern contains an anchor, such as `'$'`, so that it is specific to the 3' end.

```bash
./sumi --reference /home/user/data/miRNA --3p --umi-regex "(.{12}$)" --mismatch
```

For analyzing isomiRs in Qiagen libraries use this command, there is no need to specify the UMI is on the 3' end or regex pattern.

```bash
./sumi --reference /home/user/data/miRNA --isomir --qiagen
```
