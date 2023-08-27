# sumi: A simple small RNA umi analysis
[![Rust](https://github.com/cookienocreams/sumi/actions/workflows/rust.yml/badge.svg)](https://github.com/cookienocreams/sumi/actions/workflows/rust.yml)
[![Continuous integration](https://github.com/cookienocreams/sumi/actions/workflows/CI.yaml/badge.svg)](https://github.com/cookienocreams/sumi/actions/workflows/CI.yaml)

A simple analysis for small RNA libraries with UMIs


It performs UMI error correction and deduplication using a directional graph algorithm. This script 
implements a slightly modified directional graph algorithm that allows for a Hamming Distance of 1 
between UMIs, see Fu, Y., et al, (2018). Elimination of PCR duplicates in RNA-seq and 
small RNA-seq using unique molecular identifiers. https://doi.org/10.1186/s12864-018-4933-1. 
It uses a five fold threshold, while the original algorithm uses a two fold count threshold.

## Script dependencies
- bowtie2
- samtools
- cutadapt

## Installation instructions

Create executable to run on local machine using the compiled binary:

You will need to have Rust installed on your computer before starting. Rust can be installed from here: [Install Rust](https://www.rust-lang.org/tools/install)

Download the git repository using git or manually and change into the repository folder.
```bash
git clone https://github.com/cookienocreams/sumi.git sumi
cd sumi
```

The next step is to install all libraries and their dependencies and compile the code into a binary.

```rust
cargo build --release
```

Before the binday can be used, it must be made executable. The compliled binary will be in `./target/release/`.

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
sumi --reference /home/user/data/miRNA --umi_regex "(.{5})CCA(.{2})TCA(.{5})" --threads 24
```
