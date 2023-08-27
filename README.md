# sumi: A simple small RNA umi analysis
A simple analysis for small RNA libraries with UMIs

## Installation instructions

Create executable to run on local machine using the compiled binary:

You will need to have Rust installed on your computer before starting. Rust can be installed from here: https://julialang.org/downloads/

Download the git repository using git or manually and change into the repository folder.
```bash
git clone https://github.com/cookienocreams/sumi.git sumi
cd sumi
```

The next step is to install all libraries and their dependencies and compile the code into a binary.

```rust
cargo build --release
```

The compliled binary will be in `./target/release/`. There are numerous options that can be changed 
if desired. Use `-h` or `--help` flags to see options.

```bash
./target/release/sumi --help
```

## Basic usage

The app can be run using the `sumi` executable. To this command to analyze files with a miRNA bowtie2 reference located
in `/home/user/data/`, a 12 bp UMI on the 5' end with the structure "NNNNNCCANNTCANNNNN", with 24 threads.

```bash
cd fastqs
/home/user/sumi/target/release/sumi --reference /home/user/data/miRNA --umi_regex "(.{5})CCA(.{2})TCA(.{5})" --threads 24
```
