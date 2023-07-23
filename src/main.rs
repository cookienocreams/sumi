mod sample;
mod trim_adapters;
mod extract_umis;
mod graph;
mod rna_counts;
mod thresholds;
mod lengths;
mod common;

use crate::sample::subsample_fastqs;
use crate::trim_adapters::trim_adapters;
use crate::extract_umis::extract_umis;
use crate::extract_umis::extract_umis_qiagen;
use crate::rna_counts::rna_discovery_calculation;
use crate::lengths::calculate_read_length_distribution;
use crate::thresholds::threshold_count;
use crate::thresholds::combine_threshold_counts;
use crate::common::find_common_rnas;
use crate::common::write_common_rna_file;

extern crate bam;
// Import the `lazy_static` macro
#[macro_use]
extern crate lazy_static;

use petgraph::graph::NodeIndex;
use petgraph::Graph;
use std::collections::{HashMap, HashSet};
use std::str;
use regex::Regex;
use std::process::Command;
use std::path::Path;
use std::fs::{self, File};
use thiserror::Error;
use clap::{App, Arg};
use polars::prelude::*;
use std::io::{self, BufReader, BufWriter, Read};
use indicatif::{ProgressBar, ProgressStyle};
use std::error::Error;
use niffler::get_reader;
use csv::{Reader as csv_reader, Writer as csv_writer};
use bio::alignment::distance::simd::hamming;
use bio::io::fastq::Reader;
use flate2::bufread::MultiGzDecoder;

// Define the regex to extract index information and UMI and index information
lazy_static! {
    static ref SINGLE_INDEX_REGEX: Regex = Regex::new(r"1:N:\d:[ATCGN]+").unwrap();
    static ref DUAL_INDEX_REGEX: Regex = Regex::new(r"1:N:\d:[ATCGN]+\+[ATCGN]+").unwrap();
    static ref READ_NAME_UMI_REGEX: Regex = Regex::new(r"_([ATCG]+)").unwrap();
    static ref UMI_REGEX_QIAGEN: Regex = Regex::new(r"AACTGTAGGCACCATCAAT([ATCG]{12})AGATCGGAAG").unwrap();
}

/// List all files in the current directory that contain the target string.
///
/// # Arguments
///
/// * `files_to_capture` - The string to look for in file names.
///
/// # Returns
///
/// A vector of strings, each string being the path of a file that contains the target string.
///
/// # Example
///
/// ```
/// let target_files = capture_target_files(".txt");
/// // target_files is ["file1.txt", "file2.txt", "file3.txt"]
/// ```
fn capture_target_files(files_to_capture: &str) -> Vec<String> {
    let entries = fs::read_dir(".").expect("Failed to read directory");

    let mut files = Vec::new();

    for entry in entries {
        if let Ok(entry) = entry {
            let path = entry.path();
            if let Some(filename) = path.file_name() {
                if let Some(filename) = filename.to_str() {
                    if filename.contains(files_to_capture) {
                        files.push(filename.to_string());
                    }
                }
            }
        }
    }

    files.sort_unstable();

    files
}

/// Check if a file is gzipped by reading the first two bytes and comparing them to the gzip magic numbers
fn is_gzipped(filename: &str) -> io::Result<bool> {
    let mut file = File::open(filename)?;
    let mut buffer = [0; 2];  // Buffer to store the first two bytes
    file.read(&mut buffer)?;
    Ok(buffer == [0x1f, 0x8b])
}

/// Calcualte the mean of an input vector
fn mean(numbers: &Vec<f64>) -> f64 {
    let sum: f64 = numbers.iter().sum();
    sum as f64 / numbers.len() as f64
}

/// Calculate the average quality score of a FASTQ file.
///
/// # Arguments
///
/// * `input_fastq` - A string slice that holds the name of the FASTQ file to read.
///
/// # Returns
///
/// This function will return the average quality score of the input FASTQ file. If an
/// error occurs while reading the file, an error will be returned instead.
///
/// # Examples
///
/// ```
/// let avg_quality = average_read_quality("sample.fastq")?;
/// println!("The average quality score is {}.", avg_quality);
/// ```
fn average_read_quality(input_fastq: &str) -> Result<f64, Box<dyn std::error::Error>> {
    let file = File::open(input_fastq)?;
    let reader = Reader::new(MultiGzDecoder::new(BufReader::new(file)));

    let mut q_score_list = Vec::new();

    for record_result in reader.records() {
        let record = record_result.expect("Error during fastq record parsing");
        if record.is_empty() {
            let _check = record.check();
            break;
        }
    
        let quals = record.qual();
        let sum_qual: f64 = quals.iter().map(|&c| c as f64 - 33.0).sum();

        let q_score = sum_qual / record.seq().len() as f64;

        q_score_list.push(q_score.round());
    }

    Ok(mean(&q_score_list))
}

/// Calculate the average quality scores for a list of FASTQ files.
///
/// # Arguments
///
/// * `fastq_files` - A slice containing the names of the FASTQ files to process.
///
/// # Returns
///
/// This function will return a vector of the average quality scores of the input FASTQ 
/// files. If an error occurs while reading any of the files, an error will be returned instead.
///
/// # Examples
///
/// ```
/// let avg_qualities = average_quality_scores(&["file1.fastq", "file2.fastq", "file3.fastq"])?;
/// println!("The average quality scores are: {:?}", avg_qualities);
/// ```
fn average_quality_scores(fastq_files: &Vec<String>) -> Result<(), Box<dyn std::error::Error>> {
    let num_of_fastqs = fastq_files.len() as u64;
    let progress_br = ProgressBar::new(num_of_fastqs);
    
    progress_br.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:50.cyan/blue}] {pos}/{len} {msg} ({percent}%)")
                .expect("Progress bar error")
            .progress_chars("#>-"),
    );
    progress_br.set_message("Calculating quality scores...");

    let mut quality_map = HashMap::new();
    for file in fastq_files {
        let avg_quality = average_read_quality(file)?;
        quality_map.insert(file.clone(), avg_quality);

        progress_br.inc(1);
    }

    let mut writer = csv_writer::from_path("q_scores.csv")?;
    writer.write_record(&["sample", "quality score"])?;

    for (sample, score) in &quality_map {
        writer.write_record(&[sample, &score.to_string()])?;
    }

    writer.flush()?;

    progress_br.finish_with_message("Finished calculating quality scores");

    Ok(())
}

/// Remove all intermediate files.
///
/// # Returns
///
/// Nothing, the function deletes files.
///
/// # Example
///
/// ```
/// remove_intermediate_files();
/// ```
fn remove_intermediate_files() {
    let files_to_delete: HashSet<String> = vec![
        ".bam",
        ".cut.fastq",
        "read_lengths_",
        ".miRNA.",
        ".cutadapt_information",
        ".sam",
        ".unprocessed.cut",
        "short.fastq",
        "_subsample",
        "processed",
        "_threshold_count_sum.csv",
        "_common_RNAs"
    ]
    .into_iter()
    .flat_map(capture_target_files)
    .collect();

    for file in files_to_delete {
        let _ = fs::remove_file(&file);
    }
}

fn main() -> io::Result<()> {
    let matches = App::new("Small RNA Analysis")
        .version("1.0")
        .author("Author: Michael Hawkins")
        .about("This script can be used to analyze Small RNA UMI libraries. \
        \n\nIt performs UMI error correction and deduplication using a directional graph algorithm.\
        This script implements a slightly modified directional graph algorithm that allows for a Hamming \
        Distance of 1 between UMIs, see Fu, Y., et al, (2018). Elimination of PCR duplicates in RNA-seq and \
        small RNA-seq using unique molecular identifiers. https://doi.org/10.1186/s12864-018-4933-1. \
        It uses a five fold threshold, while the original algorithm uses a two fold count threshold.

        \nDependencies: cutadapt version >= 4, samtools version >= 1.10.1, bowtie2 version >= 2.2.1"
    )
        .arg(
            Arg::with_name("minimum_length")
                .short('m')
                .long("min-length")
                .help("The miniumum length cutoff for cutadapt. Fragments shorter than this length \
                will be discarded after adapter trimming.")
                .value_name("MINIMUM_LENGTH")
                .value_parser(clap::value_parser!(u8).range(0..75))
                .default_value("16"),
        )
        .arg(
            Arg::with_name("num_threads")
                .short('t')
                .long("threads")
                .help("The number of processors to use for alignment.")
                .value_name("NUM_THREADS")
                .value_parser(clap::value_parser!(u8))
                .default_value("12"),
        )
        .arg(
            Arg::with_name("keep_intermediate_files")
                .short('k')
                .long("keep")
                .takes_value(false)
                .help("Flag to keep intermediate files if desired. Default is to remove intermediate files."),
        )
        .arg(
            Arg::with_name("alignment_reference")
                .short('r')
                .long("reference")
                .value_name("ALIGNMENT_REFERENCE")
                .help("The bowtie2 reference location. Use full path, e.g., /home/path/to/ref. \
                Name of reference for the target RNA type should be a simple one word name such \
                as 'miRNA' or 'pfeRNA'.")
                .takes_value(true)
                .default_value("./data/miRNA"),
        )
        .arg(
            Arg::with_name("subsample_fastqs")
                .short('s')
                .long("sample")
                .value_name("NUM_READS")
                .takes_value(true)
                .help("Subsample all fastqs in the current directory to a specified number of reads.")
        )
        .arg(
            Arg::with_name("subsample_to_lowest")
                .short('l')
                .long("sample-to-lowest")
                .takes_value(false)
                .help("Subsample all fastqs in the current directory to the fastq with \
                the lowest number of reads.")
        )
        .arg(
            Arg::with_name("qiagen")
                .short('q')
                .long("qiagen")
                .takes_value(false)
                .help("Set flag if Qiagen libraries are being analyzed.")
        )
        .arg(
            Arg::with_name("quality_score")
                .short('Q')
                .long("quality_score")
                .takes_value(false)
                .help("Set flag to write output file containing each fastq's average read quality.")
        )
        .arg(
            Arg::with_name("thresholds")
                .short('T')
                .long("thresholds")
                .multiple(true) 
                .help("The thresholds for RNA species counting. Multiple values must be \
                separated by a space. Default threshold values set to 1 3 5 10.")
                .value_name("THRESHOLDS"),
        )
        .arg(
            Arg::with_name("include_unaligned_reads")
                .short('u')
                .long("unaligned")
                .takes_value(false)
                .help("Include unaligned reads in sam file post bowtie2 alignment. Can be \
                used to obtain read lengths from all fragments or for further alignment \
                 outside the program if used with the '--keep' flag. Note this can significantly \
                 increase run time.")
        )
        .arg(
            Arg::with_name("umi_regex")
            .short('p')
            .long("umi_regex")
            .value_name("UMI_REGEX")
            .help("The regular expression pattern to capture each UMI. This pattern should contain capture \
                groups for the UMI bases and can include any intermediate bases.\n\n\
                For example, to capture a UMI of 10 bases, use \"(^.{10})\". If the UMI is split by specific \
                bases, include those bases in the pattern as well. For example, if a 12-base UMI is split \
                into three groups of 4 bases each by the bases 'CCA' and 'TCA', use \"(^.{4})CCA(.{4})TCA(.{4})\".\n\n\
                Note: This is currently limited to UMIs on the 5' side of each sequence \
                unless analyzing Qiagen libraries, see '--qiagen' flag, and the double \
                quotes are required.")
            .takes_value(true)
            .default_value("(^.{12})"),
        )
        .get_matches();

    // Check if Qiagen libraries are being used
    let is_qiagen = matches.is_present("qiagen");

    // Set the library type
    let library_type: &str;
    if is_qiagen {
        library_type = "Qiagen"
    } else {
        library_type = "UMI"
    }

    // Gather all the input parameters
    let keep_intermediates: bool = matches.is_present("keep_intermediate_files");
    let num_threads: u8 = *matches.get_one("num_threads").unwrap();
    let reference: &str = matches.value_of("alignment_reference").unwrap();
    let reference_name: &str = Path::new(reference).file_name().unwrap().to_str().unwrap();
    let include_unaligned_reads = matches.is_present("include_unaligned_reads");
    let input_regex = matches.value_of("umi_regex").unwrap();

    let umi_regex_result = Regex::new(matches.value_of("umi_regex").unwrap());
    let umi_regex = match umi_regex_result {
        Ok(regex) => regex,
        Err(_) => panic!("Unexpected umi regex pattern")
    };
    
    let thresholds: Vec<usize> = match matches.values_of("thresholds") {
        Some(values) => values.map(|x| x.parse::<usize>().expect("Threshold must be a number.")).collect(),
        None => vec![1, 3, 5, 10],
    };

    let subsample: Option<u32> = match matches.value_of("subsample_fastqs") {
        Some(count) => match count.parse::<u32>() {
            // The flag was present and the value was given and is valid
            Ok(value) => Some(value),
            // The flag was present and the value was given but is not valid.
            Err(_) => None,
        },
        None => None, // The flag was not present
    };

    // Set regex patterns to set minimum length accounting for UMI length
    let num_capture_regex = Regex::new(r"\.\{(\d+)\}").unwrap();
    let base_capture_regex = Regex::new(r"[ATCG]+").unwrap();
    let mut umi_sum: u8 = 0;
        
    // Get length of UMI by combining the numbers in the caputure groups in the regex pattern
    for captures in num_capture_regex.captures_iter(input_regex) {
        if let Ok(num) = captures[1].parse::<u8>() {
            umi_sum += num;
        }
    }

    // Get length of intermediate bases between UMI groups if present
    for capture in base_capture_regex.find_iter(input_regex) {
        umi_sum += capture.as_str().len() as u8;
    }

    // Set the minimum fragment length taking UMI length into account
    let mut minimum_length = *matches.get_one::<u8>("minimum_length").unwrap();
    let umi_length = if umi_sum != 12 { umi_sum } else { 12 };

    minimum_length = 
        if minimum_length != 16 && !is_qiagen { 
            minimum_length + umi_length 
        } else if minimum_length == 16 && !is_qiagen {
            16 + umi_length
        } else {
            16 // Set minimum length for Qiagen to 16 since the UMI is on the 3' end after the adapter
        };

    // Get fastq files
    let fastq_files: Vec<String> = capture_target_files("_R1_001.fastq.gz");
    if fastq_files.is_empty() {
        panic!("No fastq files were found");
    }

    // Check to make sure all input files are gzipped
    for fastq in fastq_files.iter() {
        match is_gzipped(fastq) {
            Ok(compressed) => {
                if !compressed {
                    panic!("File {} is not gzipped. Please gzip the input files.", fastq);
                }
            },
            Err(e) => panic!("Failed to check if file {} is gzipped due to error: {}", fastq, e),
        };
    };

    // Get sample names
    let mut sample_names: Vec<String> = Vec::new();
    for filename in &fastq_files {
        let path = Path::new(&filename);
        let extensionless_name = path.file_stem().unwrap().to_str().unwrap();
        let first_name = extensionless_name.split("_").next().unwrap();
        sample_names.push(first_name.to_string());
    }

    // Create HashMap with the library type and sample names
    let mut sample_library_type: HashMap<String, String> = HashMap::new();
    for sample_name in sample_names.iter() {
        sample_library_type.insert(sample_name.clone(), library_type.to_string());
    }

    // Check quality scores if desired
    let analyze_quality_scores = matches.is_present("quality_score");
    if analyze_quality_scores {
        match average_quality_scores(&fastq_files) {
            Ok(_) => (),
            Err(err) => println!("Error: {:?}", err),
        }
    };

    // Initialize trimmed fastqs vector
    let mut trimmed_fastqs: Vec<String> = vec![];

    // Subsample fastqs if desired
    let subsample_to_lowest = matches.is_present("subsample_to_lowest");

    if subsample_to_lowest {
        // Subsample to the smallest fastq file
        let _ = subsample_fastqs(fastq_files.clone(), sample_names.clone(), None);
        let subsamples_fastq_files = capture_target_files("_subsample.fastq");

        // Trim adapters from the reads and remove UMIs
        trimmed_fastqs = trim_adapters(subsamples_fastq_files, sample_library_type.clone(), minimum_length, &umi_regex);
    } else if let Some(n_reads) = subsample {
        // Subsample to the specified number of reads
        let _ = subsample_fastqs(fastq_files.clone(), sample_names.clone(), Some(n_reads as usize));
        let subsamples_fastq_files = capture_target_files("_subsample.fastq");

        // Trim adapters from the reads and remove UMIs
        trimmed_fastqs = trim_adapters(subsamples_fastq_files, sample_library_type.clone(), minimum_length, &umi_regex);
    } else {
        // No subsampling
        trimmed_fastqs = trim_adapters(fastq_files.clone(), sample_library_type.clone(), minimum_length, &umi_regex);
    }

    // Calculate number of RNA in each sample 
    // Error correct UMIs in sample bam file
    let Ok(sam_files) = rna_discovery_calculation(trimmed_fastqs
                                                                , sample_names.clone()
                                                                , reference, num_threads
                                                                , include_unaligned_reads) 
                                                                else {panic!("An error occurred during the RNA counting calculation")};

    let rna_counts_files = capture_target_files(&format!("_{}_counts", reference_name));

    // Create read length distribution
    let _ = calculate_read_length_distribution(sam_files);

    let _ = threshold_count(rna_counts_files.clone(), thresholds, sample_names.clone(), reference);

    let threshold_files = capture_target_files("_threshold_count_sum.csv");

    let _ = combine_threshold_counts(threshold_files, "RNA_threshold_counts.csv".to_string(), sample_names.clone());

    let (full_rna_names_list, rna_info, rpm_info) = find_common_rnas(rna_counts_files
                                                                     , sample_names.clone()
                                                                     );
    let _ = write_common_rna_file(full_rna_names_list.clone(), rna_info, rpm_info, sample_names, reference_name);

    // Remove intermediate files, such as bam and sam files, unless otherwise specified
    if !keep_intermediates {
        let _ = remove_intermediate_files();
    }

    Ok(())
}