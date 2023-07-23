use crate::{ProgressBar, ProgressStyle};
use crate::Path;
use crate::Command;
use crate::{HashSet, HashMap};
use crate::graph::find_true_umis;
use crate::graph::deduplicate_bam;
use crate::capture_target_files;
use crate::File;
use crate::BufReader;
use crate::DataFrame;
use crate::Series;
use crate::CsvWriter;
use crate::Error;
use std::io::BufRead;
use polars::prelude::{NamedFrom, SerWriter};
use std::io;

// Define a custom error type that can represent both io::Error and bio::io::fastq::Error
// Can be used to determine if an error occurred when reading the input or the fastq file
#[derive(Error, Debug)]
pub enum MyError {
    #[error("I/O error")]
    Io(#[from] io::Error),
    #[error("FASTQ error")]
    Fastq(#[from] bio::io::fastq::Error),
}

/// Checks a SAM file to determine if it contains alignment data.
///
/// The function reads through the SAM file line by line, looking for any line that does 
/// not start with '@', which would indicate it contains alignment data. If such a line 
/// is found, the function will return `Ok(())`, signaling the file does contain alignments.
///
/// If the end of the file is reached without encountering any alignment data (i.e., all 
/// lines start with '@'), an error is returned with the message "SAM file is empty".
///
/// # Arguments
/// * `sam_file` - A string slice that holds the path to the SAM file.
///
/// # Returns
/// * `Ok(())` if the SAM file contains alignment data.
/// * `Err(io::Error)` if the SAM file does not contain alignment data or if there was an 
/// error reading the file.
///
/// # Errors
/// This function will return an error if there was an issue reading the file or if the 
/// file does not contain alignment data.
///
/// # Example
/// ```no_run
/// let result = is_sam_file_empty("/path/to/your/file.sam");
/// match result {
///     Ok(_) => println!("The SAM file contains alignments."),
///     Err(e) => eprintln!("An error occurred: {}", e),
/// }
/// ```
pub fn is_sam_file_empty(sam_file: &str) -> Result<(), io::Error> {
    let file = File::open(sam_file)?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line?;
        if !line.starts_with('@') {
            return Ok(());
        }
    }

    Err(io::Error::new(io::ErrorKind::InvalidData, "SAM file is empty"))
}

/// Perform alignment of trimmed fastq files to the target RNA bowtie2 reference to calculate
/// the number of unique RNA present in each sample. UMI error correction and deduplication
/// is also done. Only aligned reads are output to SAM file in order to speed up UMI deduplication.
///
/// Returns a Vector of paths to SAM files containing RNA counts for each sample.
///
/// # Arguments
/// * `trimmed_fastqs` - Vector of paths to trimmed fastq files.
/// * `library_type` - HashMap mapping sample names to their library types.
/// * `sample_names` - Vector of sample names.
/// * `reference` - Path to the bowtie2 reference.
/// * `num_threads` - Number of threads to use for alignment and conversion processes.
/// * `include_unaligned_reads` - Boolean flag for whether to include unaligned reads or not.
///
/// # Example
/// ```
/// let trimmed_fastqs = vec!["sample1.cut.fastq", "sample2.cut.fastq", "sample3.cut.fastq"];
/// let sample_names = vec!["sample1", "sample2", "sample3"];
///
/// rna_discovery_calculation(trimmed_fastqs, sample_names, "path/to/reference", 12, false);
/// // This will generate "sample1.miRNA.sam", "sample2.miRNA.sam", "sample3.miRNA.sam" files
/// // and will return a vector of SAM file paths.
/// ```
pub fn rna_discovery_calculation(
                                trimmed_fastqs: Vec<String>, 
                                sample_names: Vec<String>,
                                reference: &str,
                                num_threads: u8,
                                include_unaligned_reads: bool
                            ) -> Result<Vec<String>, io::Error> {
    let num_of_fastqs = trimmed_fastqs.len() as u64;
    let progress_br = ProgressBar::new(num_of_fastqs);
    
    progress_br.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:50.cyan/blue}] {pos}/{len} {msg} ({percent}%)")
                .expect("Progress bar error")
            .progress_chars("#>-"),
    );
    let reference_name = Path::new(reference).file_name().unwrap().to_str().unwrap();
    progress_br.set_message(format!("Calculating the number of {} present...", reference_name));

    let unaligned = if include_unaligned_reads {
        ""
    } else {
        "--no-unal"
    };

    for (fastq_file, sample_name) in trimmed_fastqs.iter().zip(sample_names.iter()) {
        // Align to bowtie2 reference
        // Make mapped SAM file with header for faster analysis
        let _ = Command::new("bowtie2")
            .args(&["--norc"
                    , "--threads"
                    , &num_threads.to_string()
                    , unaligned
                    , "-x"
                    , reference
                    , "-U"
                    , fastq_file
                    , "-S"
                    , &format!("{}.UMI.sam", sample_name)]
                )
            .output()
            .expect("Failed to run bowtie2");

        let check_result = is_sam_file_empty(&format!("{}.UMI.sam", sample_name));

        // Skip analysis of sample if no alignments are found
        match check_result {
            Ok(_) => {
                // Continue processing this sample
            }
            Err(e) => {
                println!("Error while checking {}: {}. Please check input samples and parameters."
                        , &format!("{}.UMI.sam", sample_name)
                        , e);
                // Skip this sample
                continue;
            }
        }

        // Can use SAM file for deduplication, but a BAM file uses less memory
        let _ = Command::new("samtools")
            .args(&["view"
                    , &format!("-@ {}", &num_threads).to_string()
                    , "--with-header"
                    , "-o"
                    , &format!("{}.UMI.bam", sample_name)
                    , &format!("{}.UMI.sam", sample_name)]
                )
            .output()
            .expect("Failed to run samtools to convert SAM to BAM");

        // Error correct UMIs in sample bam file
        let repr_umis: HashSet<Vec<u8>> = find_true_umis(&format!("{}.UMI.bam", sample_name))?;

        // Deduplicate mapped bam file
        let deduplication_result = deduplicate_bam(&format!("{}.UMI.bam", sample_name)
                                                                                , sample_name
                                                                                , repr_umis);

        // Check the result
        match deduplication_result {
            Ok(_) => (),
            Err(e) => eprintln!("An error occurred while deduplicating the sample {}: {}", sample_name, e),
        }

        // Convert deduplicated bam file to a sam file
        let _ = Command::new("samtools")
            .args(&["view"
                    , &format!("-@ {}", &num_threads).to_string()
                    , "--with-header"
                    , "-o"
                    , &format!("{}.{}.sam", sample_name, reference_name)
                    , &format!("{}.dedup.bam", sample_name)]
                )
            .output()
            .expect("Failed to run samtools to convert BAM to SAM");

        // Generate data for the specific RNAs captured
        let _  = generate_rna_counts(&format!("{}.{}.sam", sample_name, reference_name), sample_name, reference_name);

        progress_br.inc(1);
    }

    let sam_files: Vec<String> = capture_target_files(&format!(".{}.sam", reference_name));

    // Stop analysis if all samples lack alignment data 
    if sam_files.len() == 0 {
        panic!("No deduplicated sam files found. Analysis aborted.");
    }

    progress_br.finish_with_message(format!("Finished counting {}", reference_name));

    Ok(sam_files)
}

/// This function counts the occurrence of each RNA in an input SAM file, computes their 
/// reads per million (RPM),and writes this information to an output CSV file. The CSV 
/// file will contain the name of each unique RNA, its count,and RPM for the given sample.
///
/// The output CSV file will have the format `"{sample_name}_{reference_name}_counts.csv"` 
/// and will be located in the directory where the program is run.
///
/// # Returns
///
/// * A `Result` which is `Ok` if the function executed successfully. It is `Err` if an error occurred, 
///   with a description of the error.
///
/// # Arguments
///
/// * `input_sam_file`: The path to the input SAM file.
/// * `sample_name`: The name of the sample.
/// * `reference_name`: The name of the bowtie2 reference.
pub fn generate_rna_counts(input_sam_file: &str, sample_name: &String, reference_name: &str) 
    -> Result<(), Box<dyn std::error::Error>> {
    let file = File::open(input_sam_file)?;
    let reader = BufReader::new(file);

    // Initialize a HashMap to count the occurrence of each RNA
    let mut rna_counts: HashMap<String, i32> = HashMap::new();

    // Read the SAM file line by line
    for line in reader.lines() {
        let line = line?;

        // Skip the header lines
        if line.starts_with('@') {
            continue;
        }

        // Split each line into fields and extract the sequence name
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }
        let ref_seq_name = fields[2];

        // Increment the count for this RNA in the HashMap
        *rna_counts.entry(ref_seq_name.to_string()).or_insert(0) += 1;
    }

    // Calculate the total number of mapped reads
    let total_mapped_reads: i32 = rna_counts.values().sum();

    // Prepare vectors for DataFrame creation
    let mut names: Vec<String> = Vec::new();
    let mut counts: Vec<i32> = Vec::new();
    let mut rpms: Vec<f64> = Vec::new();

    // For each RNA, calculate its RPM and store its name, count, and RPM
    for (name, count) in rna_counts {
        names.push(name);
        counts.push(count);
        rpms.push((count as f64 / total_mapped_reads as f64) * 1_000_000.0); // Calculate RPM
    }

    // Create a DataFrame with 'name', 'count', and 'RPM' as columns
    let mut counts_df = DataFrame::new(vec![
        Series::new("name", names),
        Series::new("count", counts),
        Series::new("RPM", rpms),
    ])?;

    // Sort the DataFrame in-place by the 'count' column
    let _ = counts_df.sort_in_place(&["count"], true);

    // Create a new CSV file to store the DataFrame
    let file = File::create(format!("{}_{}_counts.csv", sample_name, reference_name)).expect("could not create file");

    // Write the DataFrame to the CSV file with headers
    CsvWriter::new(file)
        .has_header(true)
        .finish(&mut counts_df)
        .unwrap();

    Ok(())
}