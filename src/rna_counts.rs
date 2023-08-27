use crate::capture_target_files;
use crate::csv_writer;
use crate::get_read_counts;
use crate::graph::deduplicate_bam;
use crate::graph::find_true_umis;
use crate::quality::average_read_quality;
use crate::BufReader;
use crate::Command;
use crate::Config;
use crate::CsvWriter;
use crate::DataFrame;
use crate::Error;
use crate::File;
use crate::Path;
use crate::Series;
use crate::{HashMap, HashSet};
use crate::{ProgressBar, ProgressStyle};
use polars::prelude::{NamedFrom, SerWriter};
use std::io;
use std::io::BufRead;

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

    Err(io::Error::new(
        io::ErrorKind::InvalidData,
        "SAM file is empty",
    ))
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
/// * `config` - Struct containing user-specified analysis configurations, see below.
///     * `reference` - Path to the bowtie2 reference.
///     * `num_threads` - Number of threads to use for alignment and conversion processes.
///     * `include_unaligned_reads` - Boolean flag for whether to include unaligned reads or not.
///     * `write_metrics` - Boolean flag for whether to write sample metrics to an output file or not.
///
/// # Example
/// ```
/// let trimmed_fastqs = vec!["sample1.cut.fastq", "sample2.cut.fastq", "sample3.cut.fastq"];
/// let sample_names = vec!["sample1", "sample2", "sample3"];
/// let config = Config::new(matches);
///
/// rna_discovery_calculation(trimmed_fastqs, sample_names, config);
/// // This will generate "sample1.miRNA.sam", "sample2.miRNA.sam", "sample3.miRNA.sam" files
/// // and will return a vector of SAM file paths.
/// ```
pub fn rna_discovery_calculation(
    trimmed_fastqs: Vec<String>,
    sample_names: Vec<String>,
    config: &Config,
) -> Result<Vec<String>, io::Error> {
    // Calculate the number of reads in each file if metrics are enabled
    let read_counts = if config.write_metrics {
        let fastq_files: Vec<String> = capture_target_files("_R1_001.fastq.gz");
        match get_read_counts(fastq_files, &sample_names) {
            Ok(read_counts) => read_counts,
            Err(err) => panic!("An error occurred while getting read counts: {}", err),
        }
    } else {
        HashMap::new() // Return an empty HashMap when write_metrics is not true
    };

    let num_of_fastqs = trimmed_fastqs.len() as u64;
    let progress_br = ProgressBar::new(num_of_fastqs);

    progress_br.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:50.cyan/blue}] {pos}/{len} {msg} ({percent}%)")
                .expect("Progress bar error")
            .progress_chars("#>-"),
    );

    // Set analysis reference, i.e., RNA type name, based on the user provided reference name
    let reference_name = Path::new(&config.alignment_reference)
        .file_name()
        .unwrap()
        .to_str()
        .unwrap();
    progress_br.set_message(format!(
        "Calculating the number of {} present...",
        reference_name
    ));

    // Set bowtie2 flag to keep unaligned reads
    let unaligned = if config.include_unaligned_reads {
        ""
    } else {
        "--no-unal"
    };

    // Create file for writing metrics if necessary
    let mut writer: Option<csv_writer<File>> = None;
    if config.write_metrics {
        writer = Some(csv_writer::from_path("metrics.csv")?);
        writer.as_mut().unwrap().write_record([
            "sample",
            "read count",
            "quality score",
            reference_name,
        ])?;
    }

    for (fastq_file, sample_name) in trimmed_fastqs.iter().zip(sample_names.iter()) {
        // Align to bowtie2 reference
        // Make mapped SAM file with header for faster analysis
        let _ = Command::new("bowtie2")
            .args([
                "--norc",
                "--threads",
                &config.num_threads.to_string(),
                unaligned,
                "-x",
                &config.alignment_reference,
                "-U",
                fastq_file,
                "-S",
                &format!("{}.UMI.sam", sample_name),
            ])
            .output()
            .expect("Failed to run bowtie2");

        let check_result = is_sam_file_empty(&format!("{}.UMI.sam", sample_name));

        // Skip analysis of sample if no alignments are found
        match check_result {
            Ok(_) => {
                // Continue processing this sample
            }
            Err(e) => {
                println!(
                    "Error while checking {}: {}. Please check input samples and parameters.",
                    &format!("{}.UMI.sam", sample_name),
                    e
                );
                // Skip this sample
                continue;
            }
        }

        // Can use SAM file for deduplication, but a BAM file uses less memory
        let _ = Command::new("samtools")
            .args([
                "view",
                &format!("-@ {}", &config.num_threads).to_string(),
                "--with-header",
                "-o",
                &format!("{}.UMI.bam", sample_name),
                &format!("{}.UMI.sam", sample_name),
            ])
            .output()
            .expect("Failed to run samtools to convert SAM to BAM");

        // Error correct UMIs in sample bam file
        let repr_umis: HashSet<Vec<u8>> = find_true_umis(
            &format!("{}.UMI.bam", sample_name),
            config.levenshtein_distance,
        )?;

        // Deduplicate mapped bam file
        let deduplication_result =
            deduplicate_bam(&format!("{}.UMI.bam", sample_name), sample_name, repr_umis);

        // Check the result
        match deduplication_result {
            Ok(_) => (),
            Err(e) => eprintln!(
                "An error occurred while deduplicating the sample {}: {}",
                sample_name, e
            ),
        }

        // Convert deduplicated bam file to a sam file
        let _ = Command::new("samtools")
            .args([
                "view",
                &format!("-@ {}", &config.num_threads).to_string(),
                "--with-header",
                "-o",
                &format!("{}.{}.sam", sample_name, reference_name),
                &format!("{}.dedup.bam", sample_name),
            ])
            .output()
            .expect("Failed to run samtools to convert BAM to SAM");

        // Generate data for the specific RNAs captured
        let _ = generate_rna_counts(
            &format!("{}.{}.sam", sample_name, reference_name),
            sample_name,
            reference_name,
        );

        // Set correct read count depending on if subsampling was performed
        if config.write_metrics {
            let count: u64;
            if let Some(n_reads) = config.subsample_fastqs {
                count = n_reads as u64 // User chose the number of reads to subsample to
            } else if config.subsample_to_lowest {
                count = *read_counts.values().min().unwrap() // Count equals fastq with lowest read count
            } else {
                count = read_counts[sample_name] // Total read count
            };

            // Calculate the percentage of aligned reads
            if let Ok(percent_alignment) =
                get_percent_alignment(sample_name, config.num_threads, reference_name, count)
            {
                let avg_quality =
                    average_read_quality(fastq_file).expect("Failed to calculate read quality");

                writer
                    .as_mut()
                    .expect("Failed to write record")
                    .write_record([
                        sample_name,
                        &count.to_string(),
                        &avg_quality.to_string(),
                        &percent_alignment.to_string(),
                    ])?;
            };
        }

        progress_br.inc(1);
    }

    if config.write_metrics {
        writer.expect("Failed to flush writer").flush()?;
    }

    let sam_files: Vec<String> = capture_target_files(&format!(".{}.sam", reference_name));

    // Stop analysis if all samples lack alignment data
    if sam_files.is_empty() {
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
/// # Example
/// ```
/// let sam_file = "sample1.sam";
/// let reference_name = "miRNA"
///
/// let _ = generate_rna_counts(sam_file, "sample1", reference_name);
/// // This will generate a file "sample1_miRNA_counts.csv" with each miRNA found and
/// // their counts and CPM.
/// ```
pub fn generate_rna_counts(
    input_sam_file: &str,
    sample_name: &String,
    reference_name: &str,
) -> Result<(), Box<dyn std::error::Error>> {
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
    let _ = counts_df.sort_in_place(["count"], true);

    // Create a new CSV file to store the DataFrame
    let file = File::create(format!("{}_{}_counts.csv", sample_name, reference_name))
        .expect("could not create file");

    // Write the DataFrame to the CSV file with headers
    CsvWriter::new(file)
        .has_header(true)
        .finish(&mut counts_df)
        .unwrap();

    Ok(())
}

/// Get the percentage of reads that aligned to the specified reference RNA type.
///
/// Divide reads aligned with total reads to calculate percent alignment.  
///
/// /// # Arguments
///
/// * `count` - The total number of reads analyzed.
/// * `sample_name` - The name of the sample.
/// * `reference_name` - Path to the bowtie2 reference.
/// * `num_threads` - Number of threads to use for alignment and conversion processes.
///
/// # Returns
///
/// The percent each fastq aligned to the specified RNA type.
///
/// # Example
/// ```
/// let reference_name = "miRNA"
/// let count = 12345678
///
/// let percent_alignment = get_percent_alignment("sample1", 4, reference_name, count);
/// println!("The fastq contains {} percent {}.", percent_alignment, reference_name);
/// ```
pub fn get_percent_alignment(
    sample_name: &String,
    num_threads: u8,
    reference_name: &str,
    count: u64,
) -> Result<f64, Box<dyn std::error::Error>> {
    // Get the number of RNA alignments in the input fastq
    let output_result = Command::new("samtools")
        .args([
            "view",
            &format!("-@ {}", &num_threads),
            "-c",
            "-F 4",
            &format!("{}.{}.sam", sample_name, reference_name),
        ])
        .output();

    let output = output_result.expect("Failed to execute command");
    let stdout = match String::from_utf8(output.stdout) {
        Ok(v) => v,
        Err(e) => panic!("Invalid UTF-8 sequence: {}", e),
    };

    let parsed = match stdout.trim().parse::<f64>() {
        Ok(v) => v,
        Err(e) => {
            eprintln!("Failed to parse output to f64: {:?}", e);
            0.0 // Set read count to zero if there's no output
        }
    };
    let percent_alignment = 100.0 * parsed / count as f64;

    Ok(percent_alignment)
}
