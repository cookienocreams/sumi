use crate::get_reader;
use crate::create_progress_bar;
use crate::BufReader;
use crate::BufWriter;
use crate::File;
use crate::HashMap;
use crate::Path;
use rand::{Rng, SeedableRng};
use std::io::{BufRead, Write};
use bio::io::fastq;
use bio::io::fastq::Record;

/// Count the number of reads in a FASTQ file.
///
/// # Arguments
///
/// * `input_fastq` - A string slice that holds the name of the FASTQ file to read.
///
/// # Returns
///
/// This function will return the number of reads in the input FASTQ file as a `Result<u64,
/// Box<dyn std::error::Error>>`. If an error occurs while reading the file, the error
/// will be returned instead.
///
/// # Examples
///
/// ```
/// let num_reads = count_reads("my_file.fastq")?;
/// println!("The file contains {} reads.", num_reads);
/// ```
pub fn count_reads(input_fastq: &str) -> Result<u64, Box<dyn std::error::Error>> {
    let reader = File::open(input_fastq)?;
    let (mut reader, _compression) = get_reader(Box::new(reader))?;

    let line_count = BufReader::new(&mut reader)
        .lines()
        .map_while(Result::ok)
        .count();

    Ok((line_count / 4) as u64)
}

/// Find the FASTQ file with the fewest reads from a list of files.
///
/// # Arguments
///
/// * `files` - A vector of string slices where each slice holds the name of a FASTQ file.
///
/// # Returns
///
/// This function will return a tuple as a `Result<(String, u64), Box<dyn std::error::Error>>`
/// with the name of the file with the fewest reads and the count of those reads.
/// If an error occurs while reading any of the files, the error will be returned instead.
///
/// # Examples
///
/// ```
/// let files = vec!["file1.fastq", "file2.fastq", "file3.fastq"];
/// let (file, count) = find_fastq_with_fewest_reads(files)?;
/// println!("The file with the fewest reads is {} with {} reads.", file, count);
/// ```
pub fn find_fastq_with_fewest_reads(
    input_fastqs: &[String],
    read_counts: &mut HashMap<String, usize>,
) -> Result<(String, u64), Box<dyn std::error::Error>> {
    let num_of_fastqs = input_fastqs.len() as u64;
    let progress_br = create_progress_bar(num_of_fastqs, "Counting the reads in each fastq files...".to_string());

    let mut min_reads = u64::MAX;
    let mut min_index = 0;

    for (index, file) in input_fastqs.iter().enumerate() {
        // Check if the count is already stored in the Hashmap
        let count = match read_counts.get(file) {
            Some(count) => *count,
            None => {
                // If not, calculate it and store it in the Hashmap
                let count = count_reads(file)?.try_into().unwrap();
                read_counts.insert(file.to_string(), count);
                count
            }
        };

        if count < min_reads as usize {
            min_reads = count as u64;
            min_index = index;
        }

        progress_br.inc(1);
    }

    progress_br.finish_with_message("Finished counting the number of reads in the fastq files");

    Ok((input_fastqs[min_index].clone(), min_reads))
}

/// Subsamples a given FASTQ file to a target number of reads using Algorithm R for reservoir sampling.
///
/// # Arguments
///
/// * `input_file` - A string reference that holds the name of the input FASTQ file.
/// * `target_read_count` - A u32 that holds the desired number of reads in the output.
/// * `rng_seed` - An Option containing either a specified RNG for subsampling or `None`
/// which then uses the thread RNG seeded by the system.
///
/// # Description
///
/// This function randomly selects a subset of reads from an input FASTQ file using
/// Algorithm R for reservoir sampling. The number of reads selected is determined by
/// the `target_read_count` argument. The input file is read only once, and the selection
/// process does not require knowledge of the total number of reads in the input file.
///
/// # Returns
///
/// This function returns a Result containing a vector of Records. Each Record
/// represents a selected read from the input file and its corresponding lines
/// in the FASTQ file.
///
/// # Errors
///
/// This function will return an error if there is a problem reading the input file.
pub fn subsample_file(
    input_file: &str,
    target_read_count: u32,
    rng_seed: Option<u64>,
) -> Result<Vec<Record>, Box<dyn std::error::Error>> {
    // Create a random number generator
    let mut rng = match rng_seed {
        Some(rng) => rand_chacha::ChaCha8Rng::seed_from_u64(rng),
        None => rand_chacha::ChaCha8Rng::from_rng(rand::thread_rng())?,
    };

    // Open the file and get a reader that can handle potential compression
    let reader = File::open(input_file)?;
    let (mut reader, _compression) = get_reader(Box::new(reader))?;
    let reader = BufReader::new(&mut reader);

    // Create a FASTQ reader from the BufReader
    let mut fastq_reader = fastq::Reader::from_bufread(reader).records();

    // Initialize reservoir and total reads count
    let mut reservoir: Vec<Record> = Vec::with_capacity(target_read_count.try_into()?);
    let mut reads_counted = 0;

    // While there are still lines in the file...
    while let Some(Ok(record)) = fastq_reader.next() {
        reads_counted += 1; // Increment the total read count

        if reads_counted <= target_read_count {
            // If we have not yet reached the target number of reads, just store the read in the reservoir
            reservoir.push(record);
        } else {
            // Generate a random number between 0 and the total read count (inclusive)
            let random_num: usize = rng.gen_range(0..=reads_counted).try_into()?;

            // If the random number is less than the target read count, replace a read in the reservoir
            if random_num < target_read_count.try_into()? {
                reservoir[random_num] = record;
            }
        }
    }

    Ok(reservoir)
}

/// Writes selected reads to an output FASTQ file.
///
/// # Arguments
///
/// * `output_file` - A string reference that holds the name of the output FASTQ file.
/// * `reservoir` - A vector of Records. Each Record represents a selected read and its 
/// corresponding lines in the FASTQ file.
///
/// # Description
///
/// This function writes the reads selected by the `subsample_file` function to an output FASTQ file.
///
/// # Returns
///
/// This function returns a Result type indicating whether the writing was successful. If the writing is successful,
/// the function returns Ok(()). If there is a problem writing to the file, the function returns an error.
///
/// # Errors
///
/// This function will return an error if there is a problem writing to the output file.
pub fn write_to_output(
    output_file: &str,
    reservoir: Vec<Record>,
) -> Result<(), Box<dyn std::error::Error>> {
    // Write the selected reads to the output FASTQ file
    let file = File::create(output_file)?;
    let mut writer = BufWriter::new(file);

    for read in &reservoir {
        write!(writer, "{}", read)?;
    }

    Ok(())
}

/// Subsample multiple FASTQ files each to a target number of reads using Algorithm R for reservoir sampling.
///
/// # Arguments
///
/// * `input_fastqs` - A vector of Strings that hold the names of the input FASTQ files.
/// * `sample_names` - A vector of Strings that hold the names of the samples.
/// * `target_read_count` - A u32 that holds the desired number of reads in the output files.
/// * `rng_seed` - An Option containing either a specified RNG for subsampling or `None`
/// which then uses the thread RNG seeded by the system.
///
/// # Description
///
/// This function randomly selects a subset of reads from each input FASTQ file and writes them
/// to corresponding output FASTQ files, with the names based on `sample_names` followed by a `_subsample.fastq` suffix.
/// The number of reads selected is determined by the `target_read_count` argument.
/// If `target_read_count` is None, the function selects a number equal to the smallest number
/// of reads across all input files. This function uses Algorithm R for reservoir sampling to select a random
/// subset of reads from the input file. The input file is read only once, and the selection process does
/// not require knowledge of the total number of reads in the input file. Therefore, this function is
/// more memory-efficient and can handle larger files than methods that require reading the entire file
/// into memory.
///
/// Algorithm R works by filling the reservoir with the first `k` items from the input. Each remaining item is then given
/// a chance to be included in the reservoir with a probability of `k/i`, where `i` is the current item index in the input.
/// If the item is selected, it replaces a random item already in the reservoir. This process ensures each item from the
/// input has an equal chance to end up in the reservoir.
///
/// # Returns
///
/// This function does not return a value. It writes the selected reads to the output FASTQ files,
/// each named based on corresponding `sample_name` with a `_subsample.fastq` suffix.
///
/// # Example
///
/// ```
/// let input_fastqs = vec!["sample1.fastq.gz", "sample2.fastq.gz"];
/// let sample_names = vec!["sample1", "sample2"];
/// subsample_fastqs(input_fastqs, sample_names, Some(1000)).expect("Failed to subsample fastqs");
/// ```
pub fn subsample_fastqs(
    input_fastqs: Vec<String>,
    sample_names: Vec<String>,
    target_read_count: Option<u32>,
    rng_seed: Option<u64>,
) -> Result<(), Box<dyn std::error::Error>> {
    // Add a HashMap to store the read counts
    let mut read_counts: HashMap<String, usize> = HashMap::new();

    // Determine the target read count
    let (smallest_fastq_file, target_read_count) = match target_read_count {
        Some(count) => (String::new(), count), // If a target_read_count is provided, set smallest_fastq_file to an empty string
        None => {
            let (file, count) = find_fastq_with_fewest_reads(&input_fastqs, &mut read_counts)?;
            (file, count.try_into()?) // If not, find the file with the smallest read count
        }
    };

    // If smallest_fastq_file is not an empty string, it means it's the file with the minimum read count. So, copy it.
    if !smallest_fastq_file.is_empty() {
        // Copy the smallest fastq file with a new name to match the processed files
        let smallest_fastq_file_path = Path::new(&smallest_fastq_file);
        let smallest_fastq_file_stem = smallest_fastq_file_path
            .file_stem()
            .unwrap()
            .to_str()
            .unwrap();
        let smallest_fastq_file_extension = smallest_fastq_file_path
            .extension()
            .unwrap()
            .to_str()
            .unwrap();

        let smallest_fastq_file_copy = format!(
            "{}_subsample.{}",
            smallest_fastq_file_stem, smallest_fastq_file_extension
        );
        std::fs::copy(&smallest_fastq_file, smallest_fastq_file_copy)?;
    }

    // Create a progress bar
    let num_of_fastqs = input_fastqs.len() as u64;
    let progress_br = create_progress_bar(num_of_fastqs, "Subsampling fastq files...".to_string());

    // Subsample fastq files
    for (fastq, sample_name) in input_fastqs.iter().zip(sample_names.iter()) {
        let reservoir = subsample_file(fastq, target_read_count, rng_seed)?;
        let output_filename = format!("{}_subsample.fastq", sample_name);
        write_to_output(&output_filename, reservoir)?;

        progress_br.inc(1);
    }

    progress_br.finish_with_message("Finished subsampling fastq files");

    Ok(())
}
