use crate::create_progress_bar;
use crate::Command;
use crate::File;
use crate::Write;
use crate::Error;
use crate::HashMap;
use polars::prelude::*;

/// Calculate and output the distribution of read lengths in input sam files.
///
/// This function iterates over each provided sam file, executes a sequence of command line
/// operations (samtools, grep, cut) to calculate the distribution of read lengths, and
/// writes the results to a new csv file named `<sam_file_name>_read_lengths.csv`.
///
/// The csv file includes two columns: `Length` and `Reads`, where `Length` corresponds to
/// the length of the reads and `Reads` corresponds to the number of reads of that length.
///
/// # Arguments
///
/// * `sam_files`: A vector of Strings, where each String is the path to a sam file.
///
/// # Panics
///
/// The function will panic if it fails to execute the command line operations, or if it
/// fails to create or write to the output csv file.
///
/// # Examples
///
/// ```
/// calculate_read_length_distribution(vec!["file1.sam", "file2.sam", "file3.sam"]);
/// ```
pub fn calculate_read_length_distribution(sam_files: Vec<String>) {
    let num_of_fastqs = sam_files.len() as u64;
    let progress_br = create_progress_bar(num_of_fastqs, "Calculating read length distribution...".to_string());

    for sam in sam_files.iter() {
        let sample_name = sam.split('.').next().expect("Failed to split sample name");

        // Call samtools to calculate read length distribution
        let output = Command::new("samtools")
            .arg("stats")
            .arg(sam)
            .stdout(std::process::Stdio::piped())
            .spawn()
            .expect("Failed to execute samtools command");

        let output2 = Command::new("grep")
            .arg("^RL")
            .stdin(output.stdout.unwrap())
            .stdout(std::process::Stdio::piped())
            .spawn()
            .expect("Failed to execute grep command");

        let output3 = Command::new("cut")
            .args(["-f", "2-"])
            .stdin(output2.stdout.unwrap())
            .output()
            .expect("Failed to execute cut command");

        let result = String::from_utf8_lossy(&output3.stdout);

        // Write read lenghts to output file
        let file = File::create(format!("{}_read_lengths.tsv", sample_name));
        match file
            .as_ref()
            .expect("Failed to write header")
            .write_all(b"Length\tReads\n") {
                Ok(_) => (),
                Err(err) => eprintln!("Error writing length file header: {}", err)
            };
        match file
            .as_ref()
            .expect("Failed to write length data")
            .write_all(result.as_bytes()) {
                Ok(_) => (),
                Err(err) => eprintln!("Error writing length file data: {}", err)
            };

        progress_br.inc(1);
    }

    progress_br.finish_with_message("Finished counting read lengths");
}

/// Calculate and output the distribution of aligned read lengths from deduplicated isomiR sequences.
///
/// This function iterates over each isomiR sequence and uses map to calculate each sequence's 
/// read length. It then writes the results to a new csv file named either 
/// `<sample_name>_read_lengths.csv` or `<sample_name>_isomiR_read_lengths.csv`.
///
/// The csv file includes two columns: `Length` and `Reads`, where `Length` corresponds to
/// the length of the reads and `Reads` corresponds to the number of reads of that length.
///
/// # Arguments
///
/// * `dedup_seqs`: A vector of Strings, where each String is an RNA sequence.
/// * `sample_name`: The name of the sample being processed.
/// * `are_isomirs`: Boolean flag to indicate whether isomiRs or miRNA are being analyzed.
///
/// # Panics
///
/// The function will panic if it fails to execute the command line operations, or if it
/// fails to create or write to the output csv file.
///
/// # Examples
///
/// ```
/// calculate_isomer_read_length_distribution(vec!["CGTTGC", "GCATGCA"], "sample1", false);
/// ```
pub fn calculate_isomer_read_length_distribution(
    dedup_seqs: Vec<String>, 
    sample_name: &str,
    are_isomirs: bool
) -> Result<(), Box<dyn Error>> {
    let seq_lens: Vec<usize> = dedup_seqs.iter().map(|seq| seq.len()).collect();

    let mut sequence_lengths: HashMap<i32, i32> = HashMap::new();
    for seq in seq_lens {
        *sequence_lengths.entry(seq.try_into()?).or_insert(0) += 1;
    }

    let lengths: Vec<i32> = sequence_lengths.clone().into_keys().collect::<Vec<i32>>();
    let length_counts = sequence_lengths.into_values().collect::<Vec<i32>>();

    let mut lengths_df = 
        DataFrame::new(vec![
            Series::new("Length", lengths),
            Series::new("Reads", length_counts),
        ])?;
    lengths_df.sort_in_place(["Length"], false, true)?;

    let length_file = 
        if are_isomirs {
            File::create(format!("{}_isomiR_read_lengths.csv", sample_name))?
        } else {
            File::create(format!("{}_read_lengths.csv", sample_name))?
        };

    // Write the length DataFrame to the CSV file with headers
    CsvWriter::new(length_file)
        .has_header(true)
        .finish(&mut lengths_df)?;
    
    Ok(())
}