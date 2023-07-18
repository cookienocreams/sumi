use crate::ProgressBar;
use crate::ProgressStyle;
use crate::Command;
use crate::File;
use std::io::Write;

/// Calculates and outputs the distribution of read lengths in input sam files.
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
/// * `sam_files` - A vector of Strings, where each String is the path to a sam file.
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
    let progress_br = ProgressBar::new(num_of_fastqs);
    
    progress_br.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{eta_precise}] [{bar:50.cyan/blue}] {pos}/{len} {msg} ({percent}%)")
                .expect("Progress bar error")
            .progress_chars("#>-"),
    );
    progress_br.set_message("Calculating read length distribution...");

    for sam in sam_files.iter() {
        let sample_name = sam.split(".").next().unwrap();

        // Call samtools to calculate read length distribution
        let output = Command::new("samtools")
            .arg("stats")
            .arg(sam)
            .stdout(std::process::Stdio::piped())
            .spawn()
            .expect("Failed to execute command");

        let output2 = Command::new("grep")
            .arg("^RL")
            .stdin(output.stdout.unwrap())
            .stdout(std::process::Stdio::piped())
            .spawn()
            .expect("Failed to execute command");

        let output3 = Command::new("cut")
            .args(&["-f", "2-"])
            .stdin(output2.stdout.unwrap())
            .output()
            .expect("Failed to execute command");

        let result = String::from_utf8_lossy(&output3.stdout);

        // Write read lenghts to output file
        let file = File::create(format!("{}_read_lengths.csv", sample_name));
        let _ = file.as_ref().expect("Failed to write header").write_all(b"Length\tReads\n");
        let _ = file.as_ref().expect("Failed to write length data").write_all(result.as_bytes());

        progress_br.inc(1);
    }

    progress_br.finish_with_message("Finished counting read lengths");

}