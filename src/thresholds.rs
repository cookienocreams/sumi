use crate::ProgressBar;
use crate::ProgressStyle;
use crate::Path;
use crate::fs::OpenOptions;
use crate::Error;
use crate::HashMap;
use crate::csv_reader;
use crate::csv_writer;
use crate::BufReader;
use crate::File;
use std::io::Write;

/// Process a list of CSV files containing RNA counts and create output CSV files
/// which summarize the number of RNAs whose counts exceed a set of thresholds.
///
/// The input CSV files should have the following structure:
/// - The first column contains the RNA names.
/// - The second column contains the raw count of the RNA.
/// - The third column contains the RPM count.
///
/// The output CSV files have the following structure:
/// - The first row contains the name of the input file.
/// - Each subsequent row contains a threshold and the number of RNA from the
///   corresponding input file whose count exceeds that threshold.
///
/// # Arguments
///
/// * `file_names`: A vector of strings where each string is the path to an input CSV file.
/// * `thresholds`: A vector of usize values where each value is a threshold.
///
/// # Errors
///
/// This function will return an error if:
/// - There is a problem reading an input CSV file.
/// - There is a problem writing to an output CSV file.
/// - A record in an input CSV file cannot be deserialized into a tuple of a string and two numbers.
///
/// # Examples
///
/// ```
/// let file_names = vec!["sample1_miRNA_counts.csv", "sample1_miRNA_counts.csv"];
/// let thresholds = vec![1, 3, 5, 10];
/// threshold_count(file_names, thresholds).unwrap();
/// ```
pub fn threshold_count(file_names: Vec<String>, thresholds: Vec<usize>, sample_names: Vec<String>, reference: &str) 
    -> Result<(), Box<dyn Error>> {
    let num_of_files = file_names.len() as u64;
    let progress_br = ProgressBar::new(num_of_files);
        
    progress_br.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{eta_precise}] [{bar:50.cyan/blue}] {pos}/{len} {msg} ({percent}%)")
                .expect("Progress bar error")
            .progress_chars("#>-"),
    );
    let reference_name = Path::new(reference).file_name().unwrap().to_str().unwrap();
    progress_br.set_message(format!("Calculating {} threshold counts...", reference_name));

    // Loop over all files
    for (file_name, sample_name) in file_names.iter().zip(sample_names.iter()) {

        let mut output = OpenOptions::new()
            .write(true)
            .create(true)
            .open(format!("{}_threshold_count_sum.csv", sample_name))?;

        writeln!(output, ",{}", sample_name)?;

        // Loop over all thresholds
        for threshold in &thresholds {
            let mut rdr = csv::ReaderBuilder::new()
                .delimiter(b',')
                .from_path(&file_name)?;

            let mut above_threshold_count = 0;

            // Loop over all records in the current file
            for result in rdr.deserialize() {
                let record: (String, usize, f64) = result?;
                if record.1 >= *threshold {
                    above_threshold_count += 1;
                }
            }

        writeln!(output, "Threshold {},{}", threshold, above_threshold_count)?;
        }

        progress_br.inc(1);
    }

    progress_br.finish_with_message(format!("Finished counting {} thresholds", reference_name));

    Ok(())
}

/// Combine multiple threshold counts files into a single CSV file.
///
/// The output CSV file will have one row for each unique threshold across all input files,
/// and one column for each sample. Each cell contains the count for that sample at that threshold.
///
/// # Arguments
///
/// * `input_files`: A vector of strings where each string is the path to an input threshold counts CSV file.
/// * `output_file`: A string that is the path to the output CSV file.
/// * `sample_names`: A vector of strings where each string is the name of the sample.
///
/// # Errors
///
/// This function will return an error if:
/// - There is a problem reading an input CSV file.
/// - There is a problem writing to the output CSV file.
/// - A record in an input CSV file cannot be deserialized into a tuple of a string and a number.
pub fn combine_threshold_counts(input_files: Vec<String>, output_file: String, sample_names: Vec<String>) 
    -> Result<(), Box<dyn Error>> {
    // Create a HashMap to hold the count data. The outer key is the threshold, and the value is another HashMap
    // where the key is the sample name and the value is the count.
    let mut data: HashMap<usize, HashMap<&String, usize>> = HashMap::new();

    for (input_file, sample_name) in input_files.iter().zip(sample_names.iter()) {
        // Create a CSV reader for the input file
        let reader = csv_reader::from_reader(BufReader::new(File::open(input_file)?));

        for result in reader.into_records() {
            let record = result?;
            let parts: Vec<&str> = record[0].split_whitespace().collect();
            if parts.len() < 2 {
                return Err(From::from(format!("Malformed threshold record: {:?}", record)));
            }
            let threshold: usize = parts[1].parse()?;
            let count: usize = record[1].parse()?;

            // Insert the count data into the HashMap
            data.entry(threshold)
                .or_insert_with(HashMap::new)
                .insert(sample_name, count);
        }
    }

    // Create a CSV writer for the output file
    let mut writer = csv_writer::from_writer(File::create(output_file)?);

    // Write the header row to the output CSV file
    let header: Vec<String> = [" ".to_string()]
        .iter()
        .cloned()
        .chain(sample_names.iter().cloned())
        .collect();
    writer.write_record(&header)?;

    // Get a sorted list of all unique thresholds
    let mut thresholds: Vec<usize> = data.keys().cloned().collect();
    thresholds.sort_unstable();

    // Write each threshold and its count data to the output CSV file
    for threshold in thresholds {
        let mut row: Vec<String> = vec![format!("Threshold {}", threshold)];

        // For each sample, look up the count for the current threshold and add it to the row
        for sample_name in &sample_names {
            let count = data.get(&threshold).and_then(|x| x.get(sample_name)).unwrap_or(&0);
            row.push(count.to_string());
        }
        writer.write_record(&row)?;
    }

    // Flush the writer to ensure all data is written to the file
    writer.flush()?;

    Ok(())
}