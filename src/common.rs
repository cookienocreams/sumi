use crate::BufReader;
use crate::File;
use crate::HashMap;
use std::io::BufRead;
use std::io::Write;

/// Function to find RNAs that are common in all given samples, and their counts, and RPM.
///
/// # Arguments
///
/// * `rna_counts_files` - A vector of files containing sample RNA counts.
/// * `read_count_dict` - A dictionary where keys are sample names and values are the total number of reads
/// for the corresponding sample.
/// * `sample_names` - A vector of sample names.
///
/// # Returns
///
/// * A tuple containing three fields:
///     1. A list of RNAs that are common to all samples.
///     2. A list of dictionaries where each dictionary contains the RNA counts for one sample.
///     3. A list of dictionaries where each dictionary contains the RPM values for the RNAs in one sample.
///
/// # Example
///
/// ```
/// let common_rnas = find_common_rnas(
///     vec!["sample1_miRNA_counts.csv", "sample2_miRNA_counts.csv"],
///     [("sample1".to_string(), 6390314), ("sample2".to_string(), 5000000)].into_iter().collect(),
///     vec!["sample1", "sample2"]
/// );
/// ```
pub fn find_common_rnas(
    rna_counts_files: Vec<String>,
    sample_names: Vec<String>,
) -> (
    Vec<std::string::String>,
    Vec<HashMap<std::string::String, i32>>,
    Vec<HashMap<std::string::String, f64>>,
) {
    let mut rna_info = vec![];
    let mut rpm_info = vec![];
    let mut rna_names_dict: HashMap<String, usize> = HashMap::new(); // Changed i32 to usize

    for (file, _) in rna_counts_files
        .into_iter()
        .zip(sample_names.clone().into_iter())
    {
        let file = File::open(file).unwrap(); // Use File::open instead of Reader::from_path
        let reader = BufReader::new(file); // Use BufReader to read the file
        let mut sample_rna_counts_dictionary = HashMap::new();
        let mut sample_rna_rpm_dictionary = HashMap::new();

        let mut lines = reader.lines();
        lines.next(); // Skip the first line which is the header

        for line in lines {
            let record = line.unwrap();
            let fields: Vec<&str> = record.split(',').collect(); // Use split to split the record into fields
            let rna_name = fields[0].to_string();
            let rna_count: usize = fields[1].parse().unwrap(); // Parse count as usize
            let rpm: f64 = fields[2].parse().unwrap();

            // Store RNA count and RPM in the dictionaries for the current sample
            sample_rna_counts_dictionary.insert(rna_name.clone(), rna_count as i32); // Convert rna_count to i32
            sample_rna_rpm_dictionary.insert(rna_name.clone(), rpm);

            // Increment the count for the current RNA name in the global dictionary
            *rna_names_dict.entry(rna_name).or_insert(0) += 1;
        }

        // Add the dictionaries for the current sample to the storage vectors
        rna_info.push(sample_rna_counts_dictionary);
        rpm_info.push(sample_rna_rpm_dictionary);
    }

    // Find RNA present in all samples, i.e., have counts equal to the number of samples analyzed
    let common_rnas: Vec<String> = rna_names_dict
        .into_iter()
        .filter(|(_, count)| *count == sample_names.len())
        .map(|(name, _)| name)
        .collect();

    // Return the common rnas, and the info about counts and RPM
    (common_rnas, rna_info, rpm_info)
}

/// Function to write RNAs common to all samples to two output files: one for counts and one for RPM values.
/// Each row in the output files represents a RNA and each subsequent column represents the counts
/// (or RPM values) for that RNA in one sample.
///
/// # Arguments
///
/// * `rna_names` - A vector of RNA names that are common to all samples.
/// * `rna_info` - A vector of HashMaps where each HashMap contains the RNA counts for one sample.
/// * `rpm_info` - A vector of HashMaps where each HashMap contains the RPM values for the RNAs in one sample.
/// * `sample_names` - A vector of sample names.
/// * `reference_name` - Name of the reference.
///
/// # Returns
///
/// Nothing, the function writes to output files.
///
/// # Example
///
/// ```
/// let rna_names = vec!["mirna1", "mirna2", "mirna3"];
/// let rna_info = vec![
///     hashmap!{ "mirna1" => 10, "mirna2" => 20, "mirna3" => 30 },
///     hashmap!{ "mirna1" => 40, "mirna2" => 50, "mirna3" => 60 },
/// ];
/// let rpm_info = vec![
///     hashmap!{ "mirna1" => 0.1, "mirna2" => 0.2, "mirna3" => 0.3 },
///     hashmap!{ "mirna1" => 0.4, "mirna2" => 0.5, "mirna3" => 0.6 },
/// ];
/// let sample_names = vec!["sample1", "sample2"];
///
/// write_common_rna_file(
///     rna_names,
///     rna_info,
///     rpm_info,
///     sample_names,
///     "miRNA",
/// );
/// ```
pub fn write_common_rna_file(
    rna_names: Vec<String>,
    rna_info: Vec<HashMap<String, i32>>,
    rpm_info: Vec<HashMap<String, f64>>,
    sample_names: Vec<String>,
    reference_name: &str,
) {
    let file_path = &format!("Common_{}s.tsv", reference_name);
    let file_path_rpm = &format!("Common_{}s_RPM.tsv", reference_name);

    // Open the output files
    let mut common_rna_file = File::create(file_path).unwrap();
    let mut common_rna_file_rpm = File::create(file_path_rpm).unwrap();

    // Write the headers to the output files
    let header: String = sample_names.join("\t");
    writeln!(common_rna_file, "{}\t{}", reference_name, header).unwrap();
    writeln!(common_rna_file_rpm, "{}\t{}", reference_name, header).unwrap();

    // For each RNA, write its counts and RPM values in all samples to the output files
    for rna in &rna_names {
        // Start the output lines with the RNA name
        let mut counts_line = format!("{}\t", rna);
        let mut rpms_line = format!("{}\t", rna);

        // Add the count or RPM value for the current RNA in each sample to the output lines
        for (index, _sample) in sample_names.iter().enumerate() {
            counts_line.push_str(&format!("{}\t", rna_info[index][rna]));
            rpms_line.push_str(&format!("{}\t", rpm_info[index][rna]));
        }

        // Write the output lines to the output files, removing the trailing tab character if present
        writeln!(common_rna_file, "{}", counts_line.trim_end_matches('\t')).unwrap();
        writeln!(common_rna_file_rpm, "{}", rpms_line.trim_end_matches('\t')).unwrap();
    }
}
