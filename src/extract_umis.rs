use crate::is_gzipped;
use crate::Read;
use crate::Regex;
use crate::DUAL_INDEX_REGEX;
use crate::SINGLE_INDEX_REGEX;
use crate::UMI_REGEX_QIAGEN;
use bio::io::fastq::{Reader, Record, Writer};
use flate2::bufread::MultiGzDecoder;
use std::fs::File;
use std::io::BufReader;

/// Extracts Unique Molecular Identifiers (UMIs) from a given fastq file and writes the
/// UMI-trimmed sequences to a new fastq file. This function is flexible in the definition of
/// the UMI, as the user can provide a regular expression with multiple capture groups to specify
/// the location and format of the UMI within the sequence.
///
/// # Arguments
/// * `sample_name` - The name of the sample.
/// * `umi_delineator` - Delimiter used between the appended UMI and the read name. The
/// default is "_" as is used in UMI-tools.
/// * `umi_regex` - A regular expression used to capture the UMI from the read sequence.
/// The expression may contain multiple capture groups to allow for flexible UMI definitions.
///
/// # Example
/// ```rust
/// use regex::Regex;
///
/// let umi_regex = Regex::new(r"(^.{12})").unwrap();
/// extract_umis("sample1", "_", &umi_regex);
/// // This will create an output file named "sample1.cut.fastq" with UMI-trimmed sequences.
/// ```
pub fn extract_umis(
    sample_name: &str,
    umi_delineator: &str,
    umi_regex: &Regex,
) -> Result<(), Box<dyn std::error::Error>> {
    let filename = format!("{}.{}.{}", sample_name, "unprocessed.cut", "fastq");

    // Attempt to open the input file
    let reader = Reader::from_file(filename)?;
    let mut writer = Writer::to_file(format!("{}.cut.fastq", sample_name))?;

    // Loop over all the records in the input file
    for record_result in reader.records() {
        let record = record_result?;

        // Convert the sequence to a string
        let read_sequence = std::str::from_utf8(record.seq())?;

        // Check if the sequence matches the UMI regex
        if let Some(captures) = umi_regex.captures(read_sequence) {
            let umi_str: String;
            let umi_start: usize;
            let umi_len: usize;

            // Check if there is more than one capture group in the regex
            if captures.len() > 2 {
                // Need length of full matched sequence, i.e., the UMI length + other bases
                let umi_match = captures.get(0).unwrap();
                let mut ids = vec![];
                for i in 1..captures.len() {
                    let id = captures[i].to_owned();
                    ids.push(id);
                }
                umi_str = ids.join("");
                umi_start = umi_match.start();
                umi_len = umi_match.len();
            } else if captures.len() == 2 {
                // Only one capture group in the regex
                let umi_match = captures.get(1).unwrap();
                umi_str = umi_match.as_str().to_string();
                umi_start = umi_match.start();
                umi_len = umi_match.len();
            } else {
                // No UMI found in this sequence
                continue;
            };

            // Try to find the index information in the record description
            let indexes = match DUAL_INDEX_REGEX.find(record.desc().unwrap_or_default()) {
                Some(matched) => matched.as_str(),
                None => match SINGLE_INDEX_REGEX.find(record.desc().unwrap_or_default()) {
                    Some(matched) => matched.as_str(),
                    None => panic!("Could not find index"),
                },
            };

            // Create a new description with the UMI appended
            let read_info_with_umi =
                format!("{}{}{} {}", record.id(), umi_delineator, umi_str, indexes);

            // Check if the UMI is at a valid position within the sequence and quality strings
            if umi_start <= read_sequence.len() && umi_start <= record.qual().len() {
                // Trim the UMI from the sequence and quality strings
                let umi_trimmed_sequence = &record.seq()[umi_len..];
                let trimmed_quality = &record.qual()[umi_len..];

                // Create a new record with the UMI-trimmed sequence and quality
                let new_record = Record::with_attrs(
                    &read_info_with_umi,
                    None,
                    umi_trimmed_sequence,
                    trimmed_quality,
                );

                // Write the new record to the output file
                if let Err(err) = writer.write_record(&new_record) {
                    eprintln!("Error writing record: {}", err);
                }
            }
        }
    }

    // Return Ok if everything went well
    Ok(())
}

/// Extracts Unique Molecular Identifiers (UMIs) from a given gzipped fastq file and writes the
/// UMI-trimmed sequences to a new uncompressed fastq file. The UMIs are searched using a regex
/// pattern related to known adapter sequences and are appended to the 3' adapter.
///
/// # Arguments
/// * `input_fastq` - The path to the input gzipped fastq file.
/// * `sample_name` - The name of the sample.
/// * `umi_delineator` - Delimiter used between appended UMI and read name. The
/// default is "_" as is used in UMI-tools.
///
/// # Example
/// ```rust
/// extract_umis_qiagen("sample.unprocessed.fastq.gz", "sample1", "_", true);
/// // This will create an output file named "sample1.processed.fastq" with UMI-trimmed sequences.
/// ```
///
/// Note: The UMI is extracted using the regex pattern "AACTGTAGGCACCATCAAT(.{12})AGATCGGAAG"
/// where the UMI is the 12 bases following the adapter sequence "AACTGTAGGCACCATCAAT" and before
/// the sequence "AGATCGGAAG". This is different from the previous version where UMIs were
/// expected to be the first 12 characters of the sequence.
pub fn extract_umis_qiagen(
    input_fastq: &str,
    sample_name: &str,
    umi_delineator: &str,
) -> std::io::Result<()> {
    let input_file = File::open(input_fastq)?;

    // Check if the file is gzipped or not
    let reader: Box<dyn Read> = if is_gzipped(input_fastq)? {
        Box::new(MultiGzDecoder::new(BufReader::new(input_file)))
    } else {
        Box::new(BufReader::new(input_file))
    };

    let mut writer = Writer::new(File::create(format!("{}.processed.fastq", sample_name))?);

    for record_result in Reader::new(reader).records() {
        let record = match record_result {
            Ok(record) => record,
            Err(err) => {
                eprintln!("Error reading record: {}", err);
                continue;
            }
        };

        let read_sequence = match String::from_utf8(record.seq().to_vec()) {
            Ok(s) => s,
            Err(err) => {
                eprintln!("Error decoding sequence as UTF-8: {}", err);
                continue;
            }
        };

        if let Some(captures) = UMI_REGEX_QIAGEN.captures(&read_sequence) {
            if let Some(umi_match) = captures.get(1) {
                let umi_str = umi_match.as_str();
                let umi_start = umi_match.start();

                let indexes = match DUAL_INDEX_REGEX.find(record.desc().unwrap_or_default()) {
                    Some(matched) => matched.as_str(),
                    None => match SINGLE_INDEX_REGEX.find(record.desc().unwrap_or_default()) {
                        Some(matched) => matched.as_str(),
                        None => panic!("Could not find index"),
                    },
                };

                let read_info_with_umi =
                    format!("{}{}{} {}", record.id(), umi_delineator, umi_str, indexes);

                if umi_start <= read_sequence.len() && umi_start <= record.qual().len() {
                    let umi_trimmed_sequence = read_sequence[..umi_start].to_string();
                    let trimmed_quality = &record.qual()[..umi_start];

                    let new_record = Record::with_attrs(
                        &read_info_with_umi,
                        None,
                        umi_trimmed_sequence.as_bytes(),
                        trimmed_quality,
                    );

                    if let Err(err) = writer.write_record(&new_record) {
                        eprintln!("Error writing record: {}", err);
                    }
                }
            }
        }
    }

    Ok(())
}
