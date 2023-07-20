use bio::io::fastq::{Reader, Writer, Record};
use std::io::BufReader;
use flate2::bufread::MultiGzDecoder;
use std::fs::File;
use crate::SINGLE_INDEX_REGEX;
use crate::DUAL_INDEX_REGEX;
use crate::UMI_REGEX_QIAGEN;
use crate::rna_counts::MyError;

/// Extract Unique Molecular Identifiers (UMIs) from a given fastq file and write the 
/// UMI-trimmed sequences to a new fastq file. The UMIs are assumed to be the first 12 characters
/// of the sequence.
///
/// # Arguments
/// * `sample_name` - The name of the sample.
/// * `umi_delineator` - Delimiter used between appended UMI and read name. The
/// default is "_" as is used in UMI-tools.
///
/// # Example
/// ```rust
/// extract_umis("sample1", "_");
/// // This will create an output file named "sample1.cut.fastq" with UMI-trimmed sequences.
/// ```
pub fn extract_umis(sample_name: &str, umi_delineator: &str) -> std::io::Result<()> {
    let filename = format!("{}.{}.{}", sample_name, "unprocessed.cut", "fastq");

    // Open the input and output files
    let reader = match Reader::from_file(&filename) {
        Ok(reader) => reader,
        Err(err) => return Err(std::io::Error::new(std::io::ErrorKind::Other, err.to_string())),
    };

    let mut writer = Writer::to_file(format!("{}.cut.fastq", sample_name))?;

    // Go through each sequence in the FASTQ file
    for record_result in reader.records() {
        let record = match record_result {
            Ok(record) => record,
            Err(err) => return Err(std::io::Error::new(std::io::ErrorKind::Other, err.to_string())),
        };
        
        // Capture read sequence
        let read_sequence = record.seq();

        // Capture UMI in sequence
        let umi = String::from_utf8(read_sequence[0..12].to_vec()).unwrap();

        // Capture index information
        let indexes = match DUAL_INDEX_REGEX.find(record.desc().unwrap_or_default()) {
            Some(matched) => matched.as_str(),
            None => match SINGLE_INDEX_REGEX.find(record.desc().unwrap_or_default()) {
                Some(matched) => matched.as_str(),
                None => panic!("Could not find index"),  
            },
        };

        // If a UMI was found, add it to the sequence description after an underscore
        let read_info_with_umi = format!("{}{}{} {}", record.id(), umi_delineator, umi, indexes);

        // Remove the UMI sequence from the read sequence and quality
        let umi_trimmed_sequence = String::from_utf8(read_sequence[12..].to_vec()).unwrap();
        let trimmed_quality = &record.qual()[12..];
        
        // Create a new fastq record with the UMI sequence moved
        let new_record = Record::with_attrs(
            &read_info_with_umi, 
            None, 
            umi_trimmed_sequence.as_bytes(), 
            trimmed_quality,
        );

        // Write the sequence to the output file
        writer.write_record(&new_record)?;
    }

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
pub fn extract_umis_qiagen(input_fastq: &str, sample_name: &str, umi_delineator: &str) -> Result<(), MyError> {
    let mut writer = Writer::to_file(format!("{}.processed.fastq", sample_name))?;

    let input_file = File::open(input_fastq)?;
    let reader = Reader::new(MultiGzDecoder::new(BufReader::new(input_file)));

    // Go through each sequence in the FASTQ file
    for record_result in reader.records() {
        // Check if the record is as expected
        let record = match record_result {
            Ok(ref r) => r,
            Err(e) => {
                println!("Error reading record: {}", e);
                continue;
            }
        };

        // Capture read sequence as a string
        let read_sequence = String::from_utf8(record.seq().to_vec()).unwrap();
        let _read_description = record.desc().unwrap_or_default();

        // Capture UMI appended to the 3' adapter using regex
        if let Some(umi_match) = UMI_REGEX_QIAGEN.captures(&read_sequence) {
            let umi = umi_match[1].to_string();
            let umi_start = umi_match.get(0).unwrap().start();
            let umi_end = umi_match.get(0).unwrap().end();

            // Capture index information
            let indexes = match DUAL_INDEX_REGEX.find(record.desc().unwrap_or_default()) {
                Some(matched) => matched.as_str(),
                None => match SINGLE_INDEX_REGEX.find(record.desc().unwrap_or_default()) {
                    Some(matched) => matched.as_str(),
                    None => panic!("Could not find index"),  
                },
            };

            // If a UMI was found, add it to the sequence description after an underscore
            let read_info_with_umi = format!("{}{}{} {}", record.id(), umi_delineator, umi, indexes);

            // Remove the UMI sequence from the read sequence and quality
            let umi_trimmed_sequence = format!("{}{}", read_sequence[0..umi_start].to_string()
                                                             , &read_sequence[umi_end..].to_string()
                                                            );
            let mut trimmed_quality = record.qual()[0..umi_start].to_vec();
            trimmed_quality.extend_from_slice(&record.qual()[umi_end..]);

            // Create a new fastq record with the UMI sequence moved
            let new_record = Record::with_attrs(
                &read_info_with_umi,
                None,
                umi_trimmed_sequence.as_bytes(),
                &trimmed_quality,
            );

            // Write the sequence to the output file
            writer.write_record(&new_record)?;
        } else {
            continue;
        }
    } 
    Ok(())
}