use bio::io::fastq::{Reader, Writer, Record};
use std::io::BufReader;
use flate2::bufread::MultiGzDecoder;
use std::fs::File;
use crate::SINGLE_INDEX_REGEX;
use crate::DUAL_INDEX_REGEX;
use crate::UMI_REGEX_QIAGEN;

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
pub fn extract_umis_qiagen(input_fastq: &str, sample_name: &str, umi_delineator: &str) -> std::io::Result<()> {
    
    let input_file = File::open(input_fastq)?;
    let reader = Reader::new(MultiGzDecoder::new(BufReader::new(input_file)));

    let mut writer = Writer::new(File::create(format!("{}.processed.fastq", sample_name))?);

    for record_result in reader.records() {
        let record = match record_result {
            Ok(record) => record,
            Err(err) => {
                eprintln!("Error reading record: {}", err);
                continue;
            },
        };

        let read_sequence = match String::from_utf8(record.seq().to_vec()) {
            Ok(s) => s,
            Err(err) => {
                eprintln!("Error decoding sequence as UTF-8: {}", err);
                continue;
            },
        };
        
        if let Some(captures) = UMI_REGEX_QIAGEN.captures(&read_sequence) {
            if let Some(umi_match) = captures.get(1) {
                let umi_str = umi_match.as_str();
                let umi_start = umi_match.start();

                if let Some(index_match) = DUAL_INDEX_REGEX.find(record.desc().unwrap_or_default()) {
                    let indexes = index_match.as_str();
                    let read_info_with_umi = format!("{}{}{} {}", record.id(), umi_delineator, umi_str, indexes);

                    if umi_start <= read_sequence.len() && umi_start <= record.qual().len() {
                        let umi_trimmed_sequence = format!("{}", &read_sequence[..umi_start]);
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
    }
    
    Ok(())
}
