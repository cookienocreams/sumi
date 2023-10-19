use crate::is_gzipped;
use crate::Read;
use crate::Regex;
use crate::DUAL_INDEX_REGEX;
use crate::SINGLE_INDEX_REGEX;
use crate::UMI_REGEX_QIAGEN;
use crate::Error;
use bio::io::fastq::{Reader, Record, Writer};
use flate2::bufread::MultiGzDecoder;
use std::fs::File;
use std::io::BufReader;
use bio::scores::blosum62;
use bio::alignment::pairwise::*;

/// Extract Unique Molecular Identifiers (UMIs) from a given fastq file.
/// 
/// It writes the UMI-trimmed sequences to a new fastq file. The function is highly flexible in 
/// defining the UMI, as users can provide a regular expression with one or more capture groups to 
/// specify the UMI's location and format within the sequence.The function can handle potential 
/// mismatches during alignment.
///
/// # Arguments
/// * `input_fastq` - The path to the input fastq file, which can be either gzipped or uncompressed.
/// * `sample_name` - The name of the sample.
/// * `umi_delineator` - Delimiter used between the appended UMI and the read name. The
/// default is "_" as is used in UMI-tools.
/// * `umi_regex` - A regular expression used to capture the UMI from the read sequence. The 
/// expression can contain multiple capture groups for flexibility in UMI definitions.
/// * `adapter` - The adapter sequence for alignment purposes.
/// * `is_qiagen` - A flag indicating a specific trimming behavior.
/// * `is_3p` - A flag indicating whether the UMI is on the 3' end of the sequence.
/// * `adapter_mismatch` - A flag indicating whether to allow mismatches and indels during alignment.
///
/// # Example
/// ```rust
/// use regex::Regex;
///
/// let umi_regex = Regex::new(r"(^.{12})").unwrap();
/// extract_umis("sample1.fastq", "sample1", "_", &umi_regex, &"ATCG".to_string(), false, true, false);
/// // This will create an output file named "sample1.processed.fastq" with UMI-trimmed sequences.
/// ```
pub fn extract_umis(
    input_fastq: &str,
    sample_name: &str,
    umi_delineator: &str,
    umi_regex: &Regex,
    adapter: &String,
    is_qiagen: bool,
    is_3p: bool,
    adapter_mismatch: bool
) -> Result<(), Box<dyn Error>> {
    let input_file = File::open(input_fastq)?;

    // Check if the file is gzipped or not
    let reader: Box<dyn Read> = if is_gzipped(input_fastq)? {
        Box::new(MultiGzDecoder::new(BufReader::new(input_file)))
    } else {
        Box::new(BufReader::new(input_file))
    };

    let mut writer = Writer::new(File::create(format!("{}.processed.fastq", sample_name))?);

    // Set alignment penalties for mismatch tolerant Qiagen adapter search
    let gap_open = -5;
    let gap_extend = -1;
    let qiagen_umi_regex = Regex::new(UMI_REGEX_QIAGEN.as_str()).unwrap();

    for record_result in Reader::new(reader).records() {
        let record = record_result?;

        let read_sequence = record.seq();

        // Skip aligning reads with a 5' UMI if adapter is already trimmed in an error tolerant way using cutadapt.
        if adapter_mismatch && is_qiagen {
            // Align the Qiagen 3' adapter to the read sequence with mismatches and indels allowed
            let mut aligner = 
            Aligner::with_capacity(
                read_sequence.len(), 
                adapter.len(), 
                gap_open, 
                gap_extend, 
                &blosum62
            );
            let alignment = aligner.local( read_sequence, adapter.as_bytes());

            let start = alignment.xend;
            let umi_containing_seq = &read_sequence[start..];

            let search_sequence = String::from_utf8_lossy(umi_containing_seq);

            // Extract UMIs using Regex
            regex_extraction(
                search_sequence.to_string(), record,
                umi_delineator, &mut writer,
                umi_regex, is_qiagen, is_3p, 
                start, adapter_mismatch
            )            
        } else if is_qiagen && !adapter_mismatch {
            let start = 0;
            regex_extraction(
                String::from_utf8_lossy(read_sequence).to_string(), record, 
                umi_delineator, &mut writer, 
                &qiagen_umi_regex, is_qiagen, is_3p, 
                start, adapter_mismatch
            )
        } else {
            let start = 0;
            regex_extraction(
                String::from_utf8_lossy(read_sequence).to_string(), record, 
                umi_delineator, &mut writer, 
                umi_regex, is_qiagen, is_3p, 
                start, adapter_mismatch
            )
        }
        
    }

    // Return Ok if everything went well
    Ok(())
}

/// Extract Unique Molecular Identifiers (UMIs) from a given sequence using regex.
///
/// The UMIs are searched using a regex
/// pattern related to known adapter sequences and can be located on either the 5' or 3' end. 
/// The function is able to process sequences from both gzipped and uncompressed FASTQ files.
///
/// # Arguments
/// * `read_sequence` - The DNA sequence from the FASTQ record.
/// * `record` - The FASTQ record from which the UMI and other details are extracted.
/// * `umi_delineator` - Delimiter used between appended UMI and read name. The
/// default is "_" as is used in UMI-tools.
/// * `writer` - A mutable reference to the FASTQ writer for writing out the UMI-trimmed sequences.
/// * `umi_regex` - The regex pattern used to identify and capture the UMI in the `read_sequence`.
/// * `is_qiagen` - A flag indicating a specific trimming behavior.
/// * `is_3p` - A flag indicating whether the UMI is on the 3' end of the sequence.
/// * `start` - The index where the regex first matches.
/// * `adapter_mismatch` - A flag indicating whether to allow mismatches and indels during alignment.
///
/// # Example
/// ```rust
/// let mut writer = Writer::new(File::create(format!("{}.processed.fastq", "sample_1")).unwrap());
/// 
/// regex_extraction("GCATGCTACTGCTCGAT", record, "_", &mut writer, &"(.{12}$)", false, true, 35, false);
/// // This will write to `writer` the UMI-trimmed sequences with UMI appended to the record's description.
/// ```
///
/// Note: The UMI is extracted using the provided `umi_regex`. This function is flexible to different 
/// UMI positions and patterns, unlike the previous version which expected UMIs to be the first 12 characters 
/// of the sequence.
///
/// It's important to ensure that the regex pattern provided accurately captures the UMI without unintentionally 
/// capturing other portions of the sequence.
pub fn regex_extraction (
    read_sequence: String, 
    record: Record,
    umi_delineator: &str,
    writer: &mut Writer<File>,
    umi_regex: &Regex,
    is_qiagen: bool,
    is_3p: bool,
    start: usize,
    adapter_mismatch: bool
) {
    // Check if the sequence matches the UMI regex
    if let Some(captures) = umi_regex.captures(&read_sequence) {
        let umi_str: String;
        let umi_len: usize;
        let umi_start: usize;

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
            umi_str = "".to_string(); 
            umi_start = 0;
            umi_len = 0;
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

        let umi_trimmed_sequence: Vec<u8>;
        let trimmed_quality: Vec<u8>;

        // Trim the UMI from the sequence and quality strings
        if is_qiagen && adapter_mismatch {
            umi_trimmed_sequence = record.seq()[..start].to_vec();
            trimmed_quality = record.qual()[..start].to_vec();
        } else if is_qiagen && !adapter_mismatch {
                umi_trimmed_sequence = read_sequence[..umi_start].as_bytes().to_vec();
                trimmed_quality = record.qual()[..umi_start].to_vec();
        } else if is_3p {
            umi_trimmed_sequence = read_sequence[..umi_start].as_bytes().to_vec();
            trimmed_quality = record.qual()[..umi_start].to_vec();
        } else {
            umi_trimmed_sequence = read_sequence[umi_len..].as_bytes().to_vec();
            trimmed_quality = record.qual()[umi_len..].to_vec();
        }

        let (umi_trimmed_sequence, trimmed_quality) = (umi_trimmed_sequence.as_slice(), trimmed_quality.as_slice());

        // Create new FASTQ record with the UMI removed
        let new_record = Record::with_attrs(
            &read_info_with_umi,
            None,
            umi_trimmed_sequence,
            trimmed_quality,
        );

        if let Err(err) = writer.write_record(&new_record) {
            eprintln!("Error writing record: {}", err);
        }
    }
}
