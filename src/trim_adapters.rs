use crate::HashMap;
use crate::ProgressBar;
use crate::ProgressStyle;
use crate::Command;
use crate::Path;
use bio::io::fastq::{Reader, Writer, Record};
use crate::SINGLE_INDEX_REGEX;
use crate::DUAL_INDEX_REGEX;

/// Trim the 3' adapter from each read.
///
/// UMI Read 1 Setup:
/// `5' Adapter - UMI - Insert - 3' Adapter`
/// `AGATCGGAAGAGCGTCGTGTAGGGAAAGANNNNNNNNNNNN - TGTCAGTTTGTCAAATACCCCA - TGGAATTCTCGGGTGCCAAGG`
///
/// The bases on the 3' end are also quality trimmed if their quality score is below 20. Reads 
/// shorter than 16 bases or that weren't trimmed are discarded. UMIs are extracted from each
///  sequence and appended to the read name.
///
/// # Example
/// ```
/// let fastqs = vec!["sample1.fastq.gz","sample2.fastq.gz","sample3.fastq.gz"];
/// let mut library_type = HashMap::new();
/// library_type.insert("sample1", "UMI");
/// library_type.insert("sample2", "UMI");
/// library_type.insert("sample3", "UMI");
/// let minimum_length = 28;
///
/// let trimmed = trim_adapters(fastqs, library_type, minimum_length);
/// // trimmed is ["sample1.cut.fastq", "sample2.cut.fastq", "sample3.cut.fastq"]
/// ```
///
/// # Arguments
/// * `fastqs` - A vector of strings containing the file names of the FASTQ files.
/// * `library_type` - A HashMap associating each sample name to the corresponding library type.
/// * `minimum_length` - The minimum sequence length required for a read to be kept.
///
/// # Returns
/// * A vector of strings containing the names of the FASTQ files after trimming.
pub fn trim_adapters(fastqs: Vec<String>
                , library_type: HashMap<String, String>
                , minimum_length: u8)
                -> Vec<String> {
    let num_of_fastqs = fastqs.len() as u64;
    let progress_br = ProgressBar::new(num_of_fastqs);
    
    progress_br.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{eta_precise}] [{bar:50.cyan/blue}] {pos}/{len} {msg} ({percent}%)")
                .expect("Progress bar error")
            .progress_chars("#>-"),
    );
    progress_br.set_message("Trimming fastq files...");

    let mut trimmed_fastq_files: Vec<String> = vec![];
    for fastq_file in fastqs.iter() {
        let sample_name = fastq_file.split("_").next().unwrap();
        let adapter: &str;
        let output_file_name: String;
        let umi_extraction_function: Option<fn(&str, &str) -> std::io::Result<()>> = Some(extract_umis);
        let umi_extraction_arguments: Option<(&str, &str)> = Some((sample_name, "_"));

        // Set cutadapt args to remove untrimmed reads and require a quality score > 20
        let mut cutadapt_args: Vec<&str> = vec![
            "--cores=0",
            "--discard-untrimmed",
            "--quality-cutoff",
            "20",
            "--adapter",
        ];
    
        match library_type.get(&sample_name.to_string()).unwrap().as_str() {
            "UMI" => {
                adapter = "TGGAATTCTCGGGTGCCAAGG";
                output_file_name = format!("{}.unprocessed.cut.fastq", sample_name);    
            }
            _ => {
                panic!("Unknown library type.");
            }
        }        
    
        // Create the strings before passing them to the vec
        let min_length_string = minimum_length.to_string();
        let too_short_output = format!("{}.short.fastq", sample_name);
    
        cutadapt_args.extend([
            adapter,
            "--output",
            &output_file_name,
            "--minimum-length",
            &min_length_string,
            "--too-short-output",
            &too_short_output,
            fastq_file,
        ]);

        // Call cutadapt
        let output = Command::new("cutadapt")
            .args(&cutadapt_args)
            .output();

        match output {
            Ok(_) => (),
            Err(ref error) => eprintln!("Error running cutadapt: {:?}", error),
        }

        // Call function to extract UMIs from each read
        if let Some(extract) = umi_extraction_function {
            let args = umi_extraction_arguments.unwrap();
            if let Err(e) = extract(args.0, args.1) {  // handle Result
                eprintln!("Error when extracting UMIs: {:?}", e);
            }
        }

        let output_filename = format!("{}.cut.fastq", sample_name);
        if Path::new(&output_filename).exists() {
            trimmed_fastq_files.push(output_filename);
        }

        progress_br.inc(1);
    }

    progress_br.finish_with_message("Finished trimming fastq files");

    // Sort trimmed fastq files
    trimmed_fastq_files.sort_unstable();

    trimmed_fastq_files
}

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