use crate::extract_umis;
use crate::create_progress_bar;
use crate::Command;
use crate::HashMap;
use crate::Path;
use crate::Regex;

/// Trim the 3' adapter from each read.
///
/// The bases on the 3' end are also quality trimmed if their quality score is below 20. Reads
/// shorter than set minimum or that weren't trimmed are discarded. UMIs are extracted from each
/// sequence and appended to the read name.
///
/// # Example
/// ```
/// use regex::Regex;
/// 
/// let fastqs = vec!["sample1_R1.fastq.gz","sample2_R1.fastq.gz","sample3_R1.fastq.gz"];
/// let mut library_type = HashMap::new();
/// library_type.insert("sample1", "standard");
/// library_type.insert("sample2", "standard");
/// library_type.insert("sample3", "standard");
///
/// let umi_regex = Regex::new(r"(^.{12})").unwrap();
///
/// let trimmed = trim_adapters(fastqs, library_type, 16, 30, umi_regex, &"ATCG".to_string(), false, true, false);
/// // trimmed is ["sample1.cut.fastq", "sample2.cut.fastq", "sample3.cut.fastq"]
/// ```
///
/// # Arguments
/// * `fastqs` - A vector of strings containing the file names of the FASTQ files.
/// * `library_type` - A HashMap associating each sample name to the corresponding library type.
/// * `minimum_length` - The minimum sequence length required for a read to be kept.
/// * `maximum_length` - The maximum sequence length acceptable for a read to be kept.
/// * `umi_regex` - A regular expression used to match Unique Molecular Identifiers (UMIs) in the reads.
/// * `adapter` - The 3' adapter sequence.
/// * `is_qiagen` - A Boolean indicating whether the analysis is using Qiagen data.
/// * `is_3p` - A Boolean indicating whether the UMI is on the 3' end of each read.
/// * `mismatch` - A Boolean indicating if 1 bp mismatches in the adapter sequence are allowed.
///
/// # Returns
/// * A vector of strings containing the names of the FASTQ files after trimming.
pub fn trim_adapters(
    fastqs: Vec<String>,
    library_type: HashMap<String, String>,
    minimum_length: u8,
    maximum_length: u8,
    umi_regex: &Regex,
    adapter: &str,
    is_qiagen: bool,
    is_3p: bool,
    mismatch: bool
) -> Vec<String> {
    let num_of_fastqs = fastqs.len() as u64;
    let progress_br = create_progress_bar(num_of_fastqs, "Trimming fastq files...".to_string());

    let mut trimmed_fastq_files: Vec<String> = vec![];
    for fastq_file in fastqs.iter() {
        let sample_name = fastq_file.split('_').next().unwrap();

        // Create the strings before passing them to the vec
        let min_length_string = minimum_length.to_string();
        let max_length_string = maximum_length.to_string();
        let too_short_output = format!("{}.short.fastq", sample_name);

        // Set cutadapt args to remove untrimmed reads and trim bases with a quality score < 25
        let mut cutadapt_args: Vec<&str> = vec![
            "--cores=0",
            "--discard-untrimmed",
            "--quality-cutoff",
            "20",
            "--adapter",
        ];

        // Determine order of UMI extraction and adapter trimming
        let umi_type = library_type.get(&sample_name.to_string()).unwrap().as_str();
        match umi_type {
            "standard" => { // Standard meaning the UMI is located between the adapters
                let output_file_name = format!("{}.unprocessed.cut.fastq", sample_name);
                cutadapt_args.extend([
                    adapter,
                    "--output",
                    &output_file_name,
                    "--minimum-length",
                    &min_length_string,
                    "--maximum-length",
                    &max_length_string,
                    "--too-short-output",
                    &too_short_output,
                    fastq_file,
                ]);
    
                // Call cutadapt
                let output = Command::new("cutadapt").args(&cutadapt_args).output();
    
                match output {
                    Ok(_) => (),
                    Err(ref err) => eprintln!("Error running cutadapt: {:?}", err),
                }
    
                // Call function to extract UMIs from each read
                if let Some(extract) = Some(extract_umis) {
                    let args = (&output_file_name, sample_name, "_", umi_regex, adapter, is_qiagen, is_3p, mismatch);
                    if let Err(err) = extract(args.0, args.1, args.2, args.3, args.4, args.5, args.6, args.7) {
                        eprintln!("Error when extracting UMIs: {:?}", err);
                    }
                }
    
                let output_filename = format!("{}.processed.fastq", sample_name);
                if Path::new(&output_filename).exists() {
                    trimmed_fastq_files.push(output_filename);
                }
            }
            "qiagen" => {
                // Extraction must come first because trimming the 3' adapter first would remove the UMI
                if let Some(extract) = Some(extract_umis) {
                    let args = (fastq_file, sample_name, "_", umi_regex, adapter, is_qiagen, is_3p, mismatch);
                    if let Err(err) = extract(args.0, args.1, args.2, args.3, args.4, args.5, args.6, args.7) {
                        eprintln!("Error when extracting UMIs: {:?}", err);
                    }
                }
    
                let fastq = format!("{}.processed.fastq", sample_name);
                let output_file_name = format!("{}.cut.fastq", sample_name);
    
                cutadapt_args.extend([
                    adapter,
                    "--output",
                    &output_file_name,
                    "--minimum-length",
                    &min_length_string,
                    "--too-short-output",
                    &too_short_output,
                    &fastq,
                ]);
    
                let output = Command::new("cutadapt").args(&cutadapt_args).output();
    
                match output {
                    Ok(_) => (),
                    Err(ref err) => eprintln!("Error running cutadapt: {:?}", err),
                }
    
                let output_filename = format!("{}.cut.fastq", sample_name);
                if Path::new(&output_filename).exists() {
                    trimmed_fastq_files.push(output_filename);
                }
            }
            _ => {
                panic!("Unknown library type.");
            }
        } 
        progress_br.inc(1);
    }

    progress_br.finish_with_message("Finished trimming fastq files");

    // Sort trimmed fastq files
    trimmed_fastq_files.sort_unstable();

    trimmed_fastq_files
}
