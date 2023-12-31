use crate::extract_umis;
use crate::Command;
use crate::HashMap;
use crate::Path;
use crate::ProgressBar;
use crate::ProgressStyle;
use crate::Regex;

/// Trim the 3' adapter from each read.
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
    let progress_br = ProgressBar::new(num_of_fastqs);

    progress_br.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:50.cyan/blue}] {pos}/{len} {msg} ({percent}%)")
                .expect("Progress bar error")
            .progress_chars("#>-"),
    );
    progress_br.set_message("Trimming fastq files...");

    let mut trimmed_fastq_files: Vec<String> = vec![];
    for fastq_file in fastqs.iter() {
        let sample_name = fastq_file.split('_').next().unwrap();

        // Create the strings before passing them to the vec
        let min_length_string = minimum_length.to_string();
        let max_length_string = maximum_length.to_string();
        let too_short_output = format!("{}.short.fastq", sample_name);

        // Set cutadapt args to remove untrimmed reads and require a quality score > 25
        let mut cutadapt_args: Vec<&str> = vec![
            "--cores=0",
            "--discard-untrimmed",
            "--quality-cutoff",
            "25,25",
            "--adapter",
        ];

        // Determine order of UMI extraction and adapter trimming
        let umi_loc = library_type.get(&sample_name.to_string()).unwrap().as_str();
        match umi_loc {
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
