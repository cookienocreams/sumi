use crate::Config;
use crate::CsvWriter;
use crate::DataFrame;
use crate::Error;
use crate::File;
use crate::Path;
use crate::Series;
use crate::{HashMap, HashSet};
use crate::READ_NAME_UMI_REGEX;
use crate::find_true_isomir_umis;
use polars::prelude::*;
use bio::io::{fasta, fastq};
use std::hash::Hash;

// Create struct to hold the details of a canonical miRNA
pub struct MiRNA {
    mirna_name: String,
    canonical_sequence: String,
}

/// Count the occurrence of each isomiR and miRNA based on representative UMIs, and write the 
/// information to an output CSV file. The CSV file contains the name of each unique RNA, 
/// its count, and its reads per million (RPM) for the given sample.
///
/// The output CSV file will have the format `"{sample_name}_{RNA type}_counts.csv"` and will be 
/// located in the directory where the program is run.
///
/// # Returns
///
/// * A `Result` wrapping a `HashMap` where the keys are RNA names and the values are their counts.
///   If an error occurred during processing, returns an `Err` with a description of the error.
///
/// # Arguments
///
/// * `representative_umis`: A `HashSet` containing UMIs deemed to be "true" UMIs.
/// * `sample_name`: The name of the sample being processed.
/// * `umis`: A `HashMap` where the keys are UMIs and the values are their corresponding sequences.
/// * `seqs_or_mirna_counts`: A `HashMap` where the keys are miRNA names and the values are 
/// their corresponding sequences or read counts.
/// * `are_isomirs`: Boolean flag to indicate whether isomiRs or miRNA are being analyzed.
///
/// # Example
///
/// ```rust
/// let representative_umis = get_representative_umis(...);  // Function to get the representative UMIs
/// let isomir_umis = get_isomir_umis(...);  // Function to get the UMIs and their isomiRs
/// let umi = "CATGCTAGCTCA";
/// let mirna_name = "hsa-miR-640";
/// 
/// let mut mirna_counts: HashMap<String, u64> = HashMap::new();
/// *mirna_counts.entry(mirna_name.to_string()).or_insert(0) += 1;
/// 
/// let isomir_counts = calculate_read_counts(representative_umis, &"sample1".to_string(), isomir_umis, mirna_counts, true);
/// match isomir_counts {
///     Ok(counts) => println!("Successfully generated counts for sample1: {:?}", counts),
///     Err(e) => eprintln!("Error occurred: {}", e),
/// }
/// ```
///
/// # Notes
///
/// The function deduplicates UMIs and calculates the count and RPM of each isomiR and miRNA, saving the result in a CSV file.
pub fn calculate_read_counts<V>(
    representative_umis: HashSet<Vec<u8>>,
    sample_name: &String,
    umis: HashMap<Vec<u8>, String>,
    seqs_or_mirna_counts: HashMap<String, V>,
    are_isomirs: bool
) -> Result<HashMap<String, u64>, Box<dyn std::error::Error>> 
where
    V: Eq + Hash + std::fmt::Display
    {
    // Initialize a HashMap to count the occurrence of each RNA
    let mut dedup_counts: HashMap<String, u64> = HashMap::new();

    // Deduplicate UMIs and count number of isomiRs
    for (umi, rna) in umis {
        if representative_umis.contains(&umi) {
            *dedup_counts.entry(rna).or_insert(0) += 1;
        }
    }

    // Calculate the total number of mapped reads
    let total_mapped_reads: u64 = dedup_counts.values().sum();

    // Prepare vectors for DataFrame creation
    let mut names: Vec<String> = Vec::new();
    let mut counts: Vec<u64> = Vec::new();
    let mut rpms: Vec<f64> = Vec::new();
    let mut seqs: Vec<String> = Vec::new();

    // For each RNA, calculate its RPM and store its name, count, and RPM
    if are_isomirs {
        for (name, count) in &dedup_counts {
            names.push(name.to_string());
            counts.push(*count);
            rpms.push((*count as f64 / total_mapped_reads as f64) * 1_000_000.0); // Calculate RPM
            seqs.push(seqs_or_mirna_counts.get(name).unwrap().to_string())
        }
    } else {
        for (name, count) in &dedup_counts {
            names.push(name.to_string());
            counts.push(*count);
            rpms.push((*count as f64 / total_mapped_reads as f64) * 1_000_000.0);
        }
    }

    // Create a DataFrame with 'name', 'count', and 'RPM' as columns
    let mut counts_df = 
    if are_isomirs {
        DataFrame::new(vec![
            Series::new("name", names),
            Series::new("sequence", seqs), // Add sequence if analyzing isomiRs
            Series::new("count", counts),
            Series::new("RPM", rpms),
        ])?
    } else {
        DataFrame::new(vec![
            Series::new("name", names),
            Series::new("count", counts),
            Series::new("RPM", rpms),
        ])?
    };

    // Sort the DataFrame in-place by the 'count' column
    let _ = counts_df.sort_in_place(["count"], true, true);

    // Set RNA type for filename
    let rna_type = 
    if are_isomirs {
        "isomiR"
    } else {
        "miRNA"
    };

    // Create a new CSV file to store the DataFrame
    let file = File::create(format!("{}_{}_counts.csv", sample_name, rna_type))
        .expect("Could not create file");

    // Write the DataFrame to the CSV file with headers
    CsvWriter::new(file)
        .has_header(true)
        .finish(&mut counts_df)
        .unwrap();

    Ok(dedup_counts)
}

/// Generate a set of potential isomiRs derived from a given miRNA sequence.
///
/// The function produces isomiRs by truncating and extending both the 5' and 3' ends of the miRNA.
/// It creates isomiRs based on a buffer size (number of bases) specified for each end.
///
/// # Arguments
///
/// * `mirna`: The reference microRNA sequence, e.g., "TCTTTGGTTATCTAGCTGTATGA".
/// * `front_buffer_size`: Maximum number of bases to be considered for truncation and extension at the 5' end.
/// * `end_buffer_size`: Maximum number of bases to be considered for truncation and extension at the 3' end.
///
/// # Returns
///
/// Returns a `HashSet<String>` containing potential isomiR sequences derived from the miRNA.
///
/// # Examples
///
/// ```rust
/// let mirna = "TCTTTGGTTATCTAGCTGTATGA".to_string();
/// let potential_isomirs = create_isomirs(mirna, 2, 2);
/// assert!(potential_isomirs.contains("ATCTTTGGTTATCTAGCTGTATGA"));  // 5' extension example
/// assert!(potential_isomirs.contains("TCTTTGGTTATCTAGCTGTAT"));  // 3' truncation example
/// ```
pub fn create_isomirs(
    mirna: String, 
    front_buffer_size: usize, 
    end_buffer_size: usize
) -> HashSet<String> {
    let bases = vec!['A', 'T', 'C', 'G'];
    let mut isomirs: HashSet<String> = HashSet::new();

    // For each potential modification at the 5' end
    for i in 0..=front_buffer_size {
        let trimmed_5p_sequence = &mirna[i..];

        // Add 5' truncated isomiRs
        isomirs.insert(trimmed_5p_sequence.to_string());

        // Extend the 5' truncated sequence by adding bases at its 3' end
        for j in 0..=end_buffer_size {
            if j > 0 {
                let trimmed_3_sequence = &trimmed_5p_sequence[..trimmed_5p_sequence.len()-j];
                isomirs.insert(trimmed_3_sequence.to_string());
            }
            for base in &bases {
                let extended_sequence = format!("{}{}", trimmed_5p_sequence, base);
                isomirs.insert(extended_sequence);
            }
        }

        // If the current 5' modification is an addition of a base
        if i > 0 {
            for base in &bases {
                let extended_5_sequence = format!("{}{}", base, trimmed_5p_sequence);
                isomirs.insert(extended_5_sequence);
            }
        }
    }
    // Remove the reference miRNA sequence
    isomirs.remove(&mirna);

    isomirs
}

/// Create a `HashMap` representation of an RNA sequence database from a given FASTA file.
///
/// This function reads the provided FASTA file and maps each RNA sequence to its reference 
/// name (ID) from the FASTA header. This allows for efficient sequence lookup and comparison 
/// without the need for sequence alignment, making it especially useful for tracking truncated 
/// isoforms that might still align with traditional aligners like bowtie2.
///
/// # Returns
///
/// * A `Result` wrapping a `HashMap` where the keys are RNA sequences and the values are 
///   their corresponding reference names from the FASTA file. 
///   If an error occurred during processing, returns an `Err` with a description of the error.
///
/// # Arguments
///
/// * `fasta`: Path to the input FASTA file containing RNA sequences.
///
/// # Example
///
/// ```rust
/// let fasta_file = "rnas.fasta";
///
/// let rna_hashmap = create_rna_hashmap_from_fasta(&fasta_file);
/// match rna_hashmap {
///     Ok(hm) => println!("Successfully created RNA hashmap: {:?}", hm),
///     Err(e) => eprintln!("Error occurred: {}", e),
/// }
/// ```
pub fn create_rna_hashmap_from_fasta(
    fasta: &str, 
    config: &Config
) -> Result<(HashMap<String, String>, HashMap<String, String>) , Box<dyn Error>> {
    // Open the RNA fasta file
    let mut rna_file = fasta::Reader::from_file(Path::new(fasta))?.records();
    let mut mirna_hm: HashMap<String, String> = HashMap::new();
    let mut mirna_hm_mismatch: HashMap<String, String> = HashMap::new();

    while let Some(Ok(record)) = rna_file.next() {
        // Convert the byte sequence to a string
        let read_sequence = std::str::from_utf8(record.seq())?;
        let ref_name = record.id();

        mirna_hm.insert(read_sequence.to_string(), ref_name.to_string());
        mirna_hm_mismatch.insert(read_sequence.to_string(), ref_name.to_string());

        // Add 1 bp mismatched sequences if desired
        if let true = config.mismatch {
            let sequences = sequences_with_hamming_distance_of_1(read_sequence);
            for seq in sequences.iter() {
                mirna_hm_mismatch.insert(seq.to_string(), ref_name.to_string());
            }   
        }
    }

    Ok((mirna_hm, mirna_hm_mismatch))
}

/// Generate all DNA sequences with a Hamming distance of 1 from the given sequence.
///
/// # Arguments
///
/// * `seq` - A string slice that holds the DNA sequence. It should contain only the 
/// characters 'A', 'C', 'G', and 'T'.
///
/// # Returns
///
/// Returns a `Vec<String>` containing all the sequences with a Hamming distance of 1 
/// from the input sequence.
///
/// # Examples
///
/// ```
/// let seq = "CTATACAATCTACTGTCTTTC";
/// let sequences = sequences_with_hamming_distance_of_1(seq);
/// for sequence in sequences {
///     println!("{}", sequence);
/// }
/// ```
pub fn sequences_with_hamming_distance_of_1(seq: &str) -> Vec<String> {
    let mut sequences = Vec::new();
    
    for (index, nucleotide) in seq.chars().enumerate() {
        for &replacement in ['A', 'C', 'G', 'T'].iter() {
            if nucleotide != replacement {
                let mut new_seq = seq.to_string();
                new_seq.replace_range(index..index + 1, &replacement.to_string());
                sequences.push(new_seq);
            }
        }
    }
    
    sequences
}

/// Analyze a given fastq file to count canonical miRNA and isomiR sequences.
///
/// This function first determines potential isomiRs of each canonical miRNA in the `mirna_hm` 
/// using the `create_isomirs` function. For each read in the fastq file, it determines if the
/// sequence is a canonical miRNA or an isomiR. For each detected sequence, it extracts the
/// UMI from the read name and, if the sequence is an isomiR, determines the corresponding
/// canonical miRNA and the specific modifications. The function then error-corrects the UMIs 
/// and counts deduplicated isomiRs based on these representative UMIs.
///
/// # Returns
///
/// * A `Result` wrapping a `HashMap` where the keys are isomiR names and the values are their deduplicated counts.
///   If an error occurred during processing, returns an `Err` with a description of the error.
///
/// # Arguments
///
/// * `fastq_file`: Path to the input fastq file.
/// * `sample_name`: Name of the sample being processed.
/// * `config`: Configuration settings, including the Levenshtein distance for UMI error correction.
/// * `mirna_hm`: A `HashMap` where the keys are canonical miRNA sequences and the values are their names.
/// * `max_diffs`: Maximum number of differences allowed when generating potential isomiRs.
///
/// # Notes
///
/// UMIs are extracted from the read names using a regex (`READ_NAME_UMI_REGEX`).
/// The error correction of UMIs is done using the `find_true_isomir_umis` function.
///
/// # Example
///
/// ```rust
/// let fastq_file = "sample1.fastq";
/// let config = Config {
///     levenshtein_distance: 2,
///     // ... other fields ...
/// };
/// let mirna_hm = get_canonical_mirnas();  // Function to get canonical miRNAs
/// 
/// let deduplicated_counts = isomir_analysis(&fastq_file, &"sample1".to_string(), &config, mirna_hm, 2);
/// match deduplicated_counts {
///     Ok(counts) => println!("Successfully generated deduplicated counts for sample1: {:?}", counts),
///     Err(e) => eprintln!("Error occurred: {}", e),
/// }
/// ```
pub fn isomir_analysis(
    fastq_file: &String, 
    sample_name: &String, 
    config: &Config, 
    mirna_hm: HashMap<String, String>,
    mirna_hm_mismatch: HashMap<String, String>,
    max_diffs: usize
) -> Result<(HashMap<String, u64>, HashMap<String, u64>) , Box<dyn Error>> {
    // Open the RNA fastq file
    let mut file = fastq::Reader::from_file(Path::new(fastq_file))?.records();

    let mut mirna_counts: HashMap<String, u64> = HashMap::new();
    let mut mirna_umi_counts: HashMap<Vec<u8>, u32> = HashMap::new();
    let mut mirna_umis: HashMap<Vec<u8>, String> = HashMap::new();

    let mut isomir_info: HashMap<String, MiRNA> = HashMap::new();
    let mut isomir_counts: HashMap<String, u64> = HashMap::new();
    let mut isomir_umis: HashMap<Vec<u8>, String> = HashMap::new();
    let mut isomir_umi_counts: HashMap<Vec<u8>, u32> = HashMap::new();
    let mut isomir_seqs: HashMap<String, String> = HashMap::new();

    // Create potential isomiRs HashMap to compare against each read sequence
    for (mirna_seq, mirna_name) in &mirna_hm {
        let potential_isomirs = create_isomirs(mirna_seq.to_string(), max_diffs.min(2), max_diffs);

        for isomir in &potential_isomirs {

            let mirna_entry = MiRNA {
                mirna_name: mirna_name.to_string(),
                canonical_sequence: mirna_seq.to_string(),
            };
            isomir_info.insert(isomir.to_string(), mirna_entry);
        }
    }

    // Get a reference to the appropriate miRNA HashMap
    let selected_hm = 
        if config.mismatch {
            &mirna_hm_mismatch
        } else {
            &mirna_hm
        };

    while let Some(Ok(record)) = file.next() {
        // Convert the byte sequence to a string
        let read_sequence = std::str::from_utf8(record.seq())?;
        let read_name = record.id();
        
        // Check if the sequence is a canonical miRNA
        if let Some(mirna_name) = selected_hm.get(read_sequence) {
            // Get UMI sequence from the read name
            if let Some(captures) = READ_NAME_UMI_REGEX.captures(read_name) {
                let umi = captures.get(1).unwrap().as_str();
                *mirna_umi_counts.entry(umi.as_bytes().to_vec()).or_insert(0) += 1;
                mirna_umis.insert(umi.as_bytes().to_vec(), mirna_name.to_string());
            }

            *mirna_counts.entry(mirna_name.to_string()).or_insert(0) += 1;
        } 
        // Check if the sequence is an isomiR
        if let Some((isomir, mirna_entry)) = isomir_info.get_key_value(read_sequence) {
            // Define isomiR name based on additions or deletions
            let isomir_name = 
                get_isomir_name(
                mirna_entry.canonical_sequence.clone().to_string(), 
                isomir.to_string(), 
                mirna_entry.mirna_name.to_string(), 
                max_diffs
            );

            // Collect all isomiR sequences
            isomir_seqs.insert(isomir_name.clone(), isomir.to_string());

            if let Some(captures) = READ_NAME_UMI_REGEX.captures(read_name) {
                let umi = captures.get(1).unwrap().as_str();
                *isomir_umi_counts.entry(umi.as_bytes().to_vec()).or_insert(0) += 1;
                isomir_umis.insert(umi.as_bytes().to_vec(), isomir_name.to_string());
            }

            *isomir_counts.entry(isomir_name.to_string()).or_insert(0) += 1;
        }
    }

    // Error correct UMIs in sample bam file
    let representative_isomir_umis: HashSet<Vec<u8>> = find_true_isomir_umis(
        isomir_umi_counts.clone(),
        config.levenshtein_distance,
    )?;

    // Error correct UMIs in sample bam file
    let representative_mirna_umis: HashSet<Vec<u8>> = find_true_isomir_umis(
        mirna_umi_counts.clone(),
        config.levenshtein_distance,
    )?;

    // Deduplicate mapped isomirs
    let dedup_isomir_counts = 
    match calculate_read_counts(
        representative_isomir_umis,
        sample_name,
        isomir_umis,
        isomir_seqs.clone(),
        true
    ) {
        Ok(rna_names) => rna_names,
        Err(e) => panic!("{}", e)
    };

    // Deduplicate mapped mirnas
    let dedup_mirna_counts = 
    match calculate_read_counts(
        representative_mirna_umis,
        sample_name,
        mirna_umis,
        mirna_counts,
        false
    ) {
        Ok(rna_names) => rna_names,
        Err(e) => panic!("{}", e)
    };

    Ok((dedup_isomir_counts, dedup_mirna_counts))
}

/// Construct the name of an isomiR based on the given base name and detected modifications.
///
/// This function applies a naming convention to represent variations in isomiRs, 
/// such as 5' or 3' additions and deletions.
///
/// The naming convention follows the format specified in: 
/// <https://github.com/miRTop/incubator/blob/master/isomirs/isomir_naming.md>
///
/// # Arguments
/// 
/// * `base_name`: The base name of the microRNA, e.g., "hsa-miR-9-5p".
/// * `modifications`: A HashMap containing detected modifications in the 
/// format {"type of modification": "modification detail"}.
/// 
/// # Returns
/// 
/// Returns a `String` representing the complete name of the isomiR with modifications.
///
/// # Examples
///
/// ```rust
/// let mut modifications = HashMap::new();
/// modifications.insert("5' addition".to_string(), "AA".to_string());
/// let isomir_name = construct_isomir_name("hsa-miR-9-5p", &modifications);
/// assert_eq!(isomir_name, "hsa-miR-9-5p.AAs");
/// ```
pub fn construct_isomir_name(base_name: &str, modifications: &HashMap<String, String>) -> String {
    let mut name = base_name.to_string();

    // Use conditions and string concatenation to craft isomiR name
    if modifications.contains_key("5' addition") {
        name += &format!(".{}s", modifications.get("5' addition").unwrap());
    }
    if modifications.contains_key("5' deletion") {
        name += &format!(".{}s", modifications.get("5' deletion").unwrap().to_lowercase());
    }
    if modifications.contains_key("3' addition") {
        name += &format!(".{}", modifications.get("3' addition").unwrap());
    }
    if modifications.contains_key("3' non-template addition") {
        name += &format!(".{}e", modifications.get("3' non-template addition").unwrap());
    }
    if modifications.contains_key("3' deletion") {
        name += &format!(".{}", modifications.get("3' deletion").unwrap().to_lowercase());
    }
    if modifications.contains_key("5' mutation") {
        name += &format!(".{}1{}", 
            modifications.get("5' mutation").unwrap().chars().last().unwrap(), 
            modifications.get("5' mutation").unwrap().chars().next().unwrap()
        );
    }
    if modifications.contains_key("3' mutation") {
        name += &format!(".{}{}{}{}", 
            modifications.get("3' mutation").unwrap().chars().last().unwrap(), 
            modifications.get("3' mutation").unwrap().chars().nth(1).unwrap(),
            modifications.get("3' mutation").unwrap().chars().nth(2).unwrap(),
            modifications.get("3' mutation").unwrap().chars().next().unwrap()
        );
    }

    name
}

/// Retrieve the name of an isomiR based on the provided miRNA sequence, isomiR sequence, and a base miRNA name.
///
/// This function applies a naming convention to represent variations in isomiRs, such as 5' or 3' additions and deletions, 
/// and generates a full name of the isomiR considering the detected modifications.
///
/// The naming convention follows the format specified in: 
/// <https://github.com/miRTop/incubator/blob/master/isomirs/isomir_naming.md>
///
/// # Arguments
/// 
/// * `mirna`: The reference microRNA sequence, e.g., "TCTTTGGTTATCTAGCTGTATGA".
/// * `isomir`: The sequence of the isomiR under study.
/// * `mirna_name`: The base name of the microRNA, e.g., "hsa-miR-9-5p".
/// * `max_diffs`: Maximum number of allowed differences at the 5' or 3' end of the sequences for detection.
/// 
/// # Returns
/// 
/// Returns a `String` representing the complete name of the isomiR considering the modifications detected.
///
/// # Examples
///
/// ```rust
/// let mirna = "TCTTTGGTTATCTAGCTGTATGA".to_string();
/// let isomir = "AATCTTTGGTTATCTAGCTGTATGAAA".to_string();
/// let name = "hsa-miR-9-5p".to_string();
/// let isomir_name = get_isomir_name(mirna, isomir, name, 2);
/// assert_eq!(isomir_name, "hsa-miR-9-5p.AAs.AA");
/// ```
pub fn get_isomir_name(
    mirna: String, 
    isomir: String, 
    mirna_name: String, 
    max_diffs: usize
) -> String {
    let mirna_len = mirna.len();
    let isomir_len = isomir.len();

    let mut modifications: HashMap<String, String> = HashMap::new();

    // When mirna is longer than isomir
    if mirna_len > isomir_len {
        // If mirna has a deleted portion at the 5' end
        if mirna[max_diffs..] == isomir[mirna[max_diffs..].len().abs_diff(isomir_len)..] {
            modifications.insert("5' deletion".to_string(), mirna[..mirna_len - isomir_len].to_string());
        }
        // If mirna has a deleted portion at the 3' end
        else if mirna[..mirna_len - (mirna_len - isomir_len)] == isomir {
            modifications.insert("3' deletion".to_string(), mirna[mirna_len - (mirna_len - isomir_len)..].to_string());
        }
        // Case of simultaneous 5' and 3' modifications
        else {
            for i in 0..=max_diffs {
                if mirna[i..].starts_with(&isomir[..isomir.len() - i]) {
                    modifications.insert("5' deletion".to_string(), mirna[0..i].to_string());
                }
                if mirna[max_diffs..mirna_len-i].ends_with(&isomir[i..]) {
                    modifications.insert("3' deletion".to_string(), mirna[mirna_len - i..].to_string());
                }
                if isomir[i..].starts_with(&mirna[..isomir.len() - i]) {
                    modifications.insert("5' addition".to_string(), isomir[0..i].to_string());
                }
            }
        }
    }
    // When isomir is longer than mirna
    else if isomir_len > mirna_len {
        // If isomir has an added portion at the 5' end
        if isomir[max_diffs..] == mirna {
            modifications.insert("5' addition".to_string(), isomir[..isomir_len - mirna_len].to_string());
        }
        // If isomir has an added portion at the 3' end
        else if isomir[..isomir_len - (isomir_len - mirna_len)] == mirna {
            // Check if additional bases are from the template
            let last_base = mirna.chars().last().unwrap();
            let extra_bases = &isomir[isomir_len - (isomir_len - mirna_len)..];
            let template_base: Vec<bool> = extra_bases.chars().map(|base| base == last_base).collect();
            if template_base.iter().any(|&base| !base) {
                modifications.insert("3' non-template addition".to_string(), isomir[isomir_len - (isomir_len - mirna_len)..].to_string());
            } else {
                modifications.insert("3' addition".to_string(), isomir[isomir_len - (isomir_len - mirna_len)..].to_string());
            }
        }
        else {
            for i in 0..=max_diffs {
                if isomir[i..].starts_with(&mirna[..mirna.len() - i]) {
                    modifications.insert("5' addition".to_string(), isomir[0..i].to_string());
                }
                if isomir[max_diffs..isomir_len-i].ends_with(&mirna[i..]) {
                    let last_base = mirna.chars().last().unwrap();
                    let extra_bases = &isomir[isomir_len - (isomir_len - mirna_len)..];
                    let template_base: Vec<bool> = extra_bases.chars().map(|base| base == last_base).collect();
                    if template_base.iter().any(|&base| !base) {
                        modifications.insert("3' non-template addition".to_string(), isomir[isomir_len - i..].to_string());
                    } else {
                        modifications.insert("3' addition".to_string(), isomir[isomir_len - i..].to_string());
                    }
                }
                if isomir[max_diffs..isomir_len].ends_with(&mirna[..mirna.len() - i]) {
                    modifications.insert("3' deletion".to_string(), mirna[0..i].to_string());
                }
            }
        }
    } // When isomir and mirna are the same length
    else if isomir_len == mirna_len {
        seed_and_extend(&mirna, &isomir, 12, &mut modifications)
    }

    construct_isomir_name(&mirna_name, &modifications)
}

/// Return a seed of specified length from the center of the reference sequence.
/// If the seed length is greater than the reference sequence length, the entire sequence is returned.
///
/// # Arguments
///
/// * `ref_seq` - A slice containing the reference sequence.
/// * `seed_len` - The desired length of the seed.
///
/// # Returns
///
/// * A slice representing the central part of the sequence with the given seed length.
pub fn get_initial_seed(ref_seq: &str, seed_len: usize) -> &str {
    let len = ref_seq.len();
    if len < seed_len {
        return ref_seq;
    }
    let start = (len - seed_len) / 2;

    &ref_seq[start..start + seed_len]
}

/// Use a seed and extend approach to find modifications on both ends (5' and 3') of isomiRs compared to the reference miRNA.
///
/// The seed and extend algorithm works as follows:
/// 1. A central seed is extracted from the reference sequence based on the specified seed length.
/// 2. This seed is then searched for in the new sequence.
/// 3. Upon locating the seed in the new sequence, the function starts extending outward in both directions.
/// 4. During extension, it captures any mismatches, additions, or deletions at both ends.
///
/// It provides a way to quickly align two sequences based on a known similar region (the seed) 
/// and then explore differences outside this aligned region. This is useful for detecting single 
/// mismatches on the end of the reference sequence or single bp shifts, i.e., a 5' deletion
/// compensated by a 3' addition.
///
/// # Arguments
///
/// * `ref_seq` - A slice containing the reference sequence.
/// * `new_seq` - A slice containing the new sequence to be compared with the reference sequence.
/// * `seed_len` - The length of the seed to be used for the initial alignment.
/// * `modifications` - A mutable reference to a HashMap to store the modifications found.
///
/// # Returns
///
/// * None. This function updates the modifications HashMap directly with entries indicating the 
/// type of modification (addition, deletion, or end mutation) and its position (5' or 3' end).
pub fn seed_and_extend(
    ref_seq: &str, 
    new_seq: &str, 
    seed_len: usize, 
    modifications: &mut HashMap<String, String>
) {
    // Set initial seed sequence for modifications search
    let ref_seed = get_initial_seed(ref_seq, seed_len);

    // Keep track of the bases that differ bewteen the sequences
    let mut front_diff_bases = String::new();
    let mut rear_diff_bases = String::new();

    // Set initial seed sequence
    let mut seed = String::new();
    seed.push_str(ref_seed);

    // Check for a single mutation on the 5' or 3' end of the isomiR
    let end_mutation =
        ref_seq[0..seed_len] == new_seq[0..seed_len] || ref_seq[ref_seq.len() - seed_len..] == new_seq[new_seq.len() - seed_len..];

    // Find end mutation sequence if present
    if end_mutation && ref_seq[0..seed_len] == new_seq[0..seed_len] {
        // Add single 3' mutation (base:seq_len:base, e.g. A22T)
        rear_diff_bases.push(new_seq.chars().last().unwrap());
        rear_diff_bases.push(new_seq.len().to_string().chars().next().unwrap());
        rear_diff_bases.push(new_seq.len().to_string().chars().nth(1).unwrap());
        rear_diff_bases.push(ref_seq.chars().last().unwrap());
        modifications.insert("3' mutation".to_string(), rear_diff_bases.to_string());
    } else if ref_seq[ref_seq.len() - seed_len..] == new_seq[new_seq.len() - seed_len..] {
        // Add single 5' mutation (ref_base:mut_base, e.g. AT)
        front_diff_bases.push(new_seq.chars().next().unwrap());
        front_diff_bases.push(ref_seq.chars().next().unwrap());
        modifications.insert("5' mutation".to_string(), front_diff_bases.to_string());
    }

    if !end_mutation {
        if let Some(seed_pos_in_new_seq) = new_seq.find(ref_seed) {
            for i in 1..seed_pos_in_new_seq + 1 {
                // Extend towards the 5' end
                let next_seed_base = new_seq.chars().nth(seed_pos_in_new_seq - i).unwrap();
                seed.insert_str(0, &String::from(next_seed_base));
                
                // Clear string after each round
                front_diff_bases.clear();
                
                // Continue to search backwards and keep track of 5' bases
                if let Some(pos) = new_seq.find(&seed) {
                    front_diff_bases.push_str(&ref_seq[..(seed_pos_in_new_seq - pos).min(pos + 1)]);
                    front_diff_bases.push_str(&new_seq[..(seed_pos_in_new_seq - pos).min(pos + 1)]);
                }
            }
    
            // Reset seed sequence for 3' end search
            seed.clear();
    
            for i in 0..(ref_seq.len() - seed_pos_in_new_seq - seed_len) {
                // Extend towards the 3' end
                let next_seed_base = new_seq.chars().nth(seed_pos_in_new_seq + seed_len + i).unwrap();
                seed.insert_str(seed.len(), &String::from(next_seed_base));
                
                rear_diff_bases.clear();
                
                // Continue to search forwards and keep track of 3' bases
                if let Some(pos) = new_seq.find(&seed) {
                    rear_diff_bases.push_str(&ref_seq[pos + seed.len() - 1..]);
                    rear_diff_bases.push_str(&new_seq[pos + seed.len() - 1..]);
                }
            }
    
            // Set modification type
            if ref_seq[..seed_len].starts_with(&new_seq[1..seed_len]) {
                modifications.insert("5' addition".to_string(), front_diff_bases[1..=1].to_string());
            } else {
                modifications.insert("5' deletion".to_string(), front_diff_bases[0..1].to_string());
            }
            
            if ref_seq[ref_seq.len() - seed_len..].ends_with(&new_seq[new_seq.len() - seed_len..new_seq.len() - 1]) {
                modifications.insert("3' addition".to_string(), rear_diff_bases[1..=1].to_string());
            } else {
                modifications.insert("3' deletion".to_string(), rear_diff_bases[0..1].to_string());
            }
        }
    }
    
}
