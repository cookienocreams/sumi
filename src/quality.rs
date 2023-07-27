use crate::File;
use crate::MultiGzDecoder;
use crate::BufReader;
use crate::ProgressBar;
use crate::ProgressStyle;
use crate::HashMap;
use crate::csv_writer;
use bio::io::fastq::Reader;

/// Calcualte the mean of an input vector
pub fn mean(numbers: &Vec<f32>) -> f32 {
    let sum: f32 = numbers.iter().sum();
    sum as f32 / numbers.len() as f32
}

/// Get the total error probability from a record's quality scores
///
/// This function accepts ASCII-encoded quality score and produces a scalar
/// value that estimates the 'total error probability' for that record.
/// It means that if you want to calculate average quality score, you don't just sum
/// up all the phred scores and find the average, but rather consider the 
/// probability distribution of the scores.
///
/// The Phred score for a base `Q` is calculated as `Q = ASCII value of quality score - 33`.
/// The error probability `P` for that base is then calculated as `P = 10^(-Q/10)`.
/// Then these probabilities are summed up for all the bases to derive total error.
/// 
/// # Arguments
///
/// `q_score` - The quality scores of the current record as an ASCII-encoded byte string.
///
/// # Returns
///
/// The total error probability for the record.
///
/// # References
///
/// - Illumina's explanation on quality scores and their corresponding error rates:
///   <https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/quality-scores.html>
pub fn get_q_score_probability(q_score: &[u8]) -> f32 {
    let mut probability_sum = 0.0;
    for &char in q_score.iter() {
        let phred = *&char as f32 - 33.0;
        let prob = 10.0_f32.powf(-phred / 10.0);
        probability_sum += prob
    }
    probability_sum
}

/// Calculate the average quality score of a FASTQ file
///
/// This function will calculate the average quality score for all records in a FASTQ file.
/// The function will calculate total error probability for a record from it's quality score 
/// using ` get_q_score_probability` function. This value represents the probability of error 
/// for each base, then divided by total number of sequences in fastq file. The per-base error 
/// rates are then transformed back to a Phred scale and averaged across all bases.
///
/// # Arguments
///
/// * `input_fastq` - Path of the FASTQ file to calculate average quality score.
///
/// # Returns
///
/// - The average Phred Quality Score of all bases in the FASTQ file
///
/// # Errors
/// 
/// This function will return an error if:
/// - The FASTQ file cannot be opened for any reason.
/// - There is an interruption while reading the file.
///
/// # Examples
/// ```
/// let avg_quality = average_read_quality("sample.fastq")?;
/// println!("The average quality score is {}.", avg_quality);
/// ```
pub fn average_read_quality(input_fastq: &str) -> Result<f32, Box<dyn std::error::Error>> {
    let file = File::open(input_fastq)?;
    let reader = Reader::new(MultiGzDecoder::new(BufReader::new(file)));

    let mut q_score_list = Vec::new();

    for record_result in reader.records() {
        let record = record_result?;
        let prob_q_score = get_q_score_probability(record.qual()) / record.seq().len() as f32;
        let q_score = -10.0 * prob_q_score.log10();
        q_score_list.push(q_score);
    }

    Ok(mean(&q_score_list))
}

/// Calculate the average quality scores for a list of FASTQ files.
///
/// # Arguments
///
/// * `fastq_files` - A slice containing the names of the FASTQ files to process.
///
/// # Returns
///
/// This function will return a vector of the average quality scores of the input FASTQ 
/// files. If an error occurs while reading any of the files, an error will be returned instead.
///
/// # Examples
///
/// ```
/// let avg_qualities = average_quality_scores(&["file1.fastq", "file2.fastq", "file3.fastq"])?;
/// println!("The average quality scores are: {:?}", avg_qualities);
/// ```
pub fn average_quality_scores(fastq_files: &Vec<String>) -> Result<(), Box<dyn std::error::Error>> {
    let num_of_fastqs = fastq_files.len() as u64;
    let progress_br = ProgressBar::new(num_of_fastqs);
    
    progress_br.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:50.cyan/blue}] {pos}/{len} {msg} ({percent}%)")
                .expect("Progress bar error")
            .progress_chars("#>-"),
    );
    progress_br.set_message("Calculating quality scores...");

    let mut quality_map = HashMap::new();
    for file in fastq_files {
        let avg_quality = average_read_quality(file)?;
        quality_map.insert(file.clone(), avg_quality);

        progress_br.inc(1);
    }

    let mut writer = csv_writer::from_path("q_scores.csv")?;
    writer.write_record(&["sample", "quality score"])?;

    for (sample, score) in &quality_map {
        writer.write_record(&[sample, &score.to_string()])?;
    }

    writer.flush()?;

    progress_br.finish_with_message("Finished calculating quality scores");

    Ok(())
}
