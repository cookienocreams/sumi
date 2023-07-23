use std::collections::{HashMap, HashSet};
use crate::io;
use crate::Graph;
use crate::NodeIndex;
use crate::READ_NAME_UMI_REGEX;
use crate::Error;
use crate::bam::RecordWriter;
use crate::hamming;
use crate::levenshtein;
use petgraph::visit::Bfs;
use petgraph::visit::VisitMap;
use petgraph::visit::Visitable;

// Type alias for a graph of UMIs (represented as vectors of bytes) and a mapping from 
// NodeIndex to a tuple of UMI and its count.
type UmiGraph = (Graph<Vec<u8>, i32>, HashMap<NodeIndex, (Vec<u8>, i32)>);

/// Construct a substring index from a slice of UMIs.
///
/// This function creates an index that maps each unique substring of a UMI to a set of UMIs 
/// that contain that substring. The size of the substring is determined by `slice_size`.
///
/// The index is used to quickly find potential neighbor UMIs that may be within a certain 
/// edit distance, without having to compute the distance between all pairs of UMIs.
///
/// # Arguments
///
/// * `umis` - A slice of UMIs, where each UMI is a byte array.
/// * `slice_size` - The size of the substring to use for the index. Each UMI is divided into
///   chunks of this size (from left to right, and without overlap between chunks), and each 
///   chunk is added to the index.
///
/// # Returns
///
/// A hashmap where the keys are substrings of UMIs, and the values are sets of indices. Each
/// index in a set corresponds to a UMI that contains the substring.
///
/// # Example
///
/// ```rust
/// let umis: Vec<Vec<u8>> = vec![b"ACGT".to_vec(), b"TGCA".to_vec(), b"GTAC".to_vec(), b"CATG".to_vec()];
/// let slice_size: usize = 2;
/// let index = build_substring_index(&umis, slice_size);
/// ```
///
/// In the above example, the index would map "AC" to {0}, "CG" to {0}, "GT" to {0, 2}, "TG" to {1, 3}, 
/// "CA" to {1, 3}, "GA" to {2}, and "AT" to {3}.
pub fn build_substring_index(umis: &[Vec<u8>], slice_size: usize) -> HashMap<Vec<u8>, HashSet<usize>> {
    let mut index: HashMap<Vec<u8>, HashSet<usize>> = HashMap::new();

    for (i, umi) in umis.iter().enumerate() {
        for slice in umi.chunks(slice_size) {
            index.entry(slice.to_vec()).or_insert_with(HashSet::new).insert(i);
        }
    }

    index
}

/// Create a graph from a dictionary of UMIs and their counts.
///
/// This function constructs a graph where each node represents a unique UMI and each edge represents
/// a edit distance less than or equal to `max_edits` between two UMIs. The edge creation is based 
/// on the condition that the count of one UMI should be at least five times the count of the other 
/// UMI minus 1. If the condition is met, an edge is added between the two UMIs. 
///
/// The function implements a slightly modified directional graph algorithm that allows for a edit 
/// Distance of 1 between UMIs, see Fu, Y., et al, (2018). Elimination of PCR duplicates in RNA-seq 
/// and small RNA-seq using unique molecular identifiers. <https://doi.org/10.1186/s12864-018-4933-1>. 
/// The original algorithm uses a two fold count threshold. This version of the algorithm uses a 
/// substring indexing strategy to efficiently find potential neighbors for each UMI before computing 
/// the edit distance, hence improving the performance for large datasets.
///
/// Note that Levenshtein distance could be used here, but it generally takes longer to run and 
/// has a negligible effect on UMI deduplication.
///
/// # Arguments
///
/// * `umi_count_dict` - A hashmap where the key is the UMI (as a byte array) and the value is its count.
/// * `max_edits` - The maximum edit distance allowed to consider two UMIs as neighbors 
/// (i.e., to add an edge between them).
///
/// # Returns
///
/// A tuple consisting of the UMI graph and a hashmap containing the node attributes. The graph is of 
/// type `Graph<Vec<u8>, i32>`, and the hashmap key is a `NodeIndex` which points to a tuple of UMI 
/// and its count.
///
/// # Example
///
/// ```rust
/// let umi_count_dict: HashMap<Vec<u8>, i32> = HashMap::new();
/// // populate umi_count_dict
/// let max_edits: i32 = 1;
/// let (graph, node_attributes) = umi_graph(&umi_count_dict, max_edits);
/// ```
pub fn umi_graph(umi_count_dict: &HashMap<Vec<u8>, i32>, max_edits: i32, use_levenshtein: bool) -> UmiGraph {
    let mut graph: Graph<Vec<u8>, i32> = Graph::<Vec<u8>, i32>::new();
    let mut node_attributes: HashMap<NodeIndex, (Vec<u8>, i32)> = HashMap::new();

    // Collect all unique UMIs from the given UMI-count dictionary
    let umis: Vec<Vec<u8>> = umi_count_dict.keys().cloned().collect();

    // Create nodes in the graph for each UMI
    for umi in &umis {
        let count = *umi_count_dict.get(umi).unwrap(); // get count for the current UMI
        let node = graph.add_node(umi.clone()); // add node to the graph
        node_attributes.insert(node, (umi.clone(), count));
    }

    // Slice size for substring index is calculated based on UMI length and max edits allowed
    let umi_length = umis[0].len();
    let slice_size = umi_length / (max_edits as usize + 1);

    // Build a substring index from the UMIs
    let substring_index = build_substring_index(&umis, slice_size);

    // Iterate over each UMI to find potential neighbors and add edges
    for i in 0..umis.len() {
        
        // Split current UMI into slices
        let umi_slices: Vec<Vec<u8>> = umis[i].chunks(slice_size).map(|chunk| chunk.to_vec()).collect(); 
        let mut potential_neighbors = HashSet::new();

        // Add UMIs having common substrings with the current UMI to the potential neighbors set
        for umi_slice in umi_slices {
            if let Some(neighbors) = substring_index.get(&umi_slice) {
                potential_neighbors.extend(neighbors);
            }
        }

        // For each potential neighbor, compute edit distance and add edge if conditions are met
        for j in potential_neighbors {
            // Ensure that we don't compute the distance between a UMI and itself
            if i != j {
                let edit_distance =
                if use_levenshtein {
                    levenshtein(&umis[i], &umis[j]) as i32
                } else {
                    edit(&umis[i], &umis[j]) as i32
                };

                // If edit distance is within max edits, consider for edge addition
                if edit_distance <= max_edits {
                    let node_i: NodeIndex = NodeIndex::new(i);
                    let node_j: NodeIndex = NodeIndex::new(j);
                    let (_, count_i) = node_attributes.get(&node_i).unwrap();
                    let (_, count_j) = node_attributes.get(&node_j).unwrap();

                    // Add edges based on UMI count conditions
                    if *count_i >= 5 * *count_j - 1 {
                        graph.add_edge(node_i, node_j, edit_distance);
                    }

                    if *count_j >= 5 * *count_i - 1 {
                        graph.add_edge(node_j, node_i, edit_distance);
                    }
                }
            }
        }
    }

    (graph, node_attributes)
}

/// Find representative UMIs from a UMI graph.
///
/// This function examines each connected component in the given UMI graph and identifies
/// the representative UMI for each component. A representative UMI is defined as the UMI
/// within its connected component that has the highest count.
/// 
/// This function uses a breadth-first search (BFS) to traverse each connected component, keeping 
/// track of the UMI with the highest count in the current component. The UMI with the highest 
/// count is then added to the representative UMIs HashSet.
///
/// # Arguments
///
/// * `graph` - A reference to the UMI graph which is of type `Graph<Vec<u8>, i32>`. Each node
/// represents a unique UMI, and each edge represents a edit distance less than a specified
/// threshold between two UMIs.
///
/// * `node_attributes` - A reference to a HashMap mapping node indices to a tuple of UMIs and 
/// their counts. The key is a `NodeIndex`, and the value is a tuple where the first element
/// is a UMI (as a byte array), and the second element is its count.
///
/// # Returns
///
/// * A HashSet containing the representative UMIs for each connected component in the graph. 
/// Each UMI is represented as a byte array (`Vec<u8>`).
///
/// # Example
///
/// ```rust
/// // assume `graph` and `node_attributes` are predefined
/// let representative_umis = get_representative_umis_bfs(&graph, &node_attributes);
/// ```
pub fn get_representative_umis_bfs(graph: &Graph<Vec<u8>, i32>
                                , node_attributes: &HashMap<NodeIndex
                                , (Vec<u8>, i32)>) -> HashSet<Vec<u8>> {
    // Initialize an empty HashSet to store the representative UMIs.
    let mut repr_umis: HashSet<Vec<u8>> = HashSet::new();

    // Initialize a boolean vector to keep track of visited nodes.
    let mut visited = graph.visit_map();

    // Iterate over each node in the graph.
    for start_node in graph.node_indices() {
        // If the current node has not been visited, perform a BFS from it.
        if !visited.is_visited(&start_node) {
            let mut bfs = Bfs::new(graph, start_node);
            // Initialize variables to keep track of the UMI with the maximum count in the current component.
            let mut max_count:i32 = 0;
            let mut max_umi: Vec<u8> = Vec::new();

            while let Some(node) = bfs.next(graph) {
                // Mark the current node as visited.
                visited.visit(node);

                // Get the UMI and its count for the current node.
                let cur_umi: &(Vec<u8>, i32) = node_attributes.get(&node).unwrap();

                // If the current UMI count is greater than the max count seen so far, update max_count and max_umi.
                if cur_umi.1 > max_count {
                    max_count = cur_umi.1;
                    max_umi = cur_umi.0.clone();
                }
            }

            // Insert the UMI with the maximum count in the current component into the set of representative UMIs.
            repr_umis.insert(max_umi);
        }
    }

    // Return the set of representative UMis.
    repr_umis
}

/// Process a BAM file and error correct UMIs
///
/// This function reads a BAM file, extracts UMIs from read names, and stores their counts in a dictionary.
/// It then uses a directional graph algorithm to correct UMI errors, and returns a set of corrected UMIs.
///
/// # Arguments
///
/// * `input_bam_file` - Path to an input BAM file. The BAM file should contain only mapped 
///   RNA sequences and each read's UMI in the read name, appended using an underscore, e.g. 
///   "..._CGATGATCGATG".
///
/// # Returns
///
/// * A HashSet of corrected UMIs.
///
/// # Example
///
/// ```rust
/// let corrected_umis = find_true_umis("sample1.UMI.bam");
/// ```
pub fn find_true_umis(input_bam_file: &str, use_levenshtein: bool) -> Result<HashSet<Vec<u8>>, io::Error> {
    let mut umi_counts: HashMap<Vec<u8>, i32> = HashMap::new();
    let bam_reader = bam::BamReader::from_path(input_bam_file, 11).unwrap();

    for result in bam_reader {
        let record = result.unwrap();
        let qname = record.name();

        if let Some(captures) = READ_NAME_UMI_REGEX.captures(std::str::from_utf8(qname).unwrap()) {
            let umi = captures.get(1).unwrap().as_str().as_bytes().to_vec();
            *umi_counts.entry(umi).or_insert(0) += 1;
        }
    }

    if umi_counts.len() == 0 {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "No UMIs found"));
    }

    let (graph, node_attributes) = umi_graph(&umi_counts, 1, use_levenshtein);

    // Find the representative UMIs
    let corrected_umis = get_representative_umis_bfs(&graph, &node_attributes);

    Ok(corrected_umis)
}

/// Deduplicate a BAM file using a set of unique molecular identifiers (UMIs) and their mapping positions.
///
/// This function reads a BAM file and creates a deduplicated version of it based on a provided set of
/// unique molecular identifiers (UMIs) and their mapping positions. The UMI associated with each read in the
/// BAM file is expected to be appended to the read's name, prefixed with an underscore ('_').
///
/// A read is written to the output file if its UMI is in the set of provided UMIs and its combination of
/// UMI and mapping position hasn't been observed before. If a read's UMI is not in the set or the
/// UMI-mapping position combination has been observed before, the read is discarded.
///
/// # Arguments
///
/// * `input_bam_file` - A string representing the path to the input BAM file.
/// * `sample_name` - A string representing the sample name, which is used to name the output BAM file.
/// * `true_umis` - A HashSet of UMIs (`Vec<u8>`), which are considered to be the true UMIs.
///
/// # Output
///
/// This function doesn't return any value. Instead, it writes to an output BAM file named 
/// `{sample_name}.dedup.bam`, where `{sample_name}` is the provided sample name.
///
/// # Example
///
/// ```rust
/// // assume `true_umis` is predefined
/// let result = deduplicate_bam("input.bam", "sample1", true_umis);
/// // handle the error, if any
/// match result {
///     Ok(_) => println!("Deduplication completed successfully"),
///     Err(e) => eprintln!("An error occurred: {}", e),
/// }
/// ```
pub fn deduplicate_bam(input_bam_file: &str, sample_name: &String, true_umis: HashSet<Vec<u8>>) -> Result<(), Box<dyn Error>> {
    // Open the input BAM file with 4 threads for reading
    let bam_reader = bam::BamReader::from_path(input_bam_file, 4)?;

    // Capture the header of the input BAM file for use in the output BAM file
    let header = bam_reader.header().clone();

    // Create an output BAM file with the same header as the input BAM file
    let output_file_path = format!("{}.dedup.bam", sample_name);
    let mut bam_writer = bam::BamWriter::from_path(output_file_path, header)?;

    // Initialize a HashMap to store observed combinations of UMIs and their mapping positions
    let mut umis_and_positions_added: HashMap<(Vec<u8>, i32), bool> = HashMap::new();

    // Go through each read in the BAM file
    for result in bam_reader {
        let record = result?;

        // Extract the UMI from the read's name by splitting the name by underscores and taking the last element
        let qname = std::str::from_utf8(record.name())?;
        let umi: Vec<u8> = qname.split("_").last().unwrap().bytes().collect();

        // If the UMI exists and the combination of the UMI and the read's mapping position hasn't been observed before,
        // write the read to the output file and record the UMI-mapping position combination
        if true_umis.contains(&umi) && !umis_and_positions_added.contains_key(&(umi.clone(), record.ref_id())) {
            // Write the read to the output BAM file
            bam_writer.write(&record)?;

            // Record the observed combination of UMI and mapping position
            umis_and_positions_added.insert((umi, record.ref_id()), true);
        }
    }

    // Close the output BAM file
    bam_writer.finish()?;

    Ok(())
}