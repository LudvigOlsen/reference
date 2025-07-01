use crate::cli::BigCount;
use crate::reference::kmer_codec::{DecodedCounts, KmerSpec};
use fxhash::FxHashMap;
use ndarray::Array2;
use ndarray_npy::write_npy;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::Path;

/// Write one `.npy` matrix and a companion `*_motifs.txt` file for every
/// k present in `prepared_windows`.
///
/// * `prepared_windows` – windows of decoded counts.
/// * `kmer_specs`       – validated specs: the keys determine which k values
///                        will be written, and in which order.
/// * `output_dir`       – target directory.
///
/// * For reference windows the files are named  `k<k>_counts.npy`, e.g.
///   `k3_counts.npy`.  
///
/// The matrix dimensions are **windows × motifs** with the same column order
/// used across all windows of that k-mer size.
pub fn write_decoded_counts_matrix(
    prepared_windows: &[DecodedCounts],
    kmer_specs: &HashMap<u8, KmerSpec>,
    motifs_by_k: &HashMap<u8, Vec<String>>,
    output_dir: &Path,
) -> anyhow::Result<()> {
    let n_win = prepared_windows.len();

    for &k in kmer_specs.keys() {
        // Collect reference bins for this k
        let mut ref_bins: Vec<FxHashMap<String, BigCount>> = vec![FxHashMap::default(); n_win];
        for (idx, win) in prepared_windows.iter().enumerate() {
            if let Some(bin) = win.counts.get(&k) {
                ref_bins[idx] = bin.clone();
            }
        }
        let tag = format!("k{}", k);
        write_category(&mut ref_bins, &tag, &motifs_by_k[&k], output_dir)?;
    }

    Ok(())
}

/// Write <prefix>_counts.npy and <prefix>_motifs.txt
///
/// * `motifs`  - The motifs to include for all bins in the order you want it saved in.
fn write_category(
    bins: &[FxHashMap<String, BigCount>],
    prefix: &str,
    motifs: &[String],
    out_dir: &Path,
) -> anyhow::Result<()> {
    if bins.is_empty() {
        return Ok(()); // nothing to write
    }

    // Output matrix
    let n_rows = bins.len();
    let n_cols = motifs.len();
    let mut mat = Array2::<BigCount>::zeros((n_rows, n_cols));

    // Pre-compute motif → column index once
    let col_of: FxHashMap<_, _> = motifs.iter().enumerate().map(|(c, m)| (m, c)).collect();

    for (row, hm) in bins.iter().enumerate() {
        for (motif, &cnt) in hm {
            if let Some(&col) = col_of.get(motif) {
                mat[(row, col)] = cnt; // Counts overwrite the zero
            }
        }
    }

    // Persist outputs
    write_npy(out_dir.join(format!("{prefix}_counts.npy")), &mat)?;

    let mut txt = File::create(out_dir.join(format!("{prefix}_motifs.txt")))?;
    for m in motifs {
        writeln!(txt, "{m}")?;
    }

    Ok(())
}
