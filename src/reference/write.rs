use crate::cli::BigCount;
use crate::reference::kmer_codec::{DecodedCounts, KmerSpec};
use crate::reference::process_counts::motif_order;
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
        write_category(&ref_bins, &tag, output_dir)?;
    }

    Ok(())
}

/// Write <prefix>_counts.npy and <prefix>_motifs.txt
fn write_category(
    bins: &[FxHashMap<String, BigCount>],
    prefix: &str,
    out_dir: &Path,
) -> anyhow::Result<()> {
    if bins.is_empty() {
        return Ok(()); // nothing to write
    }

    // Consistent column order = sorted keys of the first bin
    let motifs = motif_order(&bins);

    // Map motif → column index
    let col_of: FxHashMap<_, _> = motifs.iter().enumerate().map(|(i, m)| (m, i)).collect();

    let n = bins.len();
    let m = motifs.len();
    let mut mat = Array2::<BigCount>::zeros((n, m));

    for (r, hm) in bins.iter().enumerate() {
        for (motif, &cnt) in hm {
            if let Some(&c) = col_of.get(motif) {
                mat[(r, c)] = cnt;
            }
        }
    }

    // Write matrix
    write_npy(out_dir.join(format!("{}_counts.npy", prefix)), &mat)?;
    // Write motif list
    let mut f = File::create(out_dir.join(format!("{}_motifs.txt", prefix)))?;
    for mo in motifs {
        writeln!(f, "{}", mo)?;
    }
    Ok(())
}
