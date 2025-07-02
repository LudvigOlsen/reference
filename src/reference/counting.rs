use crate::{cli::BigCount, reference::kmer_codec::*};
use fxhash::FxHashMap;
use smallvec::SmallVec;

/// Count k-mers for every window on one chromosome
///
/// * `encs`       – slice of Enc {k, codes, none, n}
/// * `windows`    – (start, end, _original_idx) for every window
/// * `chrom_len`  – chromosome length (used to cap end)
///
/// Returns `Vec<FxHashMap<Kmer, BigCount>>` in the same order as `windows`.
pub fn count_kmers_by_window(
    counts_by_window: &mut Vec<FxHashMap<Kmer, BigCount>>,
    encs: &SmallVec<[Enc; 8]>,
    windows: &[(u64, u64, u64)],
    chrom_len: u64,
) {
    for (win_idx, &(win_start, mut win_end, _)) in windows.iter().enumerate() {
        let counts = &mut counts_by_window[win_idx.clone()];
        win_end = win_end.min(chrom_len as u64);

        for ref_pos in win_start..win_end {
            let remaining = win_end - ref_pos; // bp left in the window
            for enc in encs {
                let k = enc.k;
                if remaining < enc.k as u64 {
                    // k-mer would over-run
                    continue;
                }
                let code = enc.codes.get(ref_pos as usize);

                if code == enc.none || code == enc.n {
                    continue;
                }

                *counts.entry(Kmer { k, code }).or_insert(0) += 1;
            }
        }
    }
}

/// Container for storing k, codes, and sentinels
pub struct Enc<'a> {
    pub k: u8,
    pub codes: &'a KmerCodes,
    pub none: u64,
    pub n: u64,
}
