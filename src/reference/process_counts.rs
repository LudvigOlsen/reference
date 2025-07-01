use fxhash::FxHashMap;

use crate::cli::BigCount;

use crate::reference::kmer_codec::{DecodedCounts, KmerSpec};
use std::collections::{HashMap, HashSet};

fn prepare_kmer_category(
    windows: &[DecodedCounts],
    kmer_specs: &HashMap<u8, KmerSpec>,
    k: usize,
    canonical: bool,
    ensure_all: bool,
) -> (Vec<FxHashMap<String, BigCount>>, Vec<String>) {
    // Extract the raw maps
    let raw_bins = extract_bins(windows, k, canonical);

    // Build the (canonical) motif list once, if requested.
    let base_motifs: Vec<String> = if ensure_all {
        all_motifs(k, kmer_specs)
    } else {
        Vec::new()
    };

    // Build the (canonical) motif list *once* so we know what to pad with
    let mut motifs = collect_motifs(&raw_bins, base_motifs, canonical, ensure_all);
    motifs.sort_unstable();

    (raw_bins, motifs)
}

/// Prepare decoded counts for all kmer sizes in all windows.
///
/// Extracts motifs per kmer spec to allow future padding.
/// For kmers of size 1..6, this includes all possible motifs.
/// For larger kmer sizes, only the seen motifs is included as the number otherwise explodes.
///
/// * `windows`        – slice of per-window raw counts
/// * `canonical`      – canonical reverse complements when true
/// * `kmer_specs`     – validated specs for every k we want to keep
///
pub fn prepare_decoded_counts(
    windows: &[DecodedCounts],
    canonical: bool,
    kmer_specs: &HashMap<u8, KmerSpec>,
) -> (Vec<DecodedCounts>, HashMap<u8, Vec<String>>) {
    let n_windows = windows.len();

    // Initialise one empty DecodedCounts per window
    let mut out = vec![
        DecodedCounts {
            counts: HashMap::new()
        };
        n_windows
    ];

    let mut motifs_by_k: HashMap<u8, Vec<String>> = HashMap::new();

    // Loop over every k we validated
    for (&k, _) in kmer_specs {
        // Reference (match) bins for this k
        let (count_bins, motifs) =
            prepare_kmer_category(windows, kmer_specs, k as usize, canonical, k <= 6);

        // Insert into the corresponding window
        for i in 0..n_windows {
            out[i].counts.insert(k, count_bins[i].clone());
        }
        motifs_by_k.insert(k, motifs);
    }

    (out, motifs_by_k)
}

/// Collect per-window bins for the requested motif type and (optionally)
/// canonical them into strand-agnostic form.
///
/// * `windows` – slice of `DecodedCounts` (“one window” each).
/// * `k` – kmer-size to pull out of every `DecodedCounts`.
/// * `canonical` – if `true`, run the appropriate collapse_*_map helper.
///
/// Returns a fresh `Vec<FxHashMap<String, BigCount>>` – one map per window.
fn extract_bins(
    windows: &[DecodedCounts],
    k: usize, // pattern only; field values are ignored
    canonical: bool,
) -> Vec<FxHashMap<String, BigCount>> {
    windows
        .iter()
        .map(|dc| {
            // 1. Pick the raw map for this window
            let raw: FxHashMap<String, BigCount> =
                dc.counts.get(&(k as u8)).cloned().unwrap_or_default();

            // 2. Collapse if requested, otherwise return the raw map
            if canonical {
                collapse_map(&raw)
            } else {
                raw
            }
        })
        .collect()
}

/// Collect motifs for a category, optionally ensuring the full universe and filtering 'N'
fn collect_motifs(
    windows: &[FxHashMap<String, BigCount>],
    base_motifs: Vec<String>,
    canonical: bool,
    ensure_all: bool,
) -> Vec<String> {
    // Universe of motifs to keep
    let set: HashSet<String> = if ensure_all {
        base_motifs.into_iter().collect()
    } else {
        windows.iter().flat_map(|m| m.keys().cloned()).collect()
    };

    // Strand-collapse if requested
    let collapsed_set = if canonical { collapse_set(&set) } else { set };

    // Convert to sorted Vec
    let mut v: Vec<String> = collapsed_set.into_iter().collect();
    v.sort_unstable();
    v
}

/// Use the first window’s keys, sort them, and return the order.
/// Panics only if `bins` is empty.
pub fn motif_order(bins: &[FxHashMap<String, impl Copy>]) -> Vec<String> {
    assert!(
        !bins.is_empty(),
        "motif_order: received an empty slice of bins"
    );
    let mut motifs: Vec<String> = bins[0].keys().cloned().collect();
    motifs.sort_unstable();
    motifs
}

/// Return all possible reference motifs (4ᵏ) for a given k.
///
/// No motifs with 'N' are returned.
pub fn all_motifs(k: usize, specs: &HashMap<u8, KmerSpec>) -> Vec<String> {
    let spec = &specs[&(k as u8)];
    let max_code = 5u64.pow(k as u32) - 1; // no-N space
    (0..=max_code)
        .map(|c| spec.decode_kmer(c))
        .filter(|m| !m.contains('N'))
        .collect()
}

// Collapsing of motifs

/// Complement of a single nucleotide base
#[inline]
fn comp(b: char) -> char {
    match b {
        'A' | 'a' => 'T',
        'T' | 't' => 'A',
        'C' | 'c' => 'G',
        'G' | 'g' => 'C',
        'N' | 'n' => 'N',
        _ => b,
    }
}

/// Reverse-complement of a plain sequence, e.g. "AC" → "GT"
fn revcomp(seq: &str) -> String {
    seq.chars().rev().map(comp).collect()
}

/// Collapse a map of reference k-mer counts into canonical keys, summing counts
pub fn collapse_map(map: &FxHashMap<String, u64>) -> FxHashMap<String, u64> {
    let mut out: FxHashMap<String, u64> = FxHashMap::default();
    for (kmer, &count) in map {
        let canon = canonical(kmer.to_owned());
        *out.entry(canon).or_default() += count;
    }
    out
}

/// Collapse a set of motifs into canonical form
pub fn collapse_set(set: &HashSet<String>) -> HashSet<String> {
    set.iter().map(|kmer| canonical(kmer.to_owned())).collect()
}

/// Return the canonical form of `kmer`: the lexicographically smaller
/// of the k-mer and its reverse complement.
#[inline]
fn canonical(kmer: String) -> String {
    let rc = revcomp(&kmer);
    if kmer <= rc {
        kmer
    } else {
        rc
    }
}
