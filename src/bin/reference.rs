use anyhow::{Context, Result};
use clap::ArgAction;
use clap::{value_parser, ArgGroup, Parser};
use fxhash::FxHashMap;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use reference::cli::io::read_seq;
use reference::cli::BigCount;
use reference::reference::bed::load_windows;
use reference::reference::blacklist::*;
use reference::reference::kmer_codec::*;
use reference::reference::process_counts::prepare_decoded_counts;
use reference::reference::write::write_decoded_counts_matrix;
use smallvec::SmallVec;
use std::mem::drop;
use std::{
    collections::HashMap,
    fs::{create_dir_all, File},
    io::{BufWriter, Write},
    path::PathBuf,
    sync::Arc,
    time::Instant,
};

/// Command-line options for fragment length extraction tool
#[derive(Parser)]
#[command(
    name = "reference",
    about = "Count reference kmers in genomic windows",
    long_about = "Count reference kmers in genomic windows.
    

EXAMPLES:
    // Using defaults
    $ reference --ref-2bit <path/to/hg38.2bit> --output-dir <path/to/output_directory/> --kmer-sizes 3 --n-threads <N> --global -b <path/to/blacklist_1.bed> -b <path/to/blacklist_2.bed>
    ",
    author = "Ludvig Renbo Olsen",
    version = "0.0.1"
)]
#[clap(group = ArgGroup::new("windows").required(true).args(&["by_size", "by_bed", "global"]).multiple(false))]
#[clap(group = ArgGroup::new("chrom_select").args(&["chromosomes", "chromosomes_file"]).multiple(false))]
struct Cli {
    /// 2bit reference file [path]
    /// E.g., "hg38.2bit"
    #[clap(
        short = 'r',
        long,
        value_parser,
        required = true,
        help_heading = "Core"
    )]
    pub ref_2bit: PathBuf,

    /// Output directory for results [path]
    #[clap(
        short = 'o',
        long,
        value_parser,
        required = true,
        help_heading = "Core"
    )]
    pub output_dir: PathBuf,

    /// List of K-mer sizes [integer].
    #[clap(short = 'k', long, num_args = 1.., value_parser = value_parser!(u8).range(1..28), value_delimiter = ',', required=true, help_heading="Core")]
    pub kmer_sizes: Vec<u8>,

    /// Number of threads to use (increases RAM usage) [integer]
    #[clap(short = 't', long, default_value = "1", help_heading = "Core")]
    pub n_threads: usize,

    /// Use a fixed window size [integer]
    #[clap(
        long = "by-size",
        alias = "by",
        value_parser,
        group = "windows",
        help_heading = "Windows (select one)"
    )]
    pub by_size: Option<usize>,

    /// Use a BED file of windows [path]
    #[clap(
        long = "by-bed",
        value_parser,
        group = "windows",
        help_heading = "Windows (select one)"
    )]
    pub by_bed: Option<PathBuf>,

    /// Use a single genome-wide window [flag]
    #[clap(
        long = "global",
        group = "windows",
        help_heading = "Windows (select one)"
    )]
    pub global: bool,

    /// Names of chromosomes to process (comma-separated or repeated). E.g. 'chr1,chr2,chr3'.
    ///
    /// When no chromosomes are specified, it defaults to chr1..chr22.
    #[clap(long, num_args = 1.., value_parser, value_delimiter = ',', group = "chrom_select", help_heading="Chromosome Selection (select max. one)")]
    pub chromosomes: Option<Vec<String>>,

    /// File with chromosome names to process (one per line).
    #[arg(
        long,
        value_parser,
        group = "chrom_select",
        help_heading = "Chromosome Selection (select max. one)"
    )]
    pub chromosomes_file: Option<PathBuf>,

    /// Optional BED files of blacklisted regions [path]
    #[clap(short = 'b', long, value_parser, num_args = 1.., action = ArgAction::Append, help_heading="Filtering")]
    pub blacklist: Option<Vec<PathBuf>>,

    /// Minimum size of blacklist intervals to load (bp) [integer]
    #[clap(
        long,
        alias = "bl-min-size",
        default_value = "1",
        help_heading = "Filtering"
    )]
    pub blacklist_min_size: u64,

    /// Collapse each kmer with its reverse-complement. [flag]
    ///
    /// The lexicographically lowest kmer is used.
    #[clap(short = 'c', long, help_heading = "Core")]
    canonical: bool,

    /// Save counts as sparse-array. [flag]
    ///
    /// For large kmer-sizes, we cannot save dense arrays with all motifs
    /// unless we have a LOT of RAM and storage space. Enable this
    /// flag to save as a COO sparse array that can be opened in
    /// python via `scipy.sparse.load_npz()`.
    #[clap(long, help_heading = "Core")]
    pub save_sparse: bool,
}

impl Cli {
    /// Returns the final chromosome list, in priority order:
    /// 1) from `--chromosomes-file`
    /// 2) from `--chromosomes`
    /// 3) default `chr1`..`chr22`
    pub fn resolve_chromosomes(&self) -> anyhow::Result<Vec<String>> {
        if let Some(file) = &self.chromosomes_file {
            let text: String = std::fs::read_to_string(file)
                .context(format!("reading chromosome file {:?}", file))?;
            let list: Vec<String> = text
                .lines()
                .map(str::trim)
                .filter(|l| !l.is_empty() && !l.starts_with('#'))
                .map(String::from)
                .collect();
            Ok(list)
        } else if let Some(chrs) = &self.chromosomes {
            Ok(chrs.clone())
        } else {
            Ok((1..=22).map(|i| format!("chr{}", i)).collect())
        }
    }
}

fn main() {
    // Catch and handle errors
    // Ensures that tempfile has time to remove the tmp dir
    if let Err(e) = run() {
        eprintln!("{:?}", e);
        std::process::exit(1);
    }
    std::process::exit(0);
}

fn run() -> Result<()> {
    let start_time = Instant::now();
    let opt = Cli::parse();
    let chromosomes = opt.resolve_chromosomes()?;
    let pb = Arc::new(ProgressBar::new(chromosomes.len() as u64));
    pb.set_style(
        ProgressStyle::default_bar()
            .template("       {bar:40} {pos}/{len} [{elapsed_precise}] {msg}")
            .unwrap(),
    );

    // Create output directory
    create_dir_all(&opt.output_dir).context("Cannot create output_dir")?;

    // Load blacklist intervals if provided
    let blacklist_map = if let Some(beds) = &opt.blacklist {
        println!("Start: Loading blacklists");
        load_blacklists(beds, opt.blacklist_min_size, &chromosomes)?
    } else {
        HashMap::new()
    };

    let windows_map = if let Some(bed) = &opt.by_bed {
        println!("Start: Loading window coordinates");
        Some(load_windows(bed, &chromosomes)?)
    } else {
        None
    };

    let kmer_specs: HashMap<u8, KmerSpec> = build_kmer_specs(&opt.kmer_sizes)?;

    // Configure global thread‐pool size
    rayon::ThreadPoolBuilder::new()
        .num_threads(opt.n_threads as usize)
        .build_global()
        .context("building Rayon thread pool")?;

    // Prepare per-bin counts and metadata
    let mut all_bins = Vec::new();
    let mut bin_info = Vec::new();

    // Main loop: process each autosome
    println!("Start: Counting per chromosome");

    pb.set_position(0);

    let results: Vec<(
        Vec<FxHashMap<Kmer, BigCount>>,
        Vec<(String, u64, u64, u64, f64)>,
    )> = chromosomes
        .par_iter()
        .map(|chr| -> Result<(_, _)> {
            let out = process_chrom(
                &chr,
                &opt,
                &kmer_specs,
                windows_map
                    .as_ref()
                    .and_then(|m| m.get(chr).map(|v| v.as_slice())),
                //gc_bins,
                blacklist_map.get(chr).map(|v| v.as_slice()).unwrap_or(&[]),
            )?;
            pb.inc(1);
            Ok(out)
        })
        .collect::<Result<_>>()?; // short-circuits on the first Err

    pb.finish_with_message("| Finished counting");

    println!("Start: Processing counts");

    // Collect results (in chromosome order) back into the global vectors
    for (counts_by_bin, bin_vec) in results {
        let counts_decoded: Vec<DecodedCounts> = counts_by_bin
            .iter()
            .map(|c| split_and_decode_counts(c, &kmer_specs))
            .collect();
        all_bins.extend(counts_decoded);
        if !opt.global {
            bin_info.extend(bin_vec);
        }
    }

    // Convert to single hashmap for global
    // Keep wrapped in vector to simplify writer
    let all_bins = if opt.global {
        vec![merge_decoded_counts(all_bins)]
    } else {
        all_bins
    };

    // Prepare to get correct motifs (collapsed, N-filtered, etc.)
    let (mut prepared_counts, motifs_by_k) =
        prepare_decoded_counts(&all_bins, opt.canonical, &kmer_specs);

    if opt.by_bed.is_some() {
        println!("Start: Reordering counts by original window index in bed file");

        // Zip into a single Vec
        let mut paired: Vec<_> = bin_info
            .into_iter()
            .zip(prepared_counts.into_iter())
            .collect(); // (BinInfo, DecodedCounts)

        // Sort primarily by original window index
        paired.sort_unstable_by_key(|(info, _)| info.3);

        // Unzip back out if you need separate Vecs again
        (bin_info, prepared_counts) = paired.into_iter().unzip();
    }

    println!("Start: Writing counts to disk");
    write_decoded_counts_matrix(
        &prepared_counts,
        &kmer_specs,
        &motifs_by_k,
        &opt.output_dir,
        opt.save_sparse,
    )?;

    // Write bins BED file
    if !opt.global {
        println!("Start: Writing window coordinates to disk");
        let mut bed_writer = BufWriter::new(
            File::create(&opt.output_dir.join("bins.bed")).context("Create bed fail")?,
        );
        for (chr, start, end, _, overlap_perc) in &bin_info {
            writeln!(bed_writer, "{}\t{}\t{}\t{}", chr, start, end, overlap_perc)
                .context("Write bed line fail")?;
        }
    }

    // Print summary statistics and execution time
    let elapsed = start_time.elapsed();
    println!("Elapsed time: {:.2?}", elapsed);
    Ok(())
}

/* ---------- main routine -------------------------------------------- */

/// * windows  -  Optional slice of tuples with (start, end, original_idx)
fn process_chrom(
    chr: &str,
    opt: &Cli,
    kmer_specs: &HashMap<u8, KmerSpec>,
    windows: Option<&[(u64, u64, u64)]>,
    // gc_bins: usize,
    blacklist_intervals: &[(u64, u64)],
) -> anyhow::Result<(
    Vec<FxHashMap<Kmer, BigCount>>,
    Vec<(String, u64, u64, u64, f64)>,
)> {
    let mut seq_bytes = read_seq(&opt.ref_2bit, chr)?;
    apply_blacklist_mask_to_seq(&mut seq_bytes, &blacklist_intervals);
    let chrom_len = seq_bytes.len() as usize;
    let positional_codes_by_k: HashMap<u8, KmerCodes> = build_codes_per_k(&seq_bytes, kmer_specs);

    // Delete seq_bytes from memory
    drop(seq_bytes);

    // Calculate window coordinates for all windowing options
    let windows: Vec<(u64, u64, u64)> = if let Some(sz) = opt.by_size {
        // by-size
        let num_windows = ((chrom_len + sz - 1) / sz) as usize;
        (0..num_windows)
            .map(|s| ((s * sz) as u64, (sz + s * sz) as u64, s as u64))
            .collect()
    } else if opt.by_bed.is_some() {
        // by-bed
        windows.unwrap().to_owned()
    } else {
        // global
        vec![(0, chrom_len as u64, 0u64)]
    };

    let num_windows = windows.len();

    let mut counts_by_window = vec![FxHashMap::<Kmer, BigCount>::default(); num_windows];

    let mut encs: SmallVec<[Enc; 8]> = SmallVec::new();
    for (&k, spec) in kmer_specs {
        encs.push(Enc {
            k,
            codes: &positional_codes_by_k[&k],
            none: spec.sentinel_none(),
            n: spec.sentinel_n(),
        });
    }

    for (win_idx, &(win_start, mut win_end, _)) in windows.iter().enumerate() {
        let counts = &mut counts_by_window[win_idx.clone()];
        win_end = win_end.min(chrom_len as u64);

        for ref_pos in win_start..win_end {
            for enc in &encs {
                let k = enc.k;
                let code = enc.codes.get(ref_pos as usize);

                if code == enc.none || code == enc.n {
                    continue;
                }

                *counts.entry(Kmer { k, code }).or_insert(0) += 1;
            }
        }
    }

    let bin_info = {
        // build bin_info from the exact BED windows
        let mut bl_ptr = 0;
        let mut bin_info = Vec::with_capacity(num_windows);
        for (_b, (win_start, mut win_end, original_win_idx)) in windows.iter().cloned().enumerate()
        {
            win_end = win_end.min(chrom_len as u64);
            let overlap_perc =
                compute_blacklist_overlap(blacklist_intervals, win_start, win_end, &mut bl_ptr);
            bin_info.push((
                chr.to_string(),
                win_start,
                win_end,
                original_win_idx,
                overlap_perc,
            )); // total,
        }
        bin_info
    };

    Ok((counts_by_window, bin_info))
}

struct Enc<'a> {
    k: u8,
    codes: &'a KmerCodes,
    none: u64,
    n: u64,
}
