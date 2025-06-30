use clap::{value_parser, Args};
use std::path::PathBuf;

#[derive(Debug, Args)]
pub struct IOArgs {
    /// Indexed, coordinate-sorted BAM input file [path]
    #[clap(
        short = 'i',
        long,
        value_parser,
        required = true,
        help_heading = "Core"
    )]
    pub bam: PathBuf,

    /// Output directory for results [path]
    #[clap(
        short = 'o',
        long,
        value_parser,
        required = true,
        help_heading = "Core"
    )]
    pub output_dir: PathBuf,

    /// Number of threads to use (increases RAM usage) [integer]
    #[clap(short = 't', long, default_value = "1", help_heading = "Core")]
    pub n_threads: usize,
}

#[derive(Debug, Args)]
pub struct ReadFilteringArgs {
    /// Minimum mapping quality to include [integer]
    #[clap(long, alias = "mq", default_value = "60", value_parser = value_parser!(u8).range(0..), help_heading="Filtering")]
    pub min_mapq: u8,

    /// Minimum read sequence length (default: 50) [integer]
    #[clap(long, default_value = "50", value_parser = value_parser!(u8).range(0..), help_heading="Filtering")]
    pub min_seq_len: u8,

    /// Maximum fragment length (insert size) to include (default: 300) [integer]
    ///
    /// Smaller max. length --> faster and lower RAM / storage
    #[clap(long, default_value = "300", value_parser = value_parser!(u16).range(1..), help_heading="Filtering")]
    pub max_fragment_length: u16, // NOTE: Also necessary for streaming logic!

    /// Minimum number of recorded mismatches in a read [integer]
    #[clap(long, default_value = "0", value_parser = value_parser!(u16).range(0..), help_heading="Filtering")]
    pub min_nm: u16,

    /// Maximum number of recorded mismatches in a read [integer]
    #[clap(long, default_value = "5", value_parser = value_parser!(u16).range(0..), help_heading="Filtering")]
    pub max_nm: u16,
}

#[derive(Debug, Args)]
pub struct GCArgs {
    /// Enable GC binning (counts per GC-content bin) [flag]
    #[clap(short = 'g', long, alias = "bin-by-gc", requires = "ref_2bit")]
    pub bin_by_gc: bool,

    /// GC bin size (% for GC binning) [integer]
    #[clap(long, default_value = "3", requires = "bin_by_gc", value_parser = value_parser!(u8).range(1..100))]
    pub gc_bin_size_pct: u8,

    /// Minimum GC % to consider [integer]
    #[clap(long, default_value = "0", requires = "bin_by_gc", value_parser = value_parser!(u8).range(0..100))]
    pub gc_min_pct: u8,

    /// Maximum GC % to consider [integer]
    #[clap(long, default_value = "100", requires = "bin_by_gc", value_parser = value_parser!(u8).range(0..101))]
    pub gc_max_pct: u8,
}
