use anyhow::{anyhow, Context, Result};

use rust_htslib::bam::{IndexedReader, Read};
use std::path::Path;
use twobit::TwoBitFile;
// BAM

pub fn create_chromosome_reader(bam_path: &Path, chr: &str) -> Result<(IndexedReader, u32, u64)> {
    let reader = IndexedReader::from_path(bam_path).context(format!("opening BAM for {}", chr))?;
    let header = reader.header().to_owned();
    let tid = header
        .tid(chr.as_bytes())
        .ok_or_else(|| anyhow!("{} not in BAM", chr))?;
    let chrom_len = header
        .target_len(tid)
        .ok_or_else(|| anyhow!("No length for {}", chr))? as u64;
    Ok((reader, tid, chrom_len))
}

// Reference 2bit file

pub fn read_seq(path: &Path, chr: &str) -> anyhow::Result<Vec<u8>> {
    // open once
    let mut tb = TwoBitFile::open(path).context("opening 2bit")?;
    // Get reference sequence once
    let seq = tb
        .read_sequence(chr, ..)
        .context(format!("extracting reference seq for {}", chr))?;
    Ok(seq.as_bytes().to_vec())
}
