use anyhow::Context;

use std::path::Path;
use twobit::TwoBitFile;
// BAM

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
