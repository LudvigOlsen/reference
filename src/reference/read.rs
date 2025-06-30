use crate::cli::opts::ReadFilteringArgs;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Record;
use std::collections::HashSet;

pub fn filter_read(rec: &Record, opt: &ReadFilteringArgs) -> Option<()> {
    // Only get the reads of mapped and paired fragments
    // that are not secondary, supplementary, duplicate, has
    // failed QC, or has too low mapping quality
    if rec.mapq() < opt.min_mapq
        || rec.is_unmapped()
        || rec.is_secondary()
        || rec.is_supplementary()
        || rec.is_duplicate()
        || rec.is_mate_unmapped()
        || !rec.is_proper_pair()
        || rec.is_quality_check_failed()
        || rec.insert_size() == 0
        || rec.seq_len() < opt.min_seq_len as usize
        || rec.insert_size().abs() < opt.min_seq_len as i64
        || rec.insert_size().abs() > opt.max_fragment_length as i64
    // Consider this
    {
        return None;
    }

    // TODO: Move to fragment-level and check relative to overlap positions
    // as we don't need to skip if clipping or indels are outside of the overlap!
    for entry in rec.cigar().iter() {
        match entry {
            Cigar::Ins(_)
            | Cigar::Del(_)
            | Cigar::RefSkip(_)
            | Cigar::SoftClip(_)
            | Cigar::HardClip(_) => return None,
            _ => (),
        }
    }

    Some(())
}

/// Extract 'NM' aux tag from read as u16
pub fn read_nm_tag(rec: &Record) -> Option<u16> {
    // extract NM tag as u16, returning None on missing/out-of-range
    let nm: u16 = match rec.aux(b"NM") {
        // I believe it's always i32 but not sure
        Ok(Aux::I32(n)) if (0..=u16::MAX as i32).contains(&n) => n as u16,
        Ok(Aux::U8(n)) => n as u16,
        Ok(Aux::U16(n)) => n,
        Ok(Aux::I8(n)) if n >= 0 => n as u16,
        Ok(Aux::I16(n)) if n >= 0 => n as u16,
        Ok(Aux::U32(n)) if n <= u16::MAX as u32 => n as u16,
        _ => return None,
    };
    Some(nm)
}

/// Read the MD auxiliary tag from a BAM record
///
/// The MD tag (type Z) encodes mismatches between the read and reference
/// using a compact format from the SAM specification:
///   MD:Z:<string>
/// where the string consists of:
/// - numbers: counts of matching bases
/// - uppercase letters: reference bases at mismatch positions
/// - '^' followed by letters: deleted reference sequence
///
/// This function returns the raw MD string (without the "MD:Z:" prefix)
/// or None if the tag is missing or not a string
pub fn read_md_tag(rec: &Record) -> Option<String> {
    match rec.aux(b"MD") {
        Ok(Aux::String(s)) => Some(s.to_owned()),
        _ => None,
    }
}

/// Parse an MD tag string into mismatch positions
///
/// Given an MD tag (e.g. "10A5^AC2" or "100"), returns two parallel vectors:
/// - starts: zero-based reference positions where each mismatch run begins
/// - ends: position just after each mismatch run ends
///
/// The MD tag string has a structure like:
///   <match><ref_base><match>^<del_bases><match><ref_base>...
/// examples:
/// - "10A5" → 10 matches, A mismatch, 5 matches
/// - "8^AC2T3" → 8 matches, deletion of AC, 2 matches, T mismatch, 3 matches
///
/// Interpretation rules:
/// 1. A number advances the reference cursor by that many matches
/// 2. A letter A/T/C/G indicates a single-base mismatch (cursor +1)
/// 3. '^' followed by letters indicates deletions (subsequent letters are skipped)
/// 4. Mismatch runs are consecutive mismatches without intervening matches
/// 5. At each new mismatch run, record the start position; on next number >0 or end, record the end position
///
/// Returns empty vectors if no mismatches occur (purely numeric tag)
pub fn parse_md_tag(md_tag: &str, offset: u32) -> (Vec<u32>, Vec<u32>) {
    // These will hold the start and end positions of each mismatch run
    // Relative to the read
    let mut starts = Vec::new();
    let mut ends = Vec::new();

    // `pos` is our “cursor” along the reference sequence (0-based)
    let mut pos = 0u32;

    // `in_run` tells us whether we’re currently inside a stretch of mismatches
    let mut in_run = false;

    // Get a byte‐slice of the string, so we can inspect one ASCII code at a time
    let bytes = md_tag.as_bytes();
    let mut i = 0;

    // Walk through each byte in the MD tag
    while i < bytes.len() {
        let b = bytes[i];

        // ───────────────────────────────────────────────
        // 1) Digit?  (a run of matching bases)
        // ───────────────────────────────────────────────
        if b.is_ascii_digit() {
            // We’ll build up the full number (e.g. “123” → 123) in `num`
            let mut num: u32 = 0;

            // As long as the next byte is also a digit, keep consuming
            while i < bytes.len() && bytes[i].is_ascii_digit() {
                // (bytes[i] - b'0') turns the ASCII code for '0'..'9'
                // into the numeric 0..9
                //
                //   ASCII '0' = 48, '1' = 49, …, '9' = 57
                //   So '3' (51) minus '0' (48) gives 3
                let digit = (bytes[i] - b'0') as u32;

                // Accumulate into num: shift previous value *10, then add digit
                // e.g. reading "3","4","2" gives
                //   num = 0*10 + 3 = 3
                //   num = 3*10 + 4 = 34
                //   num = 34*10 + 2 = 342
                num = num * 10 + digit;

                i += 1;
            }

            // If we were inside a mismatch run, and now see a positive match-count,
            // that means the run just ended at the old `pos`
            if in_run && num > 0 {
                ends.push(pos);
                in_run = false;
            }

            // Advance our reference cursor by `num` matched bases
            pos += num;

        // ───────────────────────────────────────────────
        // 2) Caret '^'?  (deletion from the reference)
        // ───────────────────────────────────────────────
        } else if b == b'^' {
            // Skip deletion marker and deleted bases
            i += 1;
            // Count and skip each deleted reference base
            let mut del_len = 0;
            while i < bytes.len() && (bytes[i] as char).is_ascii_uppercase() {
                del_len += 1;
                i += 1;
            }
            // Deletions consume reference positions, so advance cursor
            pos += del_len;

        // ───────────────────────────────────────────────
        // 3) Uppercase letter?  (single‐base mismatch)
        // ───────────────────────────────────────────────
        } else if (b as char).is_ascii_uppercase() {
            // If we’re not already in a mismatch run, start one here:
            if !in_run {
                starts.push(pos);
                in_run = true;
            }
            // Consumed one mismatch base → advance cursor by 1
            pos += 1;
            i += 1;

        // ───────────────────────────────────────────────
        // 4) Anything else… just skip it
        // ───────────────────────────────────────────────
        } else {
            i += 1;
        }
    }

    // If we ended in the middle of a mismatch run, close it at the end of the tag
    if in_run {
        ends.push(pos);
    }

    if offset != 0 {
        let starts: Vec<u32> = starts
            .into_iter()
            .map(|s| offset.saturating_add(s))
            .collect();
        let ends: Vec<u32> = ends.into_iter().map(|e| offset.saturating_add(e)).collect();
        return (starts, ends);
    }

    (starts, ends)
}

/// Intersect mismatch runs from two reads sorted by start position
///
/// Given two slices of start and end positions for each read,
/// return a sorted Vec of (start, end) pairs present in both reads
pub fn intersect_mismatch_runs(
    starts1: &[u32],
    ends1: &[u32],
    starts2: &[u32],
    ends2: &[u32],
) -> Vec<(u32, u32)> {
    assert_eq!(
        starts1.len(),
        ends1.len(),
        "Starts1 and Ends1 must have the same length"
    );
    assert_eq!(
        starts2.len(),
        ends2.len(),
        "Starts2 and Ends2 must have the same length"
    );

    // Pair up the first read's start and end positions
    let runs1: Vec<(u32, u32)> = starts1.iter().cloned().zip(ends1.iter().cloned()).collect();

    // Build a set of runs from the second read for fast lookup
    let set2: HashSet<(u32, u32)> = starts2.iter().cloned().zip(ends2.iter().cloned()).collect();

    // Filter for runs present in both reads
    let mut shared: Vec<(u32, u32)> = runs1.into_iter().filter(|run| set2.contains(run)).collect();

    // Sort the shared runs by their start position
    shared.sort_unstable_by_key(|&(start, _)| start);
    shared
}
