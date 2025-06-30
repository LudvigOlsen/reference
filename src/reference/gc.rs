// TODO: What about N's?

/// Build a prefix-sum (cumsum) of GC counts over a byte slice
/// pref[i] = # of G/C in seq[0..i], so pref.len() == seq.len()+1
pub fn build_gc_prefix(seq: &[u8]) -> Vec<u32> {
    let mut pref = Vec::with_capacity(seq.len() + 1);
    pref.push(0);
    for &b in seq {
        let inc = match b {
            b'G' | b'g' | b'C' | b'c' => 1,
            _ => 0,
        };
        // safe to unwrap: we always have at least one element
        let last = *pref.last().unwrap();
        pref.push(last + inc);
    }
    pref
}
