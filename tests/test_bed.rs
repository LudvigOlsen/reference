#[cfg(test)]
mod tests {
    use reference::reference::bed::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    /// Helper: write a string into a temp BED file and return the handle.
    fn write_bed(contents: &str) -> NamedTempFile {
        let mut file = NamedTempFile::new().expect("create temp file");
        file.write_all(contents.as_bytes())
            .expect("write temp file");
        file
    }

    #[test]
    fn windows_are_loaded_and_sorted() -> anyhow::Result<()> {
        // BED rows intentionally out of order and with a comment
        let bed = "\
# header
chr1\t10\t20
chr1\t0\t5
chr2\t5\t15
";
        let tmp = write_bed(bed);
        let chromosomes = vec!["chr1".into(), "chr2".into()];

        let map = load_windows(tmp.path(), &chromosomes)?;

        // chr1 should hold two windows sorted by (start,end)
        let w1 = &map["chr1"];
        assert_eq!(w1.len(), 2);
        // Sorting flips them; but original indices (2-nd tuple element) remain
        assert_eq!(w1[0], (0, 5, 1)); // earlier start, original win_idx 1
        assert_eq!(w1[1], (10, 20, 0)); // later start, original win_idx 0

        // chr2 has one window with the next running index
        let w2 = &map["chr2"];
        assert_eq!(w2, &vec![(5, 15, 2)]);

        Ok(())
    }

    #[test]
    fn missing_chromosomes_emitted_empty_vec() -> anyhow::Result<()> {
        let bed = "chr1\t0\t10\n";
        let tmp = write_bed(bed);
        let chromosomes = vec!["chr1".into(), "chrX".into()];

        let map = load_windows(tmp.path(), &chromosomes)?;

        assert_eq!(map["chr1"].len(), 1);
        // chrX was requested but absent in BED → empty Vec
        assert!(map["chrX"].is_empty());

        Ok(())
    }

    #[test]
    fn invalid_coordinates_return_error() {
        let bed = "chr1\tstart\t10\n"; // non-numeric start
        let tmp = write_bed(bed);
        let chromosomes = vec!["chr1".into()];

        let err = load_windows(tmp.path(), &chromosomes).unwrap_err();
        assert!(
            err.to_string().contains("Parsing window start"),
            "unexpected error: {err}"
        );
    }
}
