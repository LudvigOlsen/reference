#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use fxhash::FxHashMap;
    use reference::reference::kmer_codec::*;
    use reference::reference::process_counts::*;

    /* --------------------------------------------------------------------- */
    /*  collapse_map / canonical                                            */
    /* --------------------------------------------------------------------- */

    #[test]
    fn revcomp_and_canonical_roundtrip() {
        // Palindrome stays identical
        let pal = "ACGT";
        assert_eq!(revcomp(pal), pal);
        assert_eq!(canonical(pal.to_string()), pal);

        // Non‑palindrome collapses to lexicographically smaller string
        let fwd = "ACG"; // rc == "CGT"
        let rc = revcomp(fwd);
        assert_eq!(rc, "CGT");
        assert!(fwd < &rc);
        assert_eq!(canonical(fwd.to_string()), fwd);
        assert_eq!(canonical(rc), fwd); // canonical of rc collapses back

        assert_eq!(canonical("AC".into()), "AC"); // AC vs GT  → AC
        assert_eq!(canonical("GT".into()), "AC"); // GT vs AC  → AC
    }

    #[test]
    fn collapse_map_sums_reverse_complements() {
        let mut m: FxHashMap<String, u64> = FxHashMap::default();
        m.insert("ACG".into(), 2);
        m.insert("CGT".into(), 3); // reverse complement of ACG
        let collapsed = collapse_map(&m);
        assert_eq!(collapsed.len(), 1);
        assert_eq!(collapsed["ACG"], 5);
    }

    /* --------------------------------------------------------------------- */
    /*  encode_base / choose_width                                          */
    /* --------------------------------------------------------------------- */

    #[test]
    fn encode_base_matches_spec() {
        assert_eq!(encode_base(b'A'), 0);
        assert_eq!(encode_base(b'C'), 1);
        assert_eq!(encode_base(b'G'), 2);
        assert_eq!(encode_base(b'T'), 3);
        assert_eq!(encode_base(b'N'), 4);
        assert_eq!(encode_base(b'X'), 4); // unknown → 4
    }

    #[test]
    fn choose_width_returns_correct_sentinals() {
        // k = 3 → 5^3 = 125 < 254 so fits in u8
        let (w, none, n) = choose_width(3).unwrap();
        assert_eq!(w, Width::U8);
        assert_eq!(none, u8::MAX as u64);
        assert_eq!(n, (u8::MAX - 1) as u64);

        // k = 10 → 5^10 ≈ 9.7e6 fits in u32
        let (w, _, _) = choose_width(10).unwrap();
        assert_eq!(w, Width::U32);

        // Large k escalates to u64
        assert_eq!(choose_width(26).unwrap().0, Width::U64);
    }

    /* --------------------------------------------------------------------- */
    /*  build_codes / decode_kmer round‑trip                                */
    /* --------------------------------------------------------------------- */

    #[test]
    fn build_and_decode_roundtrip_generated() {
        // Simple sequence without Ns, k = 2
        let seq = b"ACGTAC";
        let spec = build_kmer_specs(&[2]).unwrap().remove(&2u8).unwrap();

        // Build per‑position codes and decode back
        let codes = spec.build_codes(seq);
        for (i, &code) in codes.iter().enumerate() {
            let decoded = spec.decode_kmer(code);
            let expected = if i + 2 <= seq.len() {
                // full window fits
                std::str::from_utf8(&seq[i..i + 2]).unwrap()
            } else {
                "NN" // sentinel_none area
            };
            assert_eq!(decoded, expected);
        }
    }

    #[test]
    fn build_and_decode_roundtrip_hardcoded() {
        // Build spec for k = 3
        let spec = build_kmer_specs(&[3]).unwrap().remove(&3u8).unwrap();

        let seq = b"ACGTACN";
        let codes = spec.build_codes(seq);

        // Position-by-position expectations
        // idx : window   -> code   -> decoded
        //  0  : ACG      -> ok
        //  1  : CGT      -> ok
        //  2  : GTA      -> ok
        //  3  : TAC      -> ok
        //  4  : ACN      -> sentinel_n
        //  5  : (no full window) -> sentinel_none
        //  6  : (no full window) -> sentinel_none
        for (i, &code) in codes.iter().enumerate() {
            let decoded = spec.decode_kmer(code);
            match i {
                0 => assert_eq!(decoded, "ACG"),
                1 => assert_eq!(decoded, "CGT"),
                2 => assert_eq!(decoded, "GTA"),
                3 => assert_eq!(decoded, "TAC"),
                4 => assert_eq!(decoded, "NNN"), // contains N
                _ => assert_eq!(decoded, "NNN"), // sentinel_none tail
            }
        }
    }

    /* --------------------------------------------------------------------- */
    /*  merge_decoded_counts                                                */
    /* --------------------------------------------------------------------- */

    #[test]
    fn merge_decoded_counts_sums_bins() {
        // Build two DecodedCounts with cross‑over motifs
        let mut dc1 = DecodedCounts {
            counts: HashMap::new(),
        };
        let mut dc2 = DecodedCounts {
            counts: HashMap::new(),
        };

        dc1.counts
            .insert(3, FxHashMap::from_iter([(String::from("AAA"), 1u64)]));
        dc2.counts
            .insert(3, FxHashMap::from_iter([(String::from("AAA"), 4u64)]));
        dc2.counts
            .get_mut(&3)
            .unwrap()
            .insert(String::from("CCC"), 2u64);

        let merged = merge_decoded_counts(vec![dc1, dc2]);
        let bucket = &merged.counts[&3];
        assert_eq!(bucket["AAA"], 5);
        assert_eq!(bucket["CCC"], 2);
    }

    /* --------------------------------------------------------------------- */
    /*  all_motifs                                                           */
    /* --------------------------------------------------------------------- */

    #[test]
    fn all_motifs_returns_full_space_for_k_up_to_6() {
        let specs = build_kmer_specs(&[2]).unwrap();
        let motifs = all_motifs(2, &specs);
        // 4^2 = 16 motifs, none with N
        assert_eq!(motifs.len(), 16);
        assert!(motifs.contains(&"AA".to_string()));
        assert!(motifs.contains(&"TT".to_string()));
    }

    /* --------------------------------------------------------------------- */
    /*  prepare_decoded_counts high-level path                               */
    /* --------------------------------------------------------------------- */

    #[test]
    fn prepare_decoded_counts_outputs_expected_structure() {
        // Two windows with a single 2-mer each
        let specs = build_kmer_specs(&[7]).unwrap();
        let mut win1 = DecodedCounts {
            counts: HashMap::new(),
        };
        win1.counts
            .insert(7, FxHashMap::from_iter([(String::from("AAAAAAA"), 1u64)]));
        let mut win2 = DecodedCounts {
            counts: HashMap::new(),
        };
        // NOTE: 7-mer so it doesn't add all motifs!
        win2.counts
            .insert(7, FxHashMap::from_iter([(String::from("CCCCCCC"), 1u64)]));

        let (prepared, per_k_motifs) =
            prepare_decoded_counts(&[win1.clone(), win2.clone()], false, &specs);

        // Motifs list contains both AA and CC, sorted
        assert_eq!(
            per_k_motifs[&7],
            vec!["AAAAAAA".to_string(), "CCCCCCC".to_string()]
        );

        // Prepared bins kept original counts
        assert_eq!(prepared[0].counts[&7]["AAAAAAA"], 1);
        assert_eq!(prepared[1].counts[&7]["CCCCCCC"], 1);
    }
}
