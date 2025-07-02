#[cfg(test)]
mod counting_tests {
    use fxhash::FxHashMap;
    use reference::cli::BigCount;
    use reference::reference::counting::*;
    use reference::reference::kmer_codec::*;
    use smallvec::SmallVec;

    #[test]
    fn basic_dinucleotide_counts() {
        let seq = b"ACGTAC"; // AC CG GT TA AC

        // Build per-test environment ---------------------------------
        let specs = build_kmer_specs(&[2]).unwrap();
        let codes_by_k = build_codes_per_k(seq, &specs);

        let spec2 = &specs[&2];
        let mut encs: SmallVec<[Enc<'_>; 8]> = SmallVec::new();
        encs.push(Enc {
            k: 2,
            codes: &codes_by_k[&2],
            none: spec2.sentinel_none(),
            n: spec2.sentinel_n(),
        });
        // ----------------------------------------------------------------

        let windows = vec![(0, seq.len() as u64, 0)];
        let mut buckets = vec![FxHashMap::<Kmer, BigCount>::default(); windows.len()];

        count_kmers_by_window(&mut buckets, &encs, &windows, seq.len() as u64);

        // Decode -> human-readable
        let mut human: FxHashMap<String, u64> = FxHashMap::default();
        for (kmer, &cnt) in &buckets[0] {
            human.insert(spec2.decode_kmer(kmer.code), cnt);
        }

        assert_eq!(human["AC"], 2);
        assert_eq!(human["CG"], 1);
        assert_eq!(human["GT"], 1);
        assert_eq!(human["TA"], 1);
        assert_eq!(human.len(), 4);
    }

    #[test]
    fn windows_with_n_are_ignored() {
        let seq = b"ACNAC"; // valid AC at 0 and 3

        let specs = build_kmer_specs(&[2]).unwrap();
        let codes_by_k = build_codes_per_k(seq, &specs);
        let spec2 = &specs[&2];

        let mut encs: SmallVec<[Enc<'_>; 8]> = SmallVec::new();
        encs.push(Enc {
            k: 2,
            codes: &codes_by_k[&2],
            none: spec2.sentinel_none(),
            n: spec2.sentinel_n(),
        });

        let windows = vec![(0, seq.len() as u64, 0)];
        let mut buckets = vec![FxHashMap::<Kmer, BigCount>::default(); 1];

        count_kmers_by_window(&mut buckets, &encs, &windows, seq.len() as u64);

        assert_eq!(buckets[0].len(), 1);
        assert_eq!(buckets[0].values().copied().sum::<u64>(), 2);
    }

    #[test]
    fn multiple_windows_independent() {
        let seq = b"AAAA"; // all 2-mers = AA

        let specs = build_kmer_specs(&[2]).unwrap();
        let codes_by_k = build_codes_per_k(seq, &specs);
        let spec2 = &specs[&2];

        let mut encs: SmallVec<[Enc<'_>; 8]> = SmallVec::new();
        encs.push(Enc {
            k: 2,
            codes: &codes_by_k[&2],
            none: spec2.sentinel_none(),
            n: spec2.sentinel_n(),
        });

        let windows = vec![(0, 2, 0), (2, 4, 1)]; // two half-windows
        let mut buckets = vec![FxHashMap::<Kmer, BigCount>::default(); windows.len()];

        count_kmers_by_window(&mut buckets, &encs, &windows, seq.len() as u64);

        for bucket in buckets {
            assert_eq!(bucket.values().copied().sum::<u64>(), 1);
        }
    }

    // Window shorter than k
    #[test]
    fn window_shorter_than_k_yields_zero() {
        // Sequence is 4 bp, k-mer length is 6
        let seq = b"ACGTANA";

        let specs = build_kmer_specs(&[6]).unwrap();
        let codes_by_k = build_codes_per_k(seq, &specs);
        let spec6 = &specs[&6];

        let mut encs: SmallVec<[Enc<'_>; 8]> = SmallVec::new();
        encs.push(Enc {
            k: 6,
            codes: &codes_by_k[&6],
            none: spec6.sentinel_none(),
            n: spec6.sentinel_n(),
        });

        let windows = vec![(0, 4, 0)]; // 4-bp window
        let mut buckets = vec![FxHashMap::<Kmer, BigCount>::default(); 1];

        count_kmers_by_window(&mut buckets, &encs, &windows, seq.len() as u64);

        assert!(buckets[0].is_empty());
    }

    // Chromosome shorter than k
    #[test]
    fn chromosome_shorter_than_k_yields_zero() {
        // Whole chromosome is 2 bp, k-mer length is 3
        let seq = b"AC";

        let specs = build_kmer_specs(&[3]).unwrap();
        let codes_by_k = build_codes_per_k(seq, &specs);
        let spec3 = &specs[&3];

        let mut encs: SmallVec<[Enc<'_>; 8]> = SmallVec::new();
        encs.push(Enc {
            k: 3,
            codes: &codes_by_k[&3],
            none: spec3.sentinel_none(),
            n: spec3.sentinel_n(),
        });

        let windows = vec![(0, 2, 0)];
        let mut buckets = vec![FxHashMap::<Kmer, BigCount>::default(); 1];

        count_kmers_by_window(&mut buckets, &encs, &windows, seq.len() as u64);

        assert!(buckets[0].is_empty());
    }

    // Window length exactly k
    #[test]
    fn window_exactly_k_has_one_kmer() {
        // Sequence is 4 bp, k-mer length is 4
        let seq = b"ACGT";

        let specs = build_kmer_specs(&[4]).unwrap();
        let codes_by_k = build_codes_per_k(seq, &specs);
        let spec4 = &specs[&4];

        let mut encs: SmallVec<[Enc<'_>; 8]> = SmallVec::new();
        encs.push(Enc {
            k: 4,
            codes: &codes_by_k[&4],
            none: spec4.sentinel_none(),
            n: spec4.sentinel_n(),
        });

        let windows = vec![(0, 4, 0)];
        let mut buckets = vec![FxHashMap::<Kmer, BigCount>::default(); 1];

        count_kmers_by_window(&mut buckets, &encs, &windows, seq.len() as u64);

        // Exactly one k-mer counted
        assert_eq!(buckets[0].values().copied().sum::<u64>(), 1);
        assert_eq!(buckets[0].len(), 1);

        // Decoded motif equals the sequence
        let motif = spec4.decode_kmer(buckets[0].keys().next().unwrap().code);
        assert_eq!(motif.as_bytes(), seq);
    }

    // Tail window that contains only sentinel positions
    #[test]
    fn tail_window_with_only_sentinels_is_empty() {
        // Sequence is 6 bp, k-mer length is 3
        let seq = b"ACGTAC";

        let specs = build_kmer_specs(&[3]).unwrap();
        let codes_by_k = build_codes_per_k(seq, &specs);
        let spec3 = &specs[&3];

        let mut encs: SmallVec<[Enc<'_>; 8]> = SmallVec::new();
        encs.push(Enc {
            k: 3,
            codes: &codes_by_k[&3],
            none: spec3.sentinel_none(),
            n: spec3.sentinel_n(),
        });

        // Start inside the last k-1 bases; no full k-mer fits
        let start = seq.len() as u64 - 2;
        let windows = vec![(start, seq.len() as u64, 0)];
        let mut buckets = vec![FxHashMap::<Kmer, BigCount>::default(); 1];

        count_kmers_by_window(&mut buckets, &encs, &windows, seq.len() as u64);

        assert!(buckets[0].is_empty());
    }
}
