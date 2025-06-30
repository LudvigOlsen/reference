#[derive(Debug, Default)]
pub struct MotifExtractionCounters {
    pub total: u64,
    pub accepted: u64,
    pub left: u64,
    pub right_mate: u64,
    pub blacklisted: u64,
    pub left_clipped: u64,
    pub right_clipped: u64,
    pub left_forward: u64,
    pub left_reverse: u64,
    pub right_forward: u64,
    pub right_reverse: u64,
    pub gc_excl: u64,
    pub counted: u64,
}

impl std::ops::AddAssign for MotifExtractionCounters {
    fn add_assign(&mut self, other: Self) {
        self.total += other.total;
        self.accepted += other.accepted;
        self.left += other.left;
        self.right_mate += other.right_mate;
        self.blacklisted += other.blacklisted;
        self.left_clipped += other.left_clipped;
        self.right_clipped += other.right_clipped;
        self.left_forward += other.left_forward;
        self.left_reverse += other.left_reverse;
        self.right_forward += other.right_forward;
        self.right_reverse += other.right_reverse;
        self.gc_excl += other.gc_excl;
        self.counted += other.counted;
    }
}

#[derive(Debug, Default)]
pub struct FragsizeExtractionCounters {
    pub total: u64,
    pub accepted: u64,
    pub blacklisted: u64,
    pub gc_excl: u64,
    pub counted: u64,
}

impl std::ops::AddAssign for FragsizeExtractionCounters {
    fn add_assign(&mut self, other: Self) {
        self.total += other.total;
        self.accepted += other.accepted;
        self.blacklisted += other.blacklisted;
        self.gc_excl += other.gc_excl;
        self.counted += other.counted;
    }
}

#[derive(Debug, Default)]
pub struct RefKmerExtractionCounters {
    pub total: u64,
    pub blacklisted: u64,
    pub ambiguous: u64,
    pub counted: u64,
}

impl std::ops::AddAssign for RefKmerExtractionCounters {
    fn add_assign(&mut self, other: Self) {
        self.total += other.total;
        self.blacklisted += other.blacklisted;
        self.ambiguous += other.ambiguous;
        self.counted += other.counted;
    }
}

#[derive(Debug, Default)]
pub struct FastqMersExtractionCounters {
    pub total: u64,
    pub ambiguous: u64,
    pub counted: u64,
}

impl std::ops::AddAssign for FastqMersExtractionCounters {
    fn add_assign(&mut self, other: Self) {
        self.total += other.total;
        self.ambiguous += other.ambiguous;
        self.counted += other.counted;
    }
}

#[derive(Debug, Default)]
pub struct ConsensusDepthCounters {
    pub total: u64,
    pub accepted: u64,
    pub left: u64,
    pub right_mate: u64,
    pub gc_excl: u64,
    pub missing_md: u64,
    pub counted: u64,
}

impl std::ops::AddAssign for ConsensusDepthCounters {
    fn add_assign(&mut self, other: Self) {
        self.total += other.total;
        self.accepted += other.accepted;
        self.left += other.left;
        self.right_mate += other.right_mate;
        self.gc_excl += other.gc_excl;
        self.missing_md += other.missing_md;
        self.counted += other.counted;
    }
}
