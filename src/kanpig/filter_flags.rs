use bitflags::bitflags;

bitflags! {
    pub struct FiltFlags: u32 {
        const PASS       = 0x0;  // passing
        const GTMISMATCH = 0x1;  // genotype from AD doesn't match path genotype
        const LOWGQ      = 0x2;  // genotype quality below 5
        const LOWCOV     = 0x4;  // coverage below 5
        const LOWSQ      = 0x8;  // sample quality below below 5
        const LOWALT     = 0x16; // alt coverage below 5
        const PARTIAL    = 0x32; // best scoring path only used part of the haplotype
    }
}
