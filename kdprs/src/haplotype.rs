pub struct Haplotype {
    kfeat: Vec<f32>,
    size: i64,
    n: u64,
    coverage: u64,
}

// Can I do a ::new but also call it directly via Haplotype { a, b, c,.. }
