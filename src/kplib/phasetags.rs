/// Utility for sorting HP tags 1,2 which can sometimes be None
/// paths.sort_by(|a, b| hp_sorter(&a.hp, &b.hp));
pub fn hp_sorter(a: &Option<u8>, b: &Option<u8>) -> std::cmp::Ordering {
    match (a, b) {
        // If both are Some, compare
        (Some(va), Some(vb)) => vb.cmp(va),

        // If one is None and the other is Some(1), None comes first
        (Some(1), None) => std::cmp::Ordering::Greater,
        (None, Some(1)) => std::cmp::Ordering::Less,

        // If one is None and the other is Some(2), None comes last
        (Some(2), None) => std::cmp::Ordering::Less,
        (None, Some(2)) => std::cmp::Ordering::Greater,

        // If both are None, do nothing
        (None, None) => std::cmp::Ordering::Equal,

        // Default fallback (not strictly needed with the cases above)
        _ => std::cmp::Ordering::Equal,
    }
}
