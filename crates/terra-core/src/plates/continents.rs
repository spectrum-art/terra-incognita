//! Crust-type definitions for the rebuilt plate system.
//!
//! Note: age-based crust classification was removed in the plate system rebuild.
//! See git history before commit `eb343e4` for the previous implementation.

/// Classification of a grid cell's crust type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CrustType {
    Oceanic,
    Continental,
    ActiveMargin,
    PassiveMargin,
}

/// Returns `true` if the grid cell is any form of continental crust.
pub fn is_continental(crust: CrustType) -> bool {
    matches!(
        crust,
        CrustType::Continental | CrustType::ActiveMargin | CrustType::PassiveMargin
    )
}

/// Convenience: check the `Vec<CrustType>` directly.
pub fn is_continental_cell(crust_field: &[CrustType], idx: usize) -> bool {
    is_continental(crust_field[idx])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn is_continental_helper() {
        assert!(is_continental(CrustType::Continental));
        assert!(is_continental(CrustType::ActiveMargin));
        assert!(is_continental(CrustType::PassiveMargin));
        assert!(!is_continental(CrustType::Oceanic));
    }

    #[test]
    fn crust_variants_are_distinct() {
        assert_ne!(CrustType::Oceanic, CrustType::Continental);
        assert_ne!(CrustType::ActiveMargin, CrustType::PassiveMargin);
    }
}
