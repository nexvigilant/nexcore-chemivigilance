// Copyright © 2026 NexVigilant LLC. All Rights Reserved.

//! Crate-level error types for the chemivigilance pipeline.
//!
//! All downstream errors are wrapped in [`ChemivigilanceError`] so callers
//! interact with a single unified error surface.
//!
//! # Examples
//!
//! ```rust
//! use nexcore_chemivigilance::error::ChemivigilanceError;
//!
//! let e = ChemivigilanceError::InvalidInput("empty SMILES".to_string());
//! assert!(format!("{e}").contains("empty SMILES"));
//! ```

use nexcore_error::Error;

/// Errors produced by the chemivigilance pipeline.
#[derive(Debug, Error)]
pub enum ChemivigilanceError {
    /// The upstream `nexcore-molcore` crate returned an error.
    #[error("molcore error: {0}")]
    MolcoreError(String),

    /// The upstream `nexcore-structural-alerts` crate returned an error.
    #[error("structural alert error: {0}")]
    AlertError(String),

    /// The upstream `nexcore-qsar` crate returned an error.
    #[error("QSAR error: {0}")]
    QsarError(String),

    /// The upstream `nexcore-metabolite` crate returned an error.
    #[error("metabolite error: {0}")]
    MetaboliteError(String),

    /// The caller supplied invalid input (e.g. empty or malformed SMILES).
    #[error("invalid input: {0}")]
    InvalidInput(String),

    /// Fewer than 3 limitations were generated — ToV Axiom 4 violation.
    ///
    /// This variant should never appear when the code is correct; it exists as
    /// an explicit guard so that the invariant is type-visible.
    #[error("insufficient limitations: generated {0}, required at least 3 (ToV Axiom 4)")]
    InsufficientLimitations(usize),
}

impl From<nexcore_structural_alerts::AlertError> for ChemivigilanceError {
    fn from(e: nexcore_structural_alerts::AlertError) -> Self {
        Self::AlertError(e.to_string())
    }
}

impl From<nexcore_qsar::QsarError> for ChemivigilanceError {
    fn from(e: nexcore_qsar::QsarError) -> Self {
        Self::QsarError(e.to_string())
    }
}

impl From<nexcore_metabolite::MetaboliteError> for ChemivigilanceError {
    fn from(e: nexcore_metabolite::MetaboliteError) -> Self {
        Self::MetaboliteError(e.to_string())
    }
}

/// Convenient `Result` alias for the chemivigilance pipeline.
pub type ChemivigilanceResult<T> = Result<T, ChemivigilanceError>;

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_invalid_input_display() {
        let e = ChemivigilanceError::InvalidInput("bad smiles".to_string());
        assert!(format!("{e}").contains("bad smiles"));
    }

    #[test]
    fn test_insufficient_limitations_display() {
        let e = ChemivigilanceError::InsufficientLimitations(1);
        assert!(format!("{e}").contains('1'));
        assert!(format!("{e}").contains('3'));
    }

    #[test]
    fn test_alert_error_display() {
        let e = ChemivigilanceError::AlertError("pattern invalid".to_string());
        assert!(format!("{e}").contains("pattern invalid"));
    }

    #[test]
    fn test_molcore_error_display() {
        let e = ChemivigilanceError::MolcoreError("parse failed".to_string());
        assert!(format!("{e}").contains("parse failed"));
    }

    #[test]
    fn test_qsar_error_display() {
        let e = ChemivigilanceError::QsarError("out of range".to_string());
        assert!(format!("{e}").contains("out of range"));
    }

    #[test]
    fn test_metabolite_error_display() {
        let e = ChemivigilanceError::MetaboliteError("no sites".to_string());
        assert!(format!("{e}").contains("no sites"));
    }

    #[test]
    fn test_from_alert_error() {
        let alert_err = nexcore_structural_alerts::AlertError::SmilesParse("bad".to_string());
        let cv_err = ChemivigilanceError::from(alert_err);
        assert!(matches!(cv_err, ChemivigilanceError::AlertError(_)));
    }

    #[test]
    fn test_from_metabolite_error() {
        let meta_err = nexcore_metabolite::MetaboliteError::SmilesParse("bad".to_string());
        let cv_err = ChemivigilanceError::from(meta_err);
        assert!(matches!(cv_err, ChemivigilanceError::MetaboliteError(_)));
    }
}
