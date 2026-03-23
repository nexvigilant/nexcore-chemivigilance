// Copyright © 2026 NexVigilant LLC. All Rights Reserved.

//! Primary output types for the chemivigilance pipeline.
//!
//! [`SafetyBrief`] is the top-level report produced by
//! [`crate::pipeline::generate_safety_brief`].  It aggregates descriptor data,
//! structural alert counts, QSAR toxicity predictions, metabolite predictions,
//! regulatory flags, and the mandatory limitations required by ToV Axiom 4.
//!
//! # Examples
//!
//! ```rust
//! use nexcore_chemivigilance::brief::{AlertSummary, LimitationCategory, LimitationSeverity};
//! use serde_json;
//!
//! let summary = AlertSummary {
//!     alert_id: "M7-001".to_string(),
//!     alert_name: "Aromatic amine".to_string(),
//!     category: "Mutagenicity".to_string(),
//!     match_count: 1,
//!     confidence: 0.95,
//! };
//! let json = serde_json::to_string(&summary).unwrap_or_default();
//! assert!(json.contains("M7-001"));
//! ```

use nexcore_chrono::DateTime;
use nexcore_metabolite::types::MetaboliteTree;
use nexcore_molcore::descriptor::Descriptors;
use nexcore_qsar::types::{RiskLevel, ToxProfile};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

// ---------------------------------------------------------------------------
// DescriptorSummary
// ---------------------------------------------------------------------------

/// Serialisable mirror of [`nexcore_molcore::descriptor::Descriptors`].
///
/// [`Descriptors`] does not implement `Serialize`/`Deserialize`; this struct
/// copies every field so the values can be embedded in [`SafetyBrief`] and
/// round-tripped through JSON.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DescriptorSummary {
    /// Molecular weight in Da.
    pub molecular_weight: f64,
    /// Wildman-Crippen LogP estimate.
    pub logp: f64,
    /// Topological polar surface area (Angstrom squared).
    pub tpsa: f64,
    /// Hydrogen bond acceptor count.
    pub hba: u8,
    /// Hydrogen bond donor count.
    pub hbd: u8,
    /// Number of rotatable bonds.
    pub rotatable_bonds: u8,
    /// Total ring count (SSSR).
    pub num_rings: u8,
    /// Aromatic ring count.
    pub num_aromatic_rings: u8,
    /// Number of heavy (non-hydrogen) atoms.
    pub heavy_atom_count: usize,
}

impl From<&Descriptors> for DescriptorSummary {
    fn from(d: &Descriptors) -> Self {
        Self {
            molecular_weight: d.molecular_weight,
            logp: d.logp,
            tpsa: d.tpsa,
            hba: d.hba,
            hbd: d.hbd,
            rotatable_bonds: d.rotatable_bonds,
            num_rings: d.num_rings,
            num_aromatic_rings: d.num_aromatic_rings,
            heavy_atom_count: d.heavy_atom_count,
        }
    }
}

// ---------------------------------------------------------------------------
// AlertSummary
// ---------------------------------------------------------------------------

/// A compact, serialisable summary of a single structural alert match.
///
/// [`nexcore_structural_alerts::AlertMatch`] is not `Serialize`; this type
/// carries the fields needed for the [`SafetyBrief`] report.
///
/// # Examples
///
/// ```rust
/// use nexcore_chemivigilance::brief::AlertSummary;
///
/// let s = AlertSummary {
///     alert_id: "M7-002".to_string(),
///     alert_name: "Nitro group".to_string(),
///     category: "Mutagenicity".to_string(),
///     match_count: 1,
///     confidence: 0.9,
/// };
/// assert_eq!(s.match_count, 1);
/// ```
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlertSummary {
    /// Unique identifier of the triggered alert.
    pub alert_id: String,
    /// Human-readable alert name.
    pub alert_name: String,
    /// Toxicological category (e.g. `"Mutagenicity"`).
    pub category: String,
    /// Number of non-overlapping pattern occurrences found.
    pub match_count: usize,
    /// Confidence score in `[0.0, 1.0]` from the alert definition.
    pub confidence: f64,
}

// ---------------------------------------------------------------------------
// Limitation
// ---------------------------------------------------------------------------

/// A single methodological or regulatory limitation attached to a
/// [`SafetyBrief`].
///
/// At least 3 limitations are always generated (ToV Axiom 4: d(s) > 0).
///
/// # Examples
///
/// ```rust
/// use nexcore_chemivigilance::brief::{Limitation, LimitationCategory, LimitationSeverity};
///
/// let lim = Limitation {
///     category: LimitationCategory::ModelScope,
///     description: "Rule-based predictions only.".to_string(),
///     severity: LimitationSeverity::Medium,
/// };
/// assert_eq!(lim.severity, LimitationSeverity::Medium);
/// ```
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Limitation {
    /// Category of the limitation.
    pub category: LimitationCategory,
    /// Human-readable description.
    pub description: String,
    /// Severity of the limitation's impact on result reliability.
    pub severity: LimitationSeverity,
}

/// High-level category for a [`Limitation`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum LimitationCategory {
    /// The model does not cover this chemical class or property space.
    ModelScope,
    /// Input data quality affects prediction reliability.
    DataQuality,
    /// A fundamental constraint of the computational method.
    MethodologicalConstraint,
    /// Predictions may become stale as scientific knowledge advances.
    TemporalValidity,
    /// Formal regulatory disclaimer required for compliance.
    RegulatoryDisclaimer,
}

/// Severity of a [`Limitation`]'s impact on result reliability.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum LimitationSeverity {
    /// Minor impact; result is still highly reliable.
    Low,
    /// Moderate impact; interpretation requires care.
    Medium,
    /// Significant impact; expert review recommended.
    High,
    /// Severe impact; results should not be used without expert validation.
    Critical,
}

// ---------------------------------------------------------------------------
// RegulatoryFlag
// ---------------------------------------------------------------------------

/// A regulatory or safety flag raised during pipeline execution.
///
/// # Examples
///
/// ```rust
/// use nexcore_chemivigilance::brief::{RegulatoryFlag, RegulatoryFlagType};
///
/// let flag = RegulatoryFlag {
///     flag_type: RegulatoryFlagType::IchM7Alert,
///     description: "ICH M7 structural alert triggered.".to_string(),
///     reference: "ICH M7(R1)".to_string(),
/// };
/// assert_eq!(flag.flag_type, RegulatoryFlagType::IchM7Alert);
/// ```
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RegulatoryFlag {
    /// Classification of the regulatory concern.
    pub flag_type: RegulatoryFlagType,
    /// Human-readable description of the concern.
    pub description: String,
    /// Regulatory reference document (e.g. `"ICH M7(R1)"`).
    pub reference: String,
}

/// Type of regulatory or safety flag.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum RegulatoryFlagType {
    /// An ICH M7 structural alert was triggered.
    IchM7Alert,
    /// Mutagenicity probability is at or above the concern threshold.
    HighMutagenicity,
    /// Hepatotoxicity probability is at or above the concern threshold.
    HighHepatotoxicity,
    /// Cardiotoxicity probability is at or above the concern threshold.
    HighCardiotoxicity,
    /// Reactive intermediates were predicted during metabolite analysis.
    ReactiveIntermediate,
    /// The compound falls outside the applicability domain of the models.
    OutOfDomain,
}

// ---------------------------------------------------------------------------
// SafetyBrief
// ---------------------------------------------------------------------------

/// Complete safety report — the primary output of the chemivigilance pipeline.
///
/// Produced by [`crate::pipeline::generate_safety_brief`].
///
/// # Invariant
///
/// `limitations.len() >= 3` — enforced by [`crate::regulatory::generate_limitations`]
/// and validated in the pipeline before the brief is returned.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SafetyBrief {
    /// Original input string (SMILES or compound name).
    pub input: String,
    /// Resolved SMILES used for all calculations.
    pub smiles: String,
    /// Molecular formula (e.g. `"C9H8O4"`).
    pub molecular_formula: String,
    /// Serialisable molecular descriptor summary.
    pub descriptors: DescriptorSummary,
    /// Total number of structural alerts triggered.
    pub alert_count: usize,
    /// Compact summary of each triggered alert.
    pub alert_summary: Vec<AlertSummary>,
    /// QSAR toxicity profile.
    pub tox_profile: ToxProfile,
    /// Predicted metabolite tree.
    pub metabolite_tree: MetaboliteTree,
    /// Composite risk score in `[0.0, 1.0]`.
    pub overall_risk_score: f64,
    /// Risk level classification derived from `overall_risk_score`.
    pub risk_level: RiskLevel,
    /// Methodological and regulatory limitations (always >= 3).
    pub limitations: Vec<Limitation>,
    /// Regulatory flags raised during analysis.
    pub regulatory_flags: Vec<RegulatoryFlag>,
    /// UTC timestamp of report generation.
    pub generated_at: DateTime,
    /// Model version identifiers keyed by component name.
    pub model_versions: HashMap<String, String>,
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_alert_summary_serialization() {
        let summary = AlertSummary {
            alert_id: "M7-001".to_string(),
            alert_name: "Aromatic amine".to_string(),
            category: "Mutagenicity".to_string(),
            match_count: 2,
            confidence: 0.95,
        };
        let json = serde_json::to_string(&summary).unwrap_or_default();
        assert!(json.contains("M7-001"), "id must appear in JSON");
        assert!(
            json.contains("Mutagenicity"),
            "category must appear in JSON"
        );

        let round_tripped: AlertSummary =
            serde_json::from_str(&json).unwrap_or_else(|_| AlertSummary {
                alert_id: String::new(),
                alert_name: String::new(),
                category: String::new(),
                match_count: 0,
                confidence: 0.0,
            });
        assert_eq!(round_tripped.alert_id, "M7-001");
        assert_eq!(round_tripped.match_count, 2);
    }

    #[test]
    fn test_limitation_severity_ordering() {
        assert!(LimitationSeverity::Low < LimitationSeverity::Medium);
        assert!(LimitationSeverity::Medium < LimitationSeverity::High);
        assert!(LimitationSeverity::High < LimitationSeverity::Critical);
    }

    #[test]
    fn test_regulatory_flag_type_eq() {
        assert_eq!(
            RegulatoryFlagType::IchM7Alert,
            RegulatoryFlagType::IchM7Alert
        );
        assert_ne!(
            RegulatoryFlagType::IchM7Alert,
            RegulatoryFlagType::HighMutagenicity
        );
    }

    #[test]
    fn test_limitation_category_eq() {
        assert_eq!(
            LimitationCategory::ModelScope,
            LimitationCategory::ModelScope
        );
        assert_ne!(
            LimitationCategory::ModelScope,
            LimitationCategory::DataQuality
        );
    }

    #[test]
    fn test_descriptor_summary_from_descriptors() {
        use nexcore_molcore::descriptor::calculate_descriptors;
        use nexcore_molcore::graph::MolGraph;
        use nexcore_molcore::smiles::parse;

        let mol = parse("CCO").unwrap_or_default();
        let graph = MolGraph::from_molecule(mol);
        let desc = calculate_descriptors(&graph);
        let summary = DescriptorSummary::from(&desc);

        assert_eq!(summary.heavy_atom_count, desc.heavy_atom_count);
        assert!((summary.molecular_weight - desc.molecular_weight).abs() < f64::EPSILON);
    }

    #[test]
    fn test_descriptor_summary_serialization() {
        let summary = DescriptorSummary {
            molecular_weight: 180.16,
            logp: 1.2,
            tpsa: 63.6,
            hba: 4,
            hbd: 1,
            rotatable_bonds: 3,
            num_rings: 1,
            num_aromatic_rings: 1,
            heavy_atom_count: 13,
        };
        let json = serde_json::to_string(&summary).unwrap_or_default();
        assert!(json.contains("180.16"));
        assert!(json.contains("heavy_atom_count"));
    }
}
