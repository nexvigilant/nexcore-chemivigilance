// Copyright © 2026 NexVigilant LLC. All Rights Reserved.

//! Limitation generation for the chemivigilance pipeline.
//!
//! [`generate_limitations`] always returns at least 3 [`Limitation`] values,
//! satisfying ToV Axiom 4 (d(s) > 0 — every scientific claim carries irreducible
//! uncertainty).  Three standard limitations are unconditionally included;
//! additional conditional limitations are appended based on descriptor and
//! domain characteristics.
//!
//! # Invariant
//!
//! `generate_limitations(..).len() >= 3` — guaranteed by construction.
//!
//! # Examples
//!
//! ```rust
//! use nexcore_chemivigilance::regulatory::generate_limitations;
//! use nexcore_molcore::descriptor::calculate_descriptors;
//! use nexcore_molcore::graph::MolGraph;
//! use nexcore_molcore::smiles::parse;
//! use nexcore_qsar::types::ToxProfile;
//!
//! let mol = parse("CCO").unwrap_or_default();
//! let g = MolGraph::from_molecule(mol);
//! let desc = calculate_descriptors(&g);
//! let profile = ToxProfile::default();
//! let limitations = generate_limitations(&desc, &profile, 0);
//! assert!(limitations.len() >= 3, "must always produce at least 3 limitations");
//! ```

use nexcore_molcore::descriptor::Descriptors;
use nexcore_qsar::types::{DomainStatus, ToxProfile};

use crate::brief::{Limitation, LimitationCategory, LimitationSeverity};

// Threshold above which a heavy-atom count triggers an extra limitation.
const LARGE_MOLECULE_THRESHOLD: usize = 50;

// Threshold above which a high alert count triggers an extra limitation.
const HIGH_ALERT_COUNT_THRESHOLD: usize = 5;

/// Generate methodological and regulatory limitations for a [`crate::brief::SafetyBrief`].
///
/// # Invariant
///
/// The returned `Vec` always contains **at least 3** elements (ToV Axiom 4).
///
/// # Parameters
///
/// - `descriptors` — computed molecular descriptors.
/// - `tox_profile` — QSAR toxicity profile (used for domain status).
/// - `alert_count` — number of structural alerts triggered.
///
/// # Examples
///
/// ```rust
/// use nexcore_chemivigilance::regulatory::generate_limitations;
/// use nexcore_molcore::descriptor::calculate_descriptors;
/// use nexcore_molcore::graph::MolGraph;
/// use nexcore_molcore::smiles::parse;
/// use nexcore_qsar::types::ToxProfile;
///
/// let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap_or_default();
/// let g = MolGraph::from_molecule(mol);
/// let desc = calculate_descriptors(&g);
/// let profile = ToxProfile::default();
/// let lims = generate_limitations(&desc, &profile, 0);
/// assert!(lims.len() >= 3);
/// ```
#[must_use]
pub fn generate_limitations(
    descriptors: &Descriptors,
    tox_profile: &ToxProfile,
    alert_count: usize,
) -> Vec<Limitation> {
    let mut limitations = Vec::with_capacity(5);

    // ------------------------------------------------------------------
    // Standard limitations — always included (satisfies Axiom 4 floor)
    // ------------------------------------------------------------------

    limitations.push(Limitation {
        category: LimitationCategory::ModelScope,
        description: "Predictions are based on rule-based models (Phase 1) — not validated \
                      machine-learning models. Results are indicative, not definitive."
            .to_string(),
        severity: LimitationSeverity::Medium,
    });

    limitations.push(Limitation {
        category: LimitationCategory::DataQuality,
        description: "QSAR predictions rely on molecular descriptors computed from the SMILES \
                      representation. Stereochemistry and 3D conformations are not considered; \
                      enantiomers and diastereomers will receive identical predictions."
            .to_string(),
        severity: LimitationSeverity::Medium,
    });

    limitations.push(Limitation {
        category: LimitationCategory::RegulatoryDisclaimer,
        description: "This report is for informational purposes only and does not constitute \
                      regulatory advice. Formal ICH M7 assessment requires expert review by a \
                      qualified toxicologist and is not replaced by this automated analysis."
            .to_string(),
        severity: LimitationSeverity::High,
    });

    // ------------------------------------------------------------------
    // Conditional limitations
    // ------------------------------------------------------------------

    // Domain applicability warning.
    let domain_warning = match &tox_profile.applicability_domain {
        DomainStatus::OutOfDomain { warning, .. } => Some(format!(
            "The compound falls outside the applicability domain of the QSAR models \
                 ({warning}). Predictions for out-of-domain compounds carry substantially \
                 higher uncertainty and should be treated as preliminary.",
        )),
        DomainStatus::Borderline { warning, .. } => Some(format!(
            "The compound is at the boundary of the QSAR applicability domain ({warning}). \
                 Confidence in endpoint predictions is reduced; independent confirmation \
                 is recommended.",
        )),
        DomainStatus::InDomain { .. } => None,
    };

    if let Some(warning) = domain_warning {
        limitations.push(Limitation {
            category: LimitationCategory::MethodologicalConstraint,
            description: warning,
            severity: LimitationSeverity::High,
        });
    }

    // Large molecule warning.
    if descriptors.heavy_atom_count > LARGE_MOLECULE_THRESHOLD {
        limitations.push(Limitation {
            category: LimitationCategory::MethodologicalConstraint,
            description: format!(
                "The molecule has {} heavy atoms, which exceeds the typical drug-like range \
                 (up to ~50). Rule-based models were trained primarily on drug-like small \
                 molecules; predictions for larger structures (peptides, macrocycles, \
                 polymers) are less reliable.",
                descriptors.heavy_atom_count
            ),
            severity: LimitationSeverity::High,
        });
    }

    // High alert count warning.
    if alert_count > HIGH_ALERT_COUNT_THRESHOLD {
        limitations.push(Limitation {
            category: LimitationCategory::ModelScope,
            description: format!(
                "{alert_count} structural alerts were triggered. A high alert count may reflect \
                 over-sensitive pattern matching rather than a genuinely highly toxic molecule; \
                 each alert should be evaluated in context by a medicinal chemist."
            ),
            severity: LimitationSeverity::Medium,
        });
    }

    // Temporal validity — always add at the end to ensure freshness caveat.
    limitations.push(Limitation {
        category: LimitationCategory::TemporalValidity,
        description: "Alert libraries, toxicological thresholds, and regulatory guidelines \
                      are updated periodically. This report reflects the state of the \
                      Chemivigilance Platform at the time of generation; re-assessment \
                      may be required as knowledge evolves."
            .to_string(),
        severity: LimitationSeverity::Low,
    });

    limitations
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use nexcore_molcore::descriptor::calculate_descriptors;
    use nexcore_molcore::graph::MolGraph;
    use nexcore_molcore::smiles::parse;
    use nexcore_qsar::types::{DomainStatus, PredictionResult, RiskLevel, ToxClass, ToxProfile};

    fn descriptors_for(smiles: &str) -> Descriptors {
        let mol = parse(smiles).unwrap_or_default();
        let g = MolGraph::from_molecule(mol);
        calculate_descriptors(&g)
    }

    fn in_domain_profile() -> ToxProfile {
        let result = PredictionResult {
            probability: 0.1,
            classification: ToxClass::Negative,
            confidence: 0.9,
            in_domain: true,
            model_version: "v1".to_string(),
        };
        ToxProfile {
            mutagenicity: result.clone(),
            hepatotoxicity: result.clone(),
            cardiotoxicity: result,
            off_target_binding: vec![],
            applicability_domain: DomainStatus::InDomain { confidence: 0.9 },
            overall_risk: RiskLevel::Low,
        }
    }

    #[test]
    fn test_always_at_least_three_limitations() {
        let desc = descriptors_for("CCO");
        let profile = in_domain_profile();
        let lims = generate_limitations(&desc, &profile, 0);
        assert!(
            lims.len() >= 3,
            "must always generate at least 3 limitations, got {}",
            lims.len()
        );
    }

    #[test]
    fn test_aspirin_at_least_three_limitations() {
        let desc = descriptors_for("CC(=O)Oc1ccccc1C(=O)O");
        let profile = in_domain_profile();
        let lims = generate_limitations(&desc, &profile, 0);
        assert!(lims.len() >= 3, "aspirin must produce >= 3 limitations");
    }

    #[test]
    fn test_out_of_domain_adds_limitation() {
        let desc = descriptors_for("CCO");
        let mut profile = in_domain_profile();
        profile.applicability_domain = DomainStatus::OutOfDomain {
            distance: 2.0,
            warning: "MW < 100 Da; LogP < -3".to_string(),
        };
        let lims = generate_limitations(&desc, &profile, 0);
        let has_domain = lims
            .iter()
            .any(|l| l.category == LimitationCategory::MethodologicalConstraint);
        assert!(
            has_domain,
            "OutOfDomain status must add a MethodologicalConstraint limitation"
        );
    }

    #[test]
    fn test_borderline_adds_limitation() {
        let desc = descriptors_for("CCO");
        let mut profile = in_domain_profile();
        profile.applicability_domain = DomainStatus::Borderline {
            confidence: 0.6,
            warning: "MW slightly below threshold".to_string(),
        };
        let lims = generate_limitations(&desc, &profile, 0);
        let has_domain = lims
            .iter()
            .any(|l| l.category == LimitationCategory::MethodologicalConstraint);
        assert!(
            has_domain,
            "Borderline domain must add MethodologicalConstraint"
        );
    }

    #[test]
    fn test_high_alert_count_adds_limitation() {
        let desc = descriptors_for("c1ccccc1");
        let profile = in_domain_profile();
        let lims_few = generate_limitations(&desc, &profile, 2);
        let lims_many = generate_limitations(&desc, &profile, 6);
        assert!(
            lims_many.len() > lims_few.len(),
            "more alerts should produce more limitations"
        );
    }

    #[test]
    fn test_standard_limitation_categories_present() {
        let desc = descriptors_for("CC(=O)Oc1ccccc1C(=O)O");
        let profile = in_domain_profile();
        let lims = generate_limitations(&desc, &profile, 0);

        let has_model_scope = lims
            .iter()
            .any(|l| l.category == LimitationCategory::ModelScope);
        let has_data_quality = lims
            .iter()
            .any(|l| l.category == LimitationCategory::DataQuality);
        let has_regulatory = lims
            .iter()
            .any(|l| l.category == LimitationCategory::RegulatoryDisclaimer);

        assert!(
            has_model_scope,
            "ModelScope limitation must always be present"
        );
        assert!(
            has_data_quality,
            "DataQuality limitation must always be present"
        );
        assert!(
            has_regulatory,
            "RegulatoryDisclaimer limitation must always be present"
        );
    }

    #[test]
    fn test_temporal_validity_always_present() {
        let desc = descriptors_for("CCO");
        let profile = in_domain_profile();
        let lims = generate_limitations(&desc, &profile, 0);
        assert!(
            lims.iter()
                .any(|l| l.category == LimitationCategory::TemporalValidity),
            "TemporalValidity limitation must always be present"
        );
    }

    #[test]
    fn test_default_profile_limitations() {
        let desc = descriptors_for("CCO");
        let profile = ToxProfile::default(); // OutOfDomain by default
        let lims = generate_limitations(&desc, &profile, 0);
        assert!(
            lims.len() >= 4,
            "default (OutOfDomain) profile should add the domain constraint, giving >= 4 lims, got {}",
            lims.len()
        );
    }
}
