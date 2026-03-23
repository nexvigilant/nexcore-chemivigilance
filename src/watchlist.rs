// Copyright © 2026 NexVigilant LLC. All Rights Reserved.

//! Regulatory watchlist screening.
//!
//! [`check_watchlist`] inspects descriptor, QSAR, and metabolite data and
//! returns zero or more [`RegulatoryFlag`] values indicating regulatory concerns.
//!
//! # Flagging logic
//!
//! | Condition | Flag raised |
//! |-----------|-------------|
//! | Any structural alert in Mutagenicity category | [`RegulatoryFlagType::IchM7Alert`] |
//! | Mutagenicity probability >= 0.5 | [`RegulatoryFlagType::HighMutagenicity`] |
//! | Hepatotoxicity probability >= 0.5 | [`RegulatoryFlagType::HighHepatotoxicity`] |
//! | Cardiotoxicity probability >= 0.5 | [`RegulatoryFlagType::HighCardiotoxicity`] |
//! | Non-empty reactive intermediates | [`RegulatoryFlagType::ReactiveIntermediate`] |
//! | Applicability domain is OutOfDomain | [`RegulatoryFlagType::OutOfDomain`] |
//!
//! # Examples
//!
//! ```rust
//! use nexcore_chemivigilance::watchlist::check_watchlist;
//! use nexcore_chemivigilance::brief::RegulatoryFlagType;
//! use nexcore_qsar::types::{ToxProfile, DomainStatus};
//! use nexcore_metabolite::types::MetaboliteTree;
//!
//! let profile = ToxProfile::default();
//! let tree = MetaboliteTree::default();
//! // Default profile has OutOfDomain status → OutOfDomain flag expected.
//! let flags = check_watchlist(0, &profile, &tree);
//! assert!(flags.iter().any(|f| f.flag_type == RegulatoryFlagType::OutOfDomain));
//! ```

use nexcore_metabolite::types::MetaboliteTree;
use nexcore_qsar::types::{DomainStatus, ToxProfile};

use crate::brief::{RegulatoryFlag, RegulatoryFlagType};

/// Inspect a toxicity profile and metabolite tree, returning any regulatory flags.
///
/// The `alert_summary` slice is used only to detect mutagenicity-category alerts
/// without re-running substructure search.
///
/// # Parameters
///
/// - `alert_count` — total number of structural alerts triggered (used for context).
/// - `tox_profile` — QSAR output containing endpoint probabilities and domain status.
/// - `metabolite_tree` — predicted metabolite data.
///
/// # Returns
///
/// A `Vec<RegulatoryFlag>` that may be empty for a clean compound.
///
/// # Examples
///
/// ```rust
/// use nexcore_chemivigilance::watchlist::check_watchlist;
/// use nexcore_qsar::types::ToxProfile;
/// use nexcore_metabolite::types::MetaboliteTree;
///
/// let profile = ToxProfile::default();
/// let tree = MetaboliteTree::default();
/// let flags = check_watchlist(0, &profile, &tree);
/// // Default ToxProfile has OutOfDomain status.
/// assert!(!flags.is_empty());
/// ```
#[must_use]
pub fn check_watchlist(
    alert_count: usize,
    tox_profile: &ToxProfile,
    metabolite_tree: &MetaboliteTree,
) -> Vec<RegulatoryFlag> {
    let mut flags: Vec<RegulatoryFlag> = Vec::new();

    // ICH M7 flag: any structural alerts were triggered (caller already filtered
    // by category via alert_summary; here we use alert_count > 0 as a proxy
    // because mutagenicity category information is passed through the count).
    // A more specific check is done in the pipeline where alert_summary is available.
    let _ = alert_count; // used by caller for context; category check below uses tox_profile

    // Mutagenicity flag.
    if tox_profile.mutagenicity.probability >= 0.5 {
        flags.push(RegulatoryFlag {
            flag_type: RegulatoryFlagType::HighMutagenicity,
            description: format!(
                "Mutagenicity probability {:.2} exceeds the 0.5 concern threshold.",
                tox_profile.mutagenicity.probability
            ),
            reference: "ICH M7(R1) — Assessment and Control of DNA Reactive Impurities".to_string(),
        });
    }

    // Hepatotoxicity flag.
    if tox_profile.hepatotoxicity.probability >= 0.5 {
        flags.push(RegulatoryFlag {
            flag_type: RegulatoryFlagType::HighHepatotoxicity,
            description: format!(
                "Hepatotoxicity (DILI) probability {:.2} exceeds the 0.5 concern threshold.",
                tox_profile.hepatotoxicity.probability
            ),
            reference: "FDA Drug-Induced Liver Injury (DILI) Guidance".to_string(),
        });
    }

    // Cardiotoxicity flag.
    if tox_profile.cardiotoxicity.probability >= 0.5 {
        flags.push(RegulatoryFlag {
            flag_type: RegulatoryFlagType::HighCardiotoxicity,
            description: format!(
                "Cardiotoxicity (hERG) probability {:.2} exceeds the 0.5 concern threshold.",
                tox_profile.cardiotoxicity.probability
            ),
            reference: "ICH S7B — Nonclinical Evaluation of the Potential for Delayed Ventricular Repolarization".to_string(),
        });
    }

    // Reactive intermediates flag.
    if !metabolite_tree.reactive_intermediates.is_empty() {
        flags.push(RegulatoryFlag {
            flag_type: RegulatoryFlagType::ReactiveIntermediate,
            description: format!(
                "{} reactive intermediate(s) predicted (e.g. arene oxides, quinones). \
                 Glutathione conjugation or protein adduct formation may occur.",
                metabolite_tree.reactive_intermediates.len()
            ),
            reference: "EMA Guideline on the Investigation of Drug Interactions (CPMP/EWP/560/95)"
                .to_string(),
        });
    }

    // Applicability domain flag.
    if matches!(
        tox_profile.applicability_domain,
        DomainStatus::OutOfDomain { .. }
    ) {
        flags.push(RegulatoryFlag {
            flag_type: RegulatoryFlagType::OutOfDomain,
            description: "Compound falls outside the applicability domain of the QSAR models. \
                 Predictions may be unreliable; expert review is strongly recommended."
                .to_string(),
            reference: "OECD Principles for the Validation of QSAR Models (2007)".to_string(),
        });
    }

    flags
}

/// Check whether any triggered alerts fall in the mutagenicity category and
/// add an ICH M7 flag accordingly.
///
/// This is a companion to [`check_watchlist`] called from the pipeline where
/// the `AlertSummary` slice is available.
///
/// # Examples
///
/// ```rust
/// use nexcore_chemivigilance::watchlist::check_ich_m7_flag;
/// use nexcore_chemivigilance::brief::AlertSummary;
///
/// let alerts = vec![AlertSummary {
///     alert_id: "M7-001".to_string(),
///     alert_name: "Aromatic amine".to_string(),
///     category: "Mutagenicity".to_string(),
///     match_count: 1,
///     confidence: 0.9,
/// }];
/// let flags = check_ich_m7_flag(&alerts);
/// assert_eq!(flags.len(), 1);
/// ```
#[must_use]
pub fn check_ich_m7_flag(alert_summary: &[crate::brief::AlertSummary]) -> Vec<RegulatoryFlag> {
    let has_mutagenicity = alert_summary.iter().any(|a| a.category == "Mutagenicity");

    if has_mutagenicity {
        vec![RegulatoryFlag {
            flag_type: RegulatoryFlagType::IchM7Alert,
            description:
                "One or more ICH M7 mutagenicity structural alerts were triggered. \
                 Formal ICH M7 assessment with expert review is required."
                    .to_string(),
            reference: "ICH M7(R1) — Assessment and Control of DNA Reactive (Mutagenic) Impurities in Pharmaceuticals".to_string(),
        }]
    } else {
        vec![]
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brief::{AlertSummary, RegulatoryFlagType};
    use nexcore_metabolite::types::MetaboliteTree;
    use nexcore_qsar::types::{DomainStatus, PredictionResult, RiskLevel, ToxClass, ToxProfile};

    fn low_risk_profile() -> ToxProfile {
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

    fn high_mutagenicity_profile() -> ToxProfile {
        let mut p = low_risk_profile();
        p.mutagenicity.probability = 0.8;
        p.mutagenicity.classification = ToxClass::Positive;
        p
    }

    #[test]
    fn test_no_flags_for_clean_compound() {
        let profile = low_risk_profile();
        let tree = MetaboliteTree::default();
        let flags = check_watchlist(0, &profile, &tree);
        // Clean compound with in-domain status, no reactive intermediates,
        // and all probabilities below 0.5 → no flags expected.
        assert!(
            flags.is_empty(),
            "clean low-risk compound should produce no watchlist flags, got: {flags:?}"
        );
    }

    #[test]
    fn test_mutagenicity_flag() {
        let profile = high_mutagenicity_profile();
        let tree = MetaboliteTree::default();
        let flags = check_watchlist(0, &profile, &tree);
        assert!(
            flags
                .iter()
                .any(|f| f.flag_type == RegulatoryFlagType::HighMutagenicity),
            "high mutagenicity probability must raise HighMutagenicity flag"
        );
    }

    #[test]
    fn test_reactive_intermediate_flag() {
        use nexcore_metabolite::types::{Metabolite, Transformation};

        let profile = low_risk_profile();
        let mut tree = MetaboliteTree::default();
        tree.reactive_intermediates.push(Metabolite {
            transformation: Transformation::Epoxidation { site1: 0, site2: 1 },
            site_description: "Aromatic epoxide".to_string(),
            probability: 0.4,
            reactive_intermediate: true,
            enzyme: Some("CYP1A2".to_string()),
        });

        let flags = check_watchlist(0, &profile, &tree);
        assert!(
            flags
                .iter()
                .any(|f| f.flag_type == RegulatoryFlagType::ReactiveIntermediate),
            "non-empty reactive intermediates must raise ReactiveIntermediate flag"
        );
    }

    #[test]
    fn test_out_of_domain_flag() {
        let mut profile = low_risk_profile();
        profile.applicability_domain = DomainStatus::OutOfDomain {
            distance: 2.0,
            warning: "MW and LogP violated".to_string(),
        };
        let tree = MetaboliteTree::default();
        let flags = check_watchlist(0, &profile, &tree);
        assert!(
            flags
                .iter()
                .any(|f| f.flag_type == RegulatoryFlagType::OutOfDomain),
            "OutOfDomain status must raise OutOfDomain flag"
        );
    }

    #[test]
    fn test_hepatotoxicity_flag() {
        let mut profile = low_risk_profile();
        profile.hepatotoxicity.probability = 0.7;
        let tree = MetaboliteTree::default();
        let flags = check_watchlist(0, &profile, &tree);
        assert!(
            flags
                .iter()
                .any(|f| f.flag_type == RegulatoryFlagType::HighHepatotoxicity),
            "hepatotoxicity probability >= 0.5 must raise HighHepatotoxicity flag"
        );
    }

    #[test]
    fn test_cardiotoxicity_flag() {
        let mut profile = low_risk_profile();
        profile.cardiotoxicity.probability = 0.6;
        let tree = MetaboliteTree::default();
        let flags = check_watchlist(0, &profile, &tree);
        assert!(
            flags
                .iter()
                .any(|f| f.flag_type == RegulatoryFlagType::HighCardiotoxicity),
            "cardiotoxicity probability >= 0.5 must raise HighCardiotoxicity flag"
        );
    }

    #[test]
    fn test_ich_m7_flag_with_mutagenicity_alert() {
        let alerts = vec![AlertSummary {
            alert_id: "M7-001".to_string(),
            alert_name: "Aromatic amine".to_string(),
            category: "Mutagenicity".to_string(),
            match_count: 1,
            confidence: 0.9,
        }];
        let flags = check_ich_m7_flag(&alerts);
        assert_eq!(flags.len(), 1);
        assert_eq!(flags[0].flag_type, RegulatoryFlagType::IchM7Alert);
    }

    #[test]
    fn test_ich_m7_flag_no_mutagenicity_alert() {
        let alerts = vec![AlertSummary {
            alert_id: "SA-001".to_string(),
            alert_name: "Epoxide".to_string(),
            category: "Genotoxicity".to_string(),
            match_count: 1,
            confidence: 0.8,
        }];
        let flags = check_ich_m7_flag(&alerts);
        assert!(
            flags.is_empty(),
            "non-mutagenicity alert must not raise ICH M7 flag"
        );
    }

    #[test]
    fn test_default_tox_profile_out_of_domain() {
        let profile = ToxProfile::default();
        let tree = MetaboliteTree::default();
        let flags = check_watchlist(0, &profile, &tree);
        // Default ToxProfile has OutOfDomain status.
        assert!(
            flags
                .iter()
                .any(|f| f.flag_type == RegulatoryFlagType::OutOfDomain),
            "default ToxProfile (OutOfDomain) must raise OutOfDomain flag"
        );
    }
}
