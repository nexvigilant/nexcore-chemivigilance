// Copyright © 2026 NexVigilant LLC. All Rights Reserved.

//! Main chemivigilance pipeline orchestrator.
//!
//! [`generate_safety_brief`] is the single entry point that drives the full
//! six-stage pipeline from a SMILES string to a completed [`SafetyBrief`].
//!
//! ## Pipeline stages
//!
//! ```text
//! SMILES
//!   1. parse(smiles) → MolGraph
//!   2. calculate_descriptors(&graph) → Descriptors
//!   3. AlertLibrary::default_library() + scan_smiles(smiles) → Vec<AlertMatch>
//!   4. predict_from_descriptors(&desc, alert_count, 0) → ToxProfile
//!   5. predict_metabolites(&graph, smiles) → MetaboliteTree
//!   6. check_watchlist + check_ich_m7_flag → Vec<RegulatoryFlag>
//!   7. generate_limitations → Vec<Limitation>  (>= 3, ToV Axiom 4)
//!   8. compute risk score + RiskLevel
//!   9. assemble SafetyBrief
//! ```
//!
//! # Examples
//!
//! ```rust
//! use nexcore_chemivigilance::pipeline::{generate_safety_brief, ChemivigilanceConfig};
//!
//! let config = ChemivigilanceConfig::default();
//! let brief = generate_safety_brief("CCO", &config).unwrap_or_else(|_| {
//!     panic!("ethanol pipeline failed")
//! });
//! assert!(brief.limitations.len() >= 3);
//! ```

use std::collections::HashMap;

use nexcore_chrono::DateTime;
use nexcore_metabolite::predict::predict_metabolites;
use nexcore_molcore::descriptor::calculate_descriptors;
use nexcore_molcore::graph::MolGraph;
use nexcore_molcore::smiles::parse;
use nexcore_qsar::predict::predict_from_descriptors;
use nexcore_qsar::types::{DomainStatus, RiskLevel};
use nexcore_structural_alerts::{AlertCategory, AlertLibrary, scan_smiles};

use crate::brief::{AlertSummary, DescriptorSummary, SafetyBrief};
use crate::error::{ChemivigilanceError, ChemivigilanceResult};
use crate::regulatory::generate_limitations;
use crate::watchlist::{check_ich_m7_flag, check_watchlist};

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

/// Configuration for the chemivigilance pipeline.
///
/// # Examples
///
/// ```rust
/// use nexcore_chemivigilance::pipeline::ChemivigilanceConfig;
///
/// let config = ChemivigilanceConfig::default();
/// assert_eq!(config.alert_count_for_qsar, 0);
/// ```
#[derive(Debug, Clone, Default)]
pub struct ChemivigilanceConfig {
    /// Override the structural alert count forwarded to QSAR models.
    ///
    /// When `0` (the default) the actual count from the structural alert scan
    /// is used.  Set to a non-zero value to simulate a specific alert burden
    /// during testing.
    pub alert_count_for_qsar: usize,
}

// ---------------------------------------------------------------------------
// Pipeline
// ---------------------------------------------------------------------------

/// Run the full chemivigilance pipeline for a SMILES input.
///
/// Returns a [`SafetyBrief`] on success, or a [`ChemivigilanceError`] if any
/// upstream step fails (e.g. invalid SMILES).
///
/// # Errors
///
/// - [`ChemivigilanceError::InvalidInput`] — `smiles` is empty.
/// - [`ChemivigilanceError::MolcoreError`] — SMILES parse failure.
/// - [`ChemivigilanceError::AlertError`] — structural alert library error.
/// - [`ChemivigilanceError::InsufficientLimitations`] — internal invariant
///   violation (indicates a bug; should never occur in practice).
///
/// # Examples
///
/// ```rust
/// use nexcore_chemivigilance::pipeline::{generate_safety_brief, ChemivigilanceConfig};
///
/// let brief = generate_safety_brief("CC(=O)Oc1ccccc1C(=O)O", &ChemivigilanceConfig::default())
///     .unwrap_or_else(|e| panic!("aspirin pipeline failed: {e}"));
/// assert!(brief.limitations.len() >= 3);
/// assert!(!brief.smiles.is_empty());
/// ```
pub fn generate_safety_brief(
    smiles: &str,
    config: &ChemivigilanceConfig,
) -> ChemivigilanceResult<SafetyBrief> {
    // ------------------------------------------------------------------
    // Guard: reject empty input immediately.
    // ------------------------------------------------------------------
    if smiles.trim().is_empty() {
        return Err(ChemivigilanceError::InvalidInput(
            "SMILES string must not be empty".to_string(),
        ));
    }

    // ------------------------------------------------------------------
    // Stage 1: Parse SMILES → MolGraph
    // ------------------------------------------------------------------
    let molecule = parse(smiles).map_err(|e| ChemivigilanceError::MolcoreError(e.to_string()))?;
    let graph = MolGraph::from_molecule(molecule);

    // ------------------------------------------------------------------
    // Stage 2: Calculate molecular descriptors
    // ------------------------------------------------------------------
    let descriptors = calculate_descriptors(&graph);

    // ------------------------------------------------------------------
    // Stage 3: Structural alert scan
    // ------------------------------------------------------------------
    let alert_library = AlertLibrary::default_library();
    let alert_matches = scan_smiles(smiles, &alert_library).map_err(ChemivigilanceError::from)?;

    let alert_summary: Vec<AlertSummary> = alert_matches
        .iter()
        .map(|m| AlertSummary {
            alert_id: m.alert.id.clone(),
            alert_name: m.alert.name.clone(),
            category: format!("{:?}", m.alert.category),
            match_count: m.match_count,
            confidence: m.alert.confidence,
        })
        .collect();

    let alert_count = alert_matches.len();

    // Effective alert count for QSAR — use config override if non-zero.
    let qsar_alert_count = if config.alert_count_for_qsar > 0 {
        config.alert_count_for_qsar
    } else {
        alert_count
    };

    // ------------------------------------------------------------------
    // Stage 4: QSAR toxicity prediction
    // ------------------------------------------------------------------
    let tox_profile = predict_from_descriptors(&descriptors, qsar_alert_count, 0);

    // ------------------------------------------------------------------
    // Stage 5: Metabolite prediction
    // ------------------------------------------------------------------
    let metabolite_tree = predict_metabolites(&graph, smiles);

    // ------------------------------------------------------------------
    // Stage 6: Watchlist / regulatory flags
    // ------------------------------------------------------------------
    let mut regulatory_flags = check_watchlist(alert_count, &tox_profile, &metabolite_tree);
    let ich_m7_flags = check_ich_m7_flag(&alert_summary);
    regulatory_flags.extend(ich_m7_flags);

    // ------------------------------------------------------------------
    // Stage 7: Generate limitations (min 3, ToV Axiom 4)
    // ------------------------------------------------------------------
    let limitations = generate_limitations(&descriptors, &tox_profile, alert_count);

    if limitations.len() < 3 {
        return Err(ChemivigilanceError::InsufficientLimitations(
            limitations.len(),
        ));
    }

    // ------------------------------------------------------------------
    // Stage 8: Compute overall risk score and level
    // ------------------------------------------------------------------
    let base_score = tox_profile
        .mutagenicity
        .probability
        .max(tox_profile.hepatotoxicity.probability)
        .max(tox_profile.cardiotoxicity.probability);

    let alert_boost = if alert_count > 0 {
        (0.1_f64 * alert_count as f64).min(0.2)
    } else {
        0.0
    };

    let reactive_boost = if metabolite_tree.reactive_intermediates.is_empty() {
        0.0
    } else {
        0.1
    };

    let overall_risk_score = (base_score + alert_boost + reactive_boost).min(1.0_f64);

    let risk_level = score_to_risk_level(overall_risk_score);

    // ------------------------------------------------------------------
    // Stage 9: Molecular formula (lightweight — count atoms by element)
    // ------------------------------------------------------------------
    let molecular_formula = compute_molecular_formula(smiles);

    // ------------------------------------------------------------------
    // Stage 10: Model version metadata
    // ------------------------------------------------------------------
    let mut model_versions = HashMap::new();
    model_versions.insert(
        "mutagenicity".to_string(),
        tox_profile.mutagenicity.model_version.clone(),
    );
    model_versions.insert(
        "hepatotoxicity".to_string(),
        tox_profile.hepatotoxicity.model_version.clone(),
    );
    model_versions.insert(
        "cardiotoxicity".to_string(),
        tox_profile.cardiotoxicity.model_version.clone(),
    );
    model_versions.insert(
        "structural_alerts".to_string(),
        "nexcore-structural-alerts-v0.1.0".to_string(),
    );
    model_versions.insert(
        "metabolite".to_string(),
        "nexcore-metabolite-v0.1.0".to_string(),
    );

    // ------------------------------------------------------------------
    // Assemble SafetyBrief
    // ------------------------------------------------------------------
    let brief = SafetyBrief {
        input: smiles.to_string(),
        smiles: smiles.to_string(),
        molecular_formula,
        descriptors: DescriptorSummary::from(&descriptors),
        alert_count,
        alert_summary,
        tox_profile,
        metabolite_tree,
        overall_risk_score,
        risk_level,
        limitations,
        regulatory_flags,
        generated_at: DateTime::now(),
        model_versions,
    };

    Ok(brief)
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Derive a [`RiskLevel`] from a composite risk score.
///
/// Thresholds mirror those used inside `nexcore-qsar` for consistency:
///
/// ```text
/// score < 0.3  → Low
/// score < 0.5  → Medium
/// score < 0.7  → High
/// score >= 0.7 → VeryHigh
/// ```
fn score_to_risk_level(score: f64) -> RiskLevel {
    if score >= 0.7 {
        RiskLevel::VeryHigh
    } else if score >= 0.5 {
        RiskLevel::High
    } else if score >= 0.3 {
        RiskLevel::Medium
    } else {
        RiskLevel::Low
    }
}

/// Derive a simple molecular formula string from a SMILES by scanning the
/// parsed graph's atom list.
///
/// Returns `"unknown"` if parsing fails (the caller already checked validity
/// in stage 1, so this path is only reached in edge cases).
fn compute_molecular_formula(smiles: &str) -> String {
    let mol = match parse(smiles) {
        Ok(m) => m,
        Err(_) => return "unknown".to_string(),
    };

    // Collect element counts.
    let mut counts: std::collections::BTreeMap<u8, usize> = std::collections::BTreeMap::new();
    for atom in &mol.atoms {
        *counts.entry(atom.atomic_number).or_insert(0) += 1;
        // Add implicit hydrogens.
        *counts.entry(1).or_insert(0) += usize::from(atom.implicit_h);
    }

    // Build Hill notation: C first, H second, then alphabetical.
    let element_symbol = |an: u8| -> &'static str {
        match an {
            1 => "H",
            6 => "C",
            7 => "N",
            8 => "O",
            9 => "F",
            15 => "P",
            16 => "S",
            17 => "Cl",
            35 => "Br",
            53 => "I",
            _ => "?",
        }
    };

    let mut formula = String::new();

    // Hill order: C, H, then alphabetical remainder.
    let hill_order: &[u8] = &[6, 1, 7, 8, 9, 15, 16, 17, 35, 53];

    for &an in hill_order {
        if let Some(&count) = counts.get(&an) {
            formula.push_str(element_symbol(an));
            if count > 1 {
                formula.push_str(&count.to_string());
            }
        }
    }

    if formula.is_empty() {
        "unknown".to_string()
    } else {
        formula
    }
}

/// Indicates whether a [`DomainStatus`] represents an out-of-domain or
/// borderline assessment.
///
/// Used in tests to avoid matching on private enum variants across crate
/// boundaries.
#[must_use]
pub fn is_domain_flagged(status: &DomainStatus) -> bool {
    matches!(
        status,
        DomainStatus::OutOfDomain { .. } | DomainStatus::Borderline { .. }
    )
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use nexcore_qsar::types::RiskLevel;

    fn default_config() -> ChemivigilanceConfig {
        ChemivigilanceConfig::default()
    }

    // ------------------------------------------------------------------
    // Happy-path integration tests
    // ------------------------------------------------------------------

    #[test]
    fn test_aspirin_safety_brief() {
        let brief = generate_safety_brief("CC(=O)Oc1ccccc1C(=O)O", &default_config())
            .unwrap_or_else(|e| panic!("aspirin pipeline failed: {e}"));

        // Descriptors are populated.
        assert!(brief.descriptors.heavy_atom_count >= 13);
        assert!(brief.descriptors.molecular_weight > 100.0);
        assert_eq!(brief.descriptors.num_aromatic_rings, 1);

        // At least 3 limitations (ToV Axiom 4).
        assert!(
            brief.limitations.len() >= 3,
            "expected >= 3 limitations, got {}",
            brief.limitations.len()
        );

        // Risk level is set.
        let _ = &brief.risk_level; // just verifying it exists and is constructed

        // SMILES is preserved.
        assert_eq!(brief.smiles, "CC(=O)Oc1ccccc1C(=O)O");

        // Model versions populated.
        assert!(brief.model_versions.contains_key("mutagenicity"));
    }

    #[test]
    fn test_ethanol_safety_brief() {
        let brief = generate_safety_brief("CCO", &default_config())
            .unwrap_or_else(|e| panic!("ethanol pipeline failed: {e}"));

        // Ethanol is a small, simple molecule — expect Low or Medium risk.
        assert!(
            brief.risk_level == RiskLevel::Low || brief.risk_level == RiskLevel::Medium,
            "ethanol should be Low or Medium risk, got {:?}",
            brief.risk_level
        );

        assert!(brief.limitations.len() >= 3);
        assert_eq!(brief.smiles, "CCO");
    }

    #[test]
    fn test_benzene_safety_brief() {
        let brief = generate_safety_brief("c1ccccc1", &default_config())
            .unwrap_or_else(|e| panic!("benzene pipeline failed: {e}"));

        // Benzene has 1 aromatic ring.
        assert_eq!(brief.descriptors.num_aromatic_rings, 1);
        assert!(brief.limitations.len() >= 3);
    }

    // ------------------------------------------------------------------
    // Error path
    // ------------------------------------------------------------------

    #[test]
    fn test_invalid_smiles_returns_error() {
        let result = generate_safety_brief("INVALID$$SMILES", &default_config());
        assert!(result.is_err(), "invalid SMILES must return Err");
    }

    #[test]
    fn test_empty_smiles_returns_invalid_input_error() {
        let result = generate_safety_brief("", &default_config());
        assert!(
            matches!(result, Err(ChemivigilanceError::InvalidInput(_))),
            "empty SMILES must return InvalidInput error"
        );
    }

    #[test]
    fn test_whitespace_only_smiles_returns_invalid_input_error() {
        let result = generate_safety_brief("   ", &default_config());
        assert!(
            matches!(result, Err(ChemivigilanceError::InvalidInput(_))),
            "whitespace-only input must return InvalidInput error"
        );
    }

    // ------------------------------------------------------------------
    // ToV Axiom 4 invariant
    // ------------------------------------------------------------------

    #[test]
    fn test_minimum_three_limitations_aspirin() {
        let brief = generate_safety_brief("CC(=O)Oc1ccccc1C(=O)O", &default_config())
            .unwrap_or_else(|e| panic!("failed: {e}"));
        assert!(
            brief.limitations.len() >= 3,
            "ToV Axiom 4: min 3 limitations required"
        );
    }

    #[test]
    fn test_minimum_three_limitations_ethanol() {
        let brief = generate_safety_brief("CCO", &default_config())
            .unwrap_or_else(|e| panic!("failed: {e}"));
        assert!(
            brief.limitations.len() >= 3,
            "ToV Axiom 4: min 3 limitations required"
        );
    }

    #[test]
    fn test_minimum_three_limitations_benzene() {
        let brief = generate_safety_brief("c1ccccc1", &default_config())
            .unwrap_or_else(|e| panic!("failed: {e}"));
        assert!(
            brief.limitations.len() >= 3,
            "ToV Axiom 4: min 3 limitations required"
        );
    }

    // ------------------------------------------------------------------
    // Risk score bounds
    // ------------------------------------------------------------------

    #[test]
    fn test_risk_score_bounded() {
        let brief = generate_safety_brief("c1ccccc1[N+](=O)[O-]", &default_config())
            .unwrap_or_else(|e| panic!("nitrobenzene pipeline failed: {e}"));
        assert!(
            (0.0..=1.0).contains(&brief.overall_risk_score),
            "risk score must be in [0, 1], got {}",
            brief.overall_risk_score
        );
    }

    #[test]
    fn test_risk_level_consistent_with_score() {
        let brief = generate_safety_brief("CCO", &default_config())
            .unwrap_or_else(|e| panic!("failed: {e}"));

        let expected = score_to_risk_level(brief.overall_risk_score);
        assert_eq!(
            brief.risk_level, expected,
            "risk_level must be consistent with overall_risk_score"
        );
    }

    // ------------------------------------------------------------------
    // Molecular formula
    // ------------------------------------------------------------------

    #[test]
    fn test_aspirin_molecular_formula_contains_c_and_o() {
        let brief = generate_safety_brief("CC(=O)Oc1ccccc1C(=O)O", &default_config())
            .unwrap_or_else(|e| panic!("failed: {e}"));
        assert!(
            brief.molecular_formula.contains('C'),
            "aspirin formula must contain C, got {}",
            brief.molecular_formula
        );
        assert!(
            brief.molecular_formula.contains('O'),
            "aspirin formula must contain O, got {}",
            brief.molecular_formula
        );
    }

    // ------------------------------------------------------------------
    // Config override
    // ------------------------------------------------------------------

    #[test]
    fn test_config_alert_count_override() {
        // A config that forces 3 structural alerts into the QSAR model.
        let config = ChemivigilanceConfig {
            alert_count_for_qsar: 3,
        };
        let brief_forced =
            generate_safety_brief("c1ccccc1", &config).unwrap_or_else(|e| panic!("failed: {e}"));
        let brief_natural = generate_safety_brief("c1ccccc1", &default_config())
            .unwrap_or_else(|e| panic!("failed: {e}"));

        // With 3 forced alerts, risk should be >= natural (which may have 0).
        assert!(
            brief_forced.overall_risk_score >= brief_natural.overall_risk_score,
            "forced alerts must not reduce risk score"
        );
    }

    // ------------------------------------------------------------------
    // Score-to-risk helper
    // ------------------------------------------------------------------

    #[test]
    fn test_score_to_risk_level_boundaries() {
        assert_eq!(score_to_risk_level(0.0), RiskLevel::Low);
        assert_eq!(score_to_risk_level(0.29), RiskLevel::Low);
        assert_eq!(score_to_risk_level(0.3), RiskLevel::Medium);
        assert_eq!(score_to_risk_level(0.49), RiskLevel::Medium);
        assert_eq!(score_to_risk_level(0.5), RiskLevel::High);
        assert_eq!(score_to_risk_level(0.69), RiskLevel::High);
        assert_eq!(score_to_risk_level(0.7), RiskLevel::VeryHigh);
        assert_eq!(score_to_risk_level(1.0), RiskLevel::VeryHigh);
    }

    // ------------------------------------------------------------------
    // Metadata
    // ------------------------------------------------------------------

    #[test]
    fn test_model_versions_populated() {
        let brief = generate_safety_brief("CCO", &default_config())
            .unwrap_or_else(|e| panic!("failed: {e}"));
        assert!(brief.model_versions.contains_key("mutagenicity"));
        assert!(brief.model_versions.contains_key("hepatotoxicity"));
        assert!(brief.model_versions.contains_key("cardiotoxicity"));
        assert!(brief.model_versions.contains_key("structural_alerts"));
        assert!(brief.model_versions.contains_key("metabolite"));
    }

    #[test]
    fn test_generated_at_is_recent() {
        let before = DateTime::now();
        let brief = generate_safety_brief("CCO", &default_config())
            .unwrap_or_else(|e| panic!("failed: {e}"));
        let after = DateTime::now();
        assert!(
            brief.generated_at >= before && brief.generated_at <= after,
            "generated_at must be within the test window"
        );
    }
}
