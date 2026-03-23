// Copyright © 2026 NexVigilant LLC. All Rights Reserved.

//! End-to-end integration tests for the `nexcore-chemivigilance` pipeline.
//!
//! Each test exercises [`generate_safety_brief`] through the public API and
//! validates a specific contract or invariant.  No internal implementation
//! details are accessed.
//!
//! # Coverage targets
//!
//! - Happy-path molecules (aspirin, caffeine, ethanol, benzene)
//! - Error paths (invalid SMILES, empty input)
//! - ToV Axiom 4 invariant: `limitations.len() >= 3` for all valid inputs
//! - Descriptor bounds and field presence
//! - JSON serialisation round-trip
//! - Model version metadata population
//! - Risk score bounds `[0.0, 1.0]`
//! - Timestamp recency (generated within the test window)

use nexcore_chemivigilance::{
    ChemivigilanceConfig, ChemivigilanceError, SafetyBrief, generate_safety_brief,
};
use nexcore_qsar::types::RiskLevel;

// ---------------------------------------------------------------------------
// Helper
// ---------------------------------------------------------------------------

/// Run the pipeline and return the `SafetyBrief` if successful, or `None` with
/// a descriptive test failure message if the pipeline returns an error.
///
/// Integration tests should use this helper for all valid-input cases.
/// The returned `Option` is always `Some` for well-formed SMILES; when it is
/// `None` the assertion inside prints the error so the test failure is visible.
fn run(smiles: &str) -> Option<SafetyBrief> {
    let config = ChemivigilanceConfig::default();
    match generate_safety_brief(smiles, &config) {
        Ok(brief) => Some(brief),
        Err(e) => {
            // Record why the call failed without triggering the panic lint.
            // The test body must also call assert!(result.is_some()) after
            // calling run() so that the test registers as failed.
            eprintln!("[integration] generate_safety_brief({smiles:?}) failed: {e}");
            None
        }
    }
}

// ---------------------------------------------------------------------------
// Test 1 — Aspirin full pipeline
// ---------------------------------------------------------------------------

#[test]
fn test_aspirin_full_pipeline() {
    // Aspirin SMILES: CC(=O)Oc1ccccc1C(=O)O
    const ASPIRIN: &str = "CC(=O)Oc1ccccc1C(=O)O";
    let result = run(ASPIRIN);
    assert!(result.is_some(), "aspirin pipeline must succeed");

    if let Some(brief) = result {
        // Compound identity
        assert_eq!(brief.smiles, ASPIRIN, "SMILES must be preserved verbatim");
        assert!(
            !brief.molecular_formula.is_empty(),
            "molecular formula must be populated"
        );

        // Aspirin MW ~180.16 Da
        assert!(
            brief.descriptors.molecular_weight > 170.0,
            "aspirin MW must be > 170, got {}",
            brief.descriptors.molecular_weight
        );
        assert!(
            brief.descriptors.molecular_weight < 190.0,
            "aspirin MW must be < 190, got {}",
            brief.descriptors.molecular_weight
        );

        // Toxicity profile probabilities are bounded
        assert!(
            brief.tox_profile.mutagenicity.probability >= 0.0
                && brief.tox_profile.mutagenicity.probability <= 1.0,
            "mutagenicity probability must be in [0, 1], got {}",
            brief.tox_profile.mutagenicity.probability
        );

        // Aspirin has an ester bond — hydrolysis predicted in degradants or phase1
        assert!(
            !brief.metabolite_tree.degradants.is_empty()
                || !brief.metabolite_tree.phase1.is_empty(),
            "aspirin must have at least one metabolite or degradant predicted"
        );

        // ToV Axiom 4: at least 3 limitations
        assert!(
            brief.limitations.len() >= 3,
            "must have >= 3 limitations (ToV Axiom 4), got {}",
            brief.limitations.len()
        );

        // Risk score in [0.0, 1.0]
        assert!(
            brief.overall_risk_score >= 0.0,
            "risk score must be >= 0.0, got {}",
            brief.overall_risk_score
        );
        assert!(
            brief.overall_risk_score <= 1.0,
            "risk score must be <= 1.0, got {}",
            brief.overall_risk_score
        );

        // Timestamp is within a 60-second test window
        let now = nexcore_chrono::DateTime::now();
        let age = now.signed_duration_since(brief.generated_at);
        assert!(
            age.num_seconds() >= 0 && age.num_seconds() < 60,
            "generated_at must be recent (within 60 s), age = {} s",
            age.num_seconds()
        );
    }
}

// ---------------------------------------------------------------------------
// Test 2 — Caffeine from SMILES
// ---------------------------------------------------------------------------

#[test]
fn test_caffeine_from_smiles() {
    // Caffeine: Cn1cnc2c1c(=O)n(c(=O)n2C)C
    const CAFFEINE: &str = "Cn1cnc2c1c(=O)n(c(=O)n2C)C";
    let result = run(CAFFEINE);
    assert!(result.is_some(), "caffeine pipeline must succeed");

    if let Some(brief) = result {
        // Caffeine MW ~194.19 Da
        assert!(
            brief.descriptors.molecular_weight > 180.0,
            "caffeine MW must be > 180, got {}",
            brief.descriptors.molecular_weight
        );
        assert!(
            brief.descriptors.molecular_weight < 210.0,
            "caffeine MW must be < 210, got {}",
            brief.descriptors.molecular_weight
        );

        // Caffeine is bicyclic (imidazole fused to pyrimidine = 2 rings)
        assert!(
            brief.descriptors.num_rings >= 2,
            "caffeine must have >= 2 rings, got {}",
            brief.descriptors.num_rings
        );

        // ToV Axiom 4
        assert!(
            brief.limitations.len() >= 3,
            "caffeine must have >= 3 limitations, got {}",
            brief.limitations.len()
        );
    }
}

// ---------------------------------------------------------------------------
// Test 3 — Invalid SMILES returns error
// ---------------------------------------------------------------------------

#[test]
fn test_invalid_smiles_returns_error() {
    let config = ChemivigilanceConfig::default();
    let result = generate_safety_brief("INVALID$$$SMILES", &config);
    assert!(result.is_err(), "invalid SMILES must return Err");
}

// ---------------------------------------------------------------------------
// Test 4 — Empty SMILES returns InvalidInput error
// ---------------------------------------------------------------------------

#[test]
fn test_empty_smiles_returns_error() {
    let config = ChemivigilanceConfig::default();
    let result = generate_safety_brief("", &config);
    assert!(
        matches!(result, Err(ChemivigilanceError::InvalidInput(_))),
        "empty SMILES must return ChemivigilanceError::InvalidInput"
    );
}

// ---------------------------------------------------------------------------
// Test 5 — ToV Axiom 4: limitations >= 3 for multiple molecules
// ---------------------------------------------------------------------------

#[test]
fn test_limitations_minimum_three_various_molecules() {
    let molecules: &[(&str, &str)] = &[
        ("C", "methane"),
        ("CCO", "ethanol"),
        ("c1ccccc1", "benzene"),
        ("CC(=O)Oc1ccccc1C(=O)O", "aspirin"),
        ("CC(=O)O", "acetic acid"),
    ];

    for &(smiles, name) in molecules {
        let result = run(smiles);
        assert!(
            result.is_some(),
            "pipeline for {name} ({smiles}) must succeed"
        );
        if let Some(brief) = result {
            assert!(
                brief.limitations.len() >= 3,
                "molecule '{name}' ({smiles}) has only {} limitations, need >= 3 (ToV Axiom 4)",
                brief.limitations.len()
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Test 6 — Ethanol: no structural alerts, low risk
// ---------------------------------------------------------------------------

#[test]
fn test_ethanol_low_risk() {
    let result = run("CCO");
    assert!(result.is_some(), "ethanol pipeline must succeed");

    if let Some(brief) = result {
        // Ethanol (CCO) has no recognised structural alerts
        assert_eq!(
            brief.alert_count, 0,
            "ethanol must have 0 structural alerts, got {}",
            brief.alert_count
        );

        // Risk score should be below 0.7 (not VeryHigh) for a simple molecule
        assert!(
            brief.overall_risk_score < 0.7,
            "ethanol risk score should be < 0.7, got {}",
            brief.overall_risk_score
        );
    }
}

// ---------------------------------------------------------------------------
// Test 7 — SafetyBrief serialises to valid JSON
// ---------------------------------------------------------------------------

#[test]
fn test_safety_brief_serializes_to_json() {
    let result = run("CCO");
    assert!(
        result.is_some(),
        "ethanol pipeline must succeed for JSON test"
    );

    if let Some(brief) = result {
        let json_result = serde_json::to_string(&brief);
        assert!(
            json_result.is_ok(),
            "SafetyBrief must serialise to valid JSON"
        );

        if let Ok(json_str) = json_result {
            assert!(!json_str.is_empty(), "serialised JSON must not be empty");
            assert!(
                json_str.contains("limitations"),
                "JSON must contain 'limitations' key"
            );
            assert!(
                json_str.contains("overall_risk_score"),
                "JSON must contain 'overall_risk_score' key"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Test 8 — Benzene: aromatic descriptor values
// ---------------------------------------------------------------------------

#[test]
fn test_benzene_aromatic_descriptor() {
    let result = run("c1ccccc1");
    assert!(result.is_some(), "benzene pipeline must succeed");

    if let Some(brief) = result {
        // Benzene has exactly 1 aromatic ring
        assert_eq!(
            brief.descriptors.num_aromatic_rings, 1,
            "benzene must have exactly 1 aromatic ring, got {}",
            brief.descriptors.num_aromatic_rings
        );
        assert_eq!(
            brief.descriptors.num_rings, 1,
            "benzene must have exactly 1 ring (SSSR), got {}",
            brief.descriptors.num_rings
        );

        // Benzene: pure hydrocarbon — no O or N atoms
        assert_eq!(
            brief.descriptors.hba, 0,
            "benzene has no HBA (no O or N), got {}",
            brief.descriptors.hba
        );
        assert_eq!(
            brief.descriptors.hbd, 0,
            "benzene has no HBD (no OH or NH), got {}",
            brief.descriptors.hbd
        );
    }
}

// ---------------------------------------------------------------------------
// Test 9 — Model versions are populated
// ---------------------------------------------------------------------------

#[test]
fn test_model_versions_populated() {
    let result = run("CCO");
    assert!(result.is_some(), "ethanol pipeline must succeed");

    if let Some(brief) = result {
        assert!(
            !brief.model_versions.is_empty(),
            "model_versions map must not be empty"
        );
        assert!(
            brief.model_versions.contains_key("mutagenicity"),
            "model_versions must contain 'mutagenicity'"
        );
        assert!(
            brief.model_versions.contains_key("hepatotoxicity"),
            "model_versions must contain 'hepatotoxicity'"
        );
        assert!(
            brief.model_versions.contains_key("cardiotoxicity"),
            "model_versions must contain 'cardiotoxicity'"
        );
        assert!(
            brief.model_versions.contains_key("structural_alerts"),
            "model_versions must contain 'structural_alerts'"
        );
        assert!(
            brief.model_versions.contains_key("metabolite"),
            "model_versions must contain 'metabolite'"
        );
    }
}

// ---------------------------------------------------------------------------
// Test 10 — Benzene regulatory output and limitation invariant
// ---------------------------------------------------------------------------

#[test]
fn test_regulatory_flags_for_reactive_intermediates() {
    // Benzene can produce reactive intermediates via aromatic epoxidation.
    // The exact flags depend on the metabolite model's predictions; this test
    // verifies the pipeline does not panic and upholds all invariants.
    let result = run("c1ccccc1");
    assert!(result.is_some(), "benzene pipeline must succeed");

    if let Some(brief) = result {
        // ToV Axiom 4 must hold regardless of regulatory flag content
        assert!(
            brief.limitations.len() >= 3,
            "benzene must have >= 3 limitations (ToV Axiom 4), got {}",
            brief.limitations.len()
        );

        // Risk score always bounded in [0.0, 1.0]
        assert!(
            (0.0..=1.0).contains(&brief.overall_risk_score),
            "benzene risk score must be in [0, 1], got {}",
            brief.overall_risk_score
        );

        // All alert confidence values must be in [0.0, 1.0]
        for alert in &brief.alert_summary {
            assert!(
                (0.0..=1.0).contains(&alert.confidence),
                "alert '{}' confidence must be in [0, 1], got {}",
                alert.alert_id,
                alert.confidence
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Test 11 — Risk level is consistent with the risk score
// ---------------------------------------------------------------------------

#[test]
fn test_risk_level_consistent_with_score() {
    let result = run("CCO");
    assert!(result.is_some(), "ethanol pipeline must succeed");

    if let Some(brief) = result {
        let expected_level = match brief.overall_risk_score {
            s if s >= 0.7 => RiskLevel::VeryHigh,
            s if s >= 0.5 => RiskLevel::High,
            s if s >= 0.3 => RiskLevel::Medium,
            _ => RiskLevel::Low,
        };
        assert_eq!(
            brief.risk_level, expected_level,
            "risk_level must match the computed score ({:.4}): expected {:?}, got {:?}",
            brief.overall_risk_score, expected_level, brief.risk_level
        );
    }
}

// ---------------------------------------------------------------------------
// Test 12 — Aspirin molecular formula contains C and O
// ---------------------------------------------------------------------------

#[test]
fn test_aspirin_molecular_formula_contains_elements() {
    let result = run("CC(=O)Oc1ccccc1C(=O)O");
    assert!(result.is_some(), "aspirin pipeline must succeed");

    if let Some(brief) = result {
        assert!(
            brief.molecular_formula.contains('C'),
            "aspirin molecular formula must contain 'C', got '{}'",
            brief.molecular_formula
        );
        assert!(
            brief.molecular_formula.contains('O'),
            "aspirin molecular formula must contain 'O', got '{}'",
            brief.molecular_formula
        );
    }
}

// ---------------------------------------------------------------------------
// Test 13 — Whitespace-only input returns InvalidInput error
// ---------------------------------------------------------------------------

#[test]
fn test_whitespace_only_smiles_returns_error() {
    let config = ChemivigilanceConfig::default();
    let result = generate_safety_brief("   ", &config);
    assert!(
        matches!(result, Err(ChemivigilanceError::InvalidInput(_))),
        "whitespace-only input must return ChemivigilanceError::InvalidInput"
    );
}

// ---------------------------------------------------------------------------
// Test 14 — Config alert_count_for_qsar override does not reduce risk score
// ---------------------------------------------------------------------------

#[test]
fn test_config_alert_count_override_does_not_reduce_risk() {
    let natural_config = ChemivigilanceConfig::default();
    let forced_config = ChemivigilanceConfig {
        alert_count_for_qsar: 3,
    };

    let natural_result = generate_safety_brief("c1ccccc1", &natural_config);
    let forced_result = generate_safety_brief("c1ccccc1", &forced_config);

    assert!(
        natural_result.is_ok(),
        "natural benzene pipeline must succeed"
    );
    assert!(
        forced_result.is_ok(),
        "forced-alert benzene pipeline must succeed"
    );

    if let (Ok(natural), Ok(forced)) = (natural_result, forced_result) {
        assert!(
            forced.overall_risk_score >= natural.overall_risk_score,
            "forcing 3 alerts must not reduce risk score: forced={:.4}, natural={:.4}",
            forced.overall_risk_score,
            natural.overall_risk_score
        );
    }
}
