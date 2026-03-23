// Copyright © 2026 NexVigilant LLC. All Rights Reserved.

//! Chemivigilance pipeline orchestrator generating SafetyBrief with min-3 limitations.
//!
//! # Architecture
//!
//! The pipeline accepts a SMILES string and returns a [`SafetyBrief`] containing:
//! - Molecular descriptors
//! - Structural alert scan results (ICH M7 library)
//! - QSAR toxicity predictions (mutagenicity, hepatotoxicity, cardiotoxicity)
//! - Predicted metabolites and degradants
//! - Regulatory flags
//! - At least 3 methodological limitations (ToV Axiom 4)
//!
//! # Quick start
//!
//! ```rust
//! use nexcore_chemivigilance::{generate_safety_brief, ChemivigilanceConfig};
//!
//! let config = ChemivigilanceConfig::default();
//! let brief = generate_safety_brief("CC(=O)Oc1ccccc1C(=O)O", &config)
//!     .unwrap_or_else(|e| panic!("pipeline failed: {e}"));
//! assert!(brief.limitations.len() >= 3);
//! ```

#![forbid(unsafe_code)]
#![cfg_attr(not(test), deny(clippy::unwrap_used))]
#![cfg_attr(not(test), deny(clippy::expect_used))]
#![cfg_attr(not(test), deny(clippy::panic))]
#![warn(missing_docs)]
pub mod brief;
pub mod error;
pub mod pipeline;
pub mod regulatory;
pub mod watchlist;

pub use brief::{
    AlertSummary, DescriptorSummary, Limitation, LimitationCategory, LimitationSeverity,
    RegulatoryFlag, RegulatoryFlagType, SafetyBrief,
};
pub use error::{ChemivigilanceError, ChemivigilanceResult};
pub use pipeline::{ChemivigilanceConfig, generate_safety_brief};
