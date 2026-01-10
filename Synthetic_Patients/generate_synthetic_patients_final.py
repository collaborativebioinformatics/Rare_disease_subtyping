#!/usr/bin/env python3
"""
================================================================================
ALS FEDERATED LEARNING SYNTHETIC PATIENT GENERATOR
Contextual Penetrance Model with Population-Specific Interactions
================================================================================

This pipeline generates synthetic ALS patient cohorts by integrating:
1. ClinVar pathogenic variant annotations
2. Estimated gnomAD population-specific allele frequencies
3. Hardy-Weinberg Equilibrium genotype simulation
4. Gene-based clinical phenotype assignment
5. Contextual Interaction Model (Ancestral Modifier System)

The key innovation is the CONTEXTUAL PENETRANCE MODEL, which uses population-
specific allele frequency ratios to modulate disease severity, simulating
real-world gene-environment interactions.

Output Files:
    - variants_complete.csv           : Variants with AF, phenotypes, severity
    - patients_complete.csv           : All patients with full annotations
    - patients_carriers_only.csv      : Only patients carrying ≥1 variant
    - client_<POP>.csv                : Federated client datasets
    - client_<POP>_carriers.csv       : Carrier-only federated datasets
    - validation_report.txt           : Statistical validation metrics
    - contextual_model_examples.txt   : Worked examples of the interaction model

Usage:
    python generate_als_synthetic_complete.py [path/to/clinvar_cleaned.csv]

Author: RAIDers ALS Subtyping Project
Version: 2.0 (with Contextual Interaction Model)
"""

import pandas as pd
import numpy as np
import os
import sys
import warnings
from datetime import datetime
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

# Number of patients per population (Total = N × 5)
N_PATIENTS_PER_POP = 3000  # 3000 × 5 = 15,000 total

# Minimum carriers per variant (ensures all variants represented)
MIN_CARRIERS_PER_VARIANT = 1

# Populations (federated clients representing global biobanks)
POPULATIONS = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']

# Output directory
OUTPUT_DIR = os.path.expanduser("~/Downloads/ALS_Synthetic_Data")

# Random seed for reproducibility
RANDOM_SEED = 42

# ============================================================================
# SECTION 1: GENE-TO-PHENOTYPE KNOWLEDGE BASE
# Source: OMIM, GeneReviews, ALS Literature
# ============================================================================

GENE_PHENOTYPE_MAP = {
    # === Classical ALS Types (Adult Onset) ===
    'SOD1': {
        'condition': 'Amyotrophic lateral sclerosis type 1 (ALS1)',
        'als_subtype': 'ALS1_Classical',
        'onset_type': 'Adult',
        'onset_age': '40-60 years',
        'base_progression': 'Variable',
        'base_severity': 7,
        'base_impact': 0.80,
        'notes': 'First ALS gene discovered; 20% of familial ALS; A4V aggressive, D90A slow'
    },
    'TARDBP': {
        'condition': 'Amyotrophic lateral sclerosis type 10 (ALS10)',
        'als_subtype': 'ALS10_Classical',
        'onset_type': 'Adult',
        'onset_age': '40-60 years',
        'base_progression': 'Moderate',
        'base_severity': 7,
        'base_impact': 0.75,
        'notes': 'TDP-43 pathology; may present with FTD features'
    },
    'FUS': {
        'condition': 'Amyotrophic lateral sclerosis type 6 (ALS6)',
        'als_subtype': 'ALS6_Aggressive',
        'onset_type': 'Variable',
        'onset_age': '20-50 years',
        'base_progression': 'Rapid',
        'base_severity': 9,
        'base_impact': 0.90,
        'notes': 'Often aggressive course; juvenile forms very severe'
    },
    'OPTN': {
        'condition': 'Amyotrophic lateral sclerosis type 12 (ALS12)',
        'als_subtype': 'ALS12_Classical',
        'onset_type': 'Adult',
        'onset_age': '30-60 years',
        'base_progression': 'Moderate',
        'base_severity': 6,
        'base_impact': 0.70,
        'notes': 'May have relatively slower progression'
    },
    'VCP': {
        'condition': 'Amyotrophic lateral sclerosis type 14 (ALS14) - Multisystem Proteinopathy',
        'als_subtype': 'ALS14_Multisystem',
        'onset_type': 'Adult',
        'onset_age': '40-60 years',
        'base_progression': 'Moderate',
        'base_severity': 7,
        'base_impact': 0.75,
        'notes': 'Part of multisystem proteinopathy; may include myopathy, Paget disease, FTD'
    },
    'UBQLN2': {
        'condition': 'Amyotrophic lateral sclerosis type 15 (ALS15) - X-linked',
        'als_subtype': 'ALS15_XLinked',
        'onset_type': 'Adult',
        'onset_age': '20-50 years',
        'base_progression': 'Moderate',
        'base_severity': 7,
        'base_impact': 0.75,
        'notes': 'X-linked dominant; males more severely affected'
    },
    'PFN1': {
        'condition': 'Amyotrophic lateral sclerosis type 18 (ALS18)',
        'als_subtype': 'ALS18_Classical',
        'onset_type': 'Adult',
        'onset_age': '40-60 years',
        'base_progression': 'Moderate',
        'base_severity': 7,
        'base_impact': 0.75,
        'notes': 'Rare; actin cytoskeleton dysfunction'
    },
    'TBK1': {
        'condition': 'Amyotrophic lateral sclerosis-frontotemporal dementia (ALS-FTD4)',
        'als_subtype': 'ALS_FTD',
        'onset_type': 'Adult',
        'onset_age': '50-70 years',
        'base_progression': 'Variable',
        'base_severity': 7,
        'base_impact': 0.78,
        'notes': 'Often late onset; cognitive involvement common'
    },
    'C9orf72': {
        'condition': 'Amyotrophic lateral sclerosis-frontotemporal dementia (C9orf72)',
        'als_subtype': 'ALS_FTD',
        'onset_type': 'Adult',
        'onset_age': '50-65 years',
        'base_progression': 'Variable',
        'base_severity': 8,
        'base_impact': 0.82,
        'notes': 'Most common genetic cause; 40% of familial ALS; hexanucleotide repeat'
    },
    'SQSTM1': {
        'condition': 'Amyotrophic lateral sclerosis with Paget disease of bone',
        'als_subtype': 'ALS_Multisystem',
        'onset_type': 'Adult',
        'onset_age': '50-70 years',
        'base_progression': 'Variable',
        'base_severity': 6,
        'base_impact': 0.70,
        'notes': 'May present with Paget disease; autophagy pathway'
    },
    'CHCHD10': {
        'condition': 'Amyotrophic lateral sclerosis, mitochondrial type',
        'als_subtype': 'ALS_Mitochondrial',
        'onset_type': 'Adult',
        'onset_age': '40-60 years',
        'base_progression': 'Variable',
        'base_severity': 6,
        'base_impact': 0.70,
        'notes': 'Mitochondrial dysfunction; may have myopathy'
    },
    'NEK1': {
        'condition': 'Amyotrophic lateral sclerosis susceptibility (NEK1)',
        'als_subtype': 'ALS_Risk',
        'onset_type': 'Adult',
        'onset_age': '50-65 years',
        'base_progression': 'Variable',
        'base_severity': 6,
        'base_impact': 0.65,
        'notes': 'Risk modifier; DNA damage response pathway'
    },
    'FIG4': {
        'condition': 'Amyotrophic lateral sclerosis type 11 (ALS11)',
        'als_subtype': 'ALS11_Classical',
        'onset_type': 'Adult',
        'onset_age': '40-70 years',
        'base_progression': 'Moderate',
        'base_severity': 6,
        'base_impact': 0.70,
        'notes': 'Phosphoinositide metabolism'
    },
    'ANG': {
        'condition': 'Amyotrophic lateral sclerosis type 9 (ALS9)',
        'als_subtype': 'ALS9_Classical',
        'onset_type': 'Adult',
        'onset_age': '50-70 years',
        'base_progression': 'Variable',
        'base_severity': 6,
        'base_impact': 0.68,
        'notes': 'Rare; angiogenin pathway'
    },
    'DCTN1': {
        'condition': 'Perry syndrome / ALS-Parkinsonism-Dementia complex',
        'als_subtype': 'ALS_Parkinsonism',
        'onset_type': 'Adult',
        'onset_age': '40-60 years',
        'base_progression': 'Moderate',
        'base_severity': 7,
        'base_impact': 0.76,
        'notes': 'May present with parkinsonism, depression, weight loss, hypoventilation'
    },
    
    # === Juvenile ALS Types (Early Onset, Often Slower) ===
    'ALS2': {
        'condition': 'Amyotrophic lateral sclerosis type 2 (ALS2) - Juvenile',
        'als_subtype': 'ALS2_Juvenile',
        'onset_type': 'Juvenile',
        'onset_age': '3-20 years',
        'base_progression': 'Slow',
        'base_severity': 5,
        'base_impact': 0.55,
        'notes': 'Autosomal recessive; ascending paralysis; slower progression'
    },
    'SETX': {
        'condition': 'Amyotrophic lateral sclerosis type 4 (ALS4) - Juvenile',
        'als_subtype': 'ALS4_Juvenile',
        'onset_type': 'Juvenile',
        'onset_age': '<25 years',
        'base_progression': 'Slow',
        'base_severity': 4,
        'base_impact': 0.50,
        'notes': 'Slow progression; normal lifespan possible; distal predominant'
    },
    'SPG11': {
        'condition': 'Amyotrophic lateral sclerosis type 5 (ALS5) - Juvenile',
        'als_subtype': 'ALS5_Juvenile',
        'onset_type': 'Juvenile',
        'onset_age': '7-23 years',
        'base_progression': 'Slow',
        'base_severity': 5,
        'base_impact': 0.55,
        'notes': 'May have thin corpus callosum; cognitive involvement possible'
    },
    'SIGMAR1': {
        'condition': 'Amyotrophic lateral sclerosis type 16 (ALS16) - Juvenile',
        'als_subtype': 'ALS16_Juvenile',
        'onset_type': 'Juvenile',
        'onset_age': '1-10 years',
        'base_progression': 'Slow',
        'base_severity': 5,
        'base_impact': 0.52,
        'notes': 'Very early onset possible; ER stress response'
    },
    
    # === Other Motor Neuron Disease Genes ===
    'VAPB': {
        'condition': 'Amyotrophic lateral sclerosis type 8 (ALS8)',
        'als_subtype': 'ALS8_Atypical',
        'onset_type': 'Adult',
        'onset_age': '30-50 years',
        'base_progression': 'Slow',
        'base_severity': 5,
        'base_impact': 0.55,
        'notes': 'Brazilian founder mutation; atypical presentation; slower course'
    },
    'PLEKHG5': {
        'condition': 'Spinal muscular atrophy, lower extremity-predominant (SMA-LED)',
        'als_subtype': 'SMA_Related',
        'onset_type': 'Juvenile',
        'onset_age': '10-20 years',
        'base_progression': 'Slow',
        'base_severity': 4,
        'base_impact': 0.45,
        'notes': 'Lower motor neuron predominant; legs affected first'
    },
    'SLC52A3': {
        'condition': 'Brown-Vialetto-Van Laere syndrome 1 (BVVLS1)',
        'als_subtype': 'BVVL_Juvenile',
        'onset_type': 'Juvenile',
        'onset_age': '1-20 years',
        'base_progression': 'Variable',
        'base_severity': 6,
        'base_impact': 0.65,
        'notes': 'Bulbar palsy with hearing loss; TREATABLE with high-dose riboflavin'
    },
    'SLC52A2': {
        'condition': 'Brown-Vialetto-Van Laere syndrome 2 (BVVLS2)',
        'als_subtype': 'BVVL_Juvenile',
        'onset_type': 'Juvenile',
        'onset_age': '1-20 years',
        'base_progression': 'Variable',
        'base_severity': 6,
        'base_impact': 0.65,
        'notes': 'Similar to BVVLS1; TREATABLE with riboflavin'
    },
    'PARK7': {
        'condition': 'Parkinson disease 7 / ALS modifier',
        'als_subtype': 'ALS_Modifier',
        'onset_type': 'Adult',
        'onset_age': '30-50 years',
        'base_progression': 'Variable',
        'base_severity': 5,
        'base_impact': 0.55,
        'notes': 'Primarily Parkinson disease; may modify ALS phenotype'
    },
    'CHMP2B': {
        'condition': 'Amyotrophic lateral sclerosis type 17 (ALS17)',
        'als_subtype': 'ALS17_FTD',
        'onset_type': 'Adult',
        'onset_age': '40-60 years',
        'base_progression': 'Moderate',
        'base_severity': 7,
        'base_impact': 0.75,
        'notes': 'Often associated with frontotemporal dementia'
    },
    'ERBB4': {
        'condition': 'Amyotrophic lateral sclerosis type 19 (ALS19)',
        'als_subtype': 'ALS19_Classical',
        'onset_type': 'Adult',
        'onset_age': '40-60 years',
        'base_progression': 'Slow',
        'base_severity': 5,
        'base_impact': 0.55,
        'notes': 'May have slower progression than typical ALS'
    },
    'HNRNPA1': {
        'condition': 'Amyotrophic lateral sclerosis type 20 (ALS20)',
        'als_subtype': 'ALS20_Multisystem',
        'onset_type': 'Adult',
        'onset_age': '40-70 years',
        'base_progression': 'Variable',
        'base_severity': 6,
        'base_impact': 0.68,
        'notes': 'Part of multisystem proteinopathy spectrum'
    },
    'MATR3': {
        'condition': 'Amyotrophic lateral sclerosis type 21 (ALS21)',
        'als_subtype': 'ALS21_Classical',
        'onset_type': 'Adult',
        'onset_age': '40-60 years',
        'base_progression': 'Slow',
        'base_severity': 5,
        'base_impact': 0.55,
        'notes': 'Often slower progression; may have distal myopathy'
    },
    'TUBA4A': {
        'condition': 'Amyotrophic lateral sclerosis type 22 (ALS22)',
        'als_subtype': 'ALS22_Classical',
        'onset_type': 'Adult',
        'onset_age': '50-70 years',
        'base_progression': 'Moderate',
        'base_severity': 6,
        'base_impact': 0.68,
        'notes': 'Rare; microtubule dynamics'
    },
    'ANXA11': {
        'condition': 'Amyotrophic lateral sclerosis type 23 (ALS23)',
        'als_subtype': 'ALS23_Classical',
        'onset_type': 'Adult',
        'onset_age': '50-70 years',
        'base_progression': 'Moderate',
        'base_severity': 6,
        'base_impact': 0.68,
        'notes': 'Recently identified; RNA granule transport'
    },
}

# Default phenotype for genes not in the knowledge base
DEFAULT_PHENOTYPE = {
    'condition': 'Amyotrophic lateral sclerosis, unspecified type',
    'als_subtype': 'ALS_Unspecified',
    'onset_type': 'Adult',
    'onset_age': 'Variable',
    'base_progression': 'Variable',
    'base_severity': 6,
    'base_impact': 0.70,
    'notes': 'Gene not in curated ALS gene list; assigned default parameters'
}

# ============================================================================
# SECTION 2: CONSEQUENCE SEVERITY MODIFIERS
# Based on predicted functional impact
# ============================================================================

CONSEQUENCE_SEVERITY_MULTIPLIERS = {
    'nonsense': 1.4,           # Premature stop codon - complete LoF
    'frameshift': 1.4,         # Reading frame disruption - complete LoF
    'splice acceptor': 1.3,    # Likely exon skipping - severe LoF
    'splice donor': 1.3,       # Likely exon skipping - severe LoF
    'start lost': 1.25,        # No protein initiation
    'stop lost': 1.2,          # Extended protein, potentially toxic
    'missense': 1.0,           # Amino acid substitution - variable
    'inframe deletion': 0.9,   # May preserve some function
    'inframe insertion': 0.9,  # May preserve some function
    'deletion': 1.15,          # Structural variant
    'duplication': 1.1,        # Copy number change
}

# ============================================================================
# SECTION 3: CONTEXTUAL INTERACTION MODEL PARAMETERS
# ============================================================================

# Thresholds for Ancestral Modifier assignment
PROTECTIVE_THRESHOLD = 1.5    # Ratio > 1.5 → Protective modifier
AGGRAVATING_THRESHOLD = 0.5   # Ratio < 0.5 → Aggravating modifier

# Modifier values
PROTECTIVE_MODIFIER = 0.8     # Reduces impact by 20%
NEUTRAL_MODIFIER = 1.0        # No modification
AGGRAVATING_MODIFIER = 1.2    # Increases impact by 20%

# Progression score thresholds
SLOW_PROGRESSION_THRESHOLD = 0.6
FAST_PROGRESSION_THRESHOLD = 0.85

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def find_input_file():
    """Locate the input ClinVar file in common directories."""
    search_paths = [
        "clinvar.cleaned.csv",
        "clinvar.cleaned.csv",
        os.path.expanduser("~/Downloads/clinvar.cleaned.csv"),
        os.path.expanduser("~/Downloads/clinvar.cleaned.csv"),
        os.path.expanduser("~/Desktop/clinvar.cleaned.csv"),
    ]
    for path in search_paths:
        if os.path.exists(path):
            return path
    return None


def clean_variant_data(df):
    """
    STEP 1: Data Ingestion & Variant Standardization
    
    Converts raw ClinVar data into a structured format by:
    - Creating unique variant identifiers
    - Standardizing genomic coordinates
    - Extracting primary gene symbols
    - Preserving clinical metadata (pathogenicity, consequence)
    """
    df = df.copy()
    
    # Remove any unnamed index columns
    df = df.drop(columns=[c for c in df.columns if 'Unnamed' in c or c == ''], errors='ignore')
    
    # Create unique variant ID
    def create_variant_id(rsid, idx):
        if pd.isna(rsid) or str(rsid).strip() in [':', '', 'nan']:
            return f"var_{idx:04d}"
        return str(rsid).strip()
    
    df['variant_id'] = [create_variant_id(rsid, idx) for idx, rsid in enumerate(df['rsID'])]
    
    # Clean genomic position (handle ranges like "12345 - 12350")
    def clean_position(pos):
        if pd.isna(pos) or str(pos).strip() == '':
            return None
        pos_str = str(pos).strip()
        if ' - ' in pos_str:
            return int(pos_str.split(' - ')[0])
        try:
            return int(float(pos_str))
        except ValueError:
            return None
    
    df['position_clean'] = df['position'].apply(clean_position)
    
    # Extract primary gene symbol (handle multiple genes separated by |)
    df['gene_primary'] = df['gene'].apply(
        lambda x: str(x).split('|')[0].strip().replace('-DT', '') if pd.notna(x) else 'Unknown'
    )
    
    # Clean chromosome notation
    df['chrom_clean'] = df['chromosome'].apply(
        lambda x: str(x).replace('chr', '').strip() if pd.notna(x) else ''
    )
    
    # Filter variants with valid positions and genes
    initial_count = len(df)
    df = df[df['position_clean'].notna()].copy()
    df = df[df['gene_primary'] != ''].copy()
    df = df[df['gene_primary'] != 'Unknown'].copy()
    
    print(f"       Standardized {initial_count} → {len(df)} variants (removed {initial_count - len(df)} incomplete)")
    
    return df


def get_gene_phenotype(gene_name):
    """
    Retrieve phenotype information for a gene from the knowledge base.
    Returns default values for genes not in the curated list.
    """
    gene_clean = str(gene_name).upper().strip()
    
    for gene_key, phenotype in GENE_PHENOTYPE_MAP.items():
        if gene_key.upper() == gene_clean:
            return phenotype.copy()
    
    return DEFAULT_PHENOTYPE.copy()


def get_consequence_multiplier(consequence):
    """
    Calculate severity multiplier based on molecular consequence.
    Loss-of-function variants receive higher multipliers.
    """
    consequence_lower = str(consequence).lower()
    
    for cons_type, multiplier in CONSEQUENCE_SEVERITY_MULTIPLIERS.items():
        if cons_type in consequence_lower:
            return multiplier
    
    return 1.0  # Default for unrecognized consequences


def estimate_gnomad_af(variant_type, consequence, clinical_sig, seed_value):
    """
    STEP 2: gnomAD Allele Frequency Estimation
    
    Generates population-specific allele frequencies based on:
    A. Pathogenicity constraints (pathogenic variants are rarer)
    B. Functional impact (LoF variants under stronger selection)
    C. Population-specific multipliers (founder effects, drift)
    
    This creates the foundation for the Contextual Interaction Model.
    """
    np.random.seed(seed_value % (2**31))
    
    # A. Pathogenicity-based base frequency
    if str(clinical_sig).strip() == "Pathogenic":
        base_af = np.random.uniform(0.00005, 0.0015)  # 0.005% - 0.15%
    else:  # Likely pathogenic or other
        base_af = np.random.uniform(0.0001, 0.002)    # 0.01% - 0.2%
    
    # B. Functional impact penalty for LoF variants
    lof_consequences = ['nonsense', 'frameshift', 'splice donor', 'splice acceptor', 'start lost']
    if any(lof in str(consequence).lower() for lof in lof_consequences):
        base_af *= 0.5  # 50% reduction for LoF (stronger purifying selection)
    elif 'missense' in str(consequence).lower():
        base_af *= 1.2  # Missense may be slightly more tolerated
    
    # C. Population-specific multipliers (simulating founder effects & drift)
    # These create the variation needed for the Contextual Interaction Model
    population_multipliers = {
        'AFR': np.random.uniform(0.7, 1.6),   # Highest genetic diversity
        'AMR': np.random.uniform(0.6, 1.4),   # Admixed population
        'EAS': np.random.uniform(0.4, 1.3),   # Strong founder effects
        'EUR': np.random.uniform(0.5, 1.5),   # Known ALS founder mutations (e.g., SOD1)
        'SAS': np.random.uniform(0.5, 1.4),   # Consanguinity effects
    }
    
    # Calculate population-specific AFs
    population_afs = {}
    for pop, multiplier in population_multipliers.items():
        pop_af = base_af * multiplier * np.random.uniform(0.85, 1.15)
        pop_af = max(0.000005, min(0.05, pop_af))  # Bound: 0.0005% to 5%
        population_afs[f'gnomAD_AF_{pop}'] = pop_af
    
    # Calculate global average AF
    population_afs['gnomAD_AF'] = np.mean(list(population_afs.values()))
    
    return population_afs


def simulate_genotypes_hwe(allele_frequency, n_samples):
    """
    STEP 3: Hardy-Weinberg Equilibrium Genotype Simulation
    
    Translates allele frequencies into individual genotypes using HWE:
    - Homozygous Reference (0): q² where q = 1 - p
    - Heterozygous (1): 2pq
    - Homozygous Alternate (2): p² where p = allele frequency
    
    For rare pathogenic variants, most carriers will be heterozygous.
    """
    if allele_frequency <= 0 or allele_frequency >= 1:
        allele_frequency = 0.0001
    
    p = allele_frequency      # Alternate allele frequency
    q = 1 - p                 # Reference allele frequency
    
    # HWE genotype probabilities
    prob_hom_ref = q ** 2           # 0/0
    prob_het = 2 * p * q            # 0/1
    prob_hom_alt = p ** 2           # 1/1
    
    # Normalize probabilities
    total = prob_hom_ref + prob_het + prob_hom_alt
    probabilities = [prob_hom_ref/total, prob_het/total, prob_hom_alt/total]
    
    return np.random.choice([0, 1, 2], size=n_samples, p=probabilities)


def calculate_relative_allelic_ratio(pop_af, global_af):
    """
    Calculate the Relative Allelic Ratio for the Contextual Interaction Model.
    
    Ratio = Population AF / Global AF
    
    - Ratio > 1.5: Variant more common in this population (potential tolerance)
    - Ratio < 0.5: Variant rare in this population (potential sensitivity)
    - Otherwise: Neutral
    """
    if global_af <= 0:
        return 1.0
    return pop_af / global_af


def assign_ancestral_modifier(allelic_ratio):
    """
    STEP 5: Contextual Interaction Model - Ancestral Modifier Assignment
    
    Based on the Relative Allelic Ratio, assigns a clinical modifier:
    
    - PROTECTIVE (0.8): Ratio > 1.5
      Rationale: Higher frequency suggests population may have evolved
      genetic modifiers that buffer the pathogenic effect.
    
    - AGGRAVATING (1.2): Ratio < 0.5
      Rationale: Lower frequency suggests population is more sensitive
      to the variant's deleterious effects (stronger purifying selection).
    
    - NEUTRAL (1.0): 0.5 ≤ Ratio ≤ 1.5
      Rationale: Frequency similar to global average; no modifier effect.
    """
    if allelic_ratio > PROTECTIVE_THRESHOLD:
        return PROTECTIVE_MODIFIER, 'Protective'
    elif allelic_ratio < AGGRAVATING_THRESHOLD:
        return AGGRAVATING_MODIFIER, 'Aggravating'
    else:
        return NEUTRAL_MODIFIER, 'Neutral'


def calculate_interaction_score(base_impact, ancestral_modifier, consequence_modifier):
    """
    Calculate the final Interaction Score.
    
    Interaction Score = Base Impact × Ancestral Modifier × Consequence Modifier
    
    This score determines the predicted disease progression.
    """
    return base_impact * ancestral_modifier * consequence_modifier


def assign_progression_from_score(interaction_score):
    """
    Convert Interaction Score to clinical progression category.
    
    - Score < 0.6: Slow Progression
    - Score 0.6-0.85: Moderate Progression
    - Score > 0.85: Fast Progression
    """
    if interaction_score < SLOW_PROGRESSION_THRESHOLD:
        return 'Slow'
    elif interaction_score > FAST_PROGRESSION_THRESHOLD:
        return 'Fast'
    else:
        return 'Moderate'


def calculate_severity_score(base_severity, consequence, clinical_sig):
    """
    Calculate numeric severity score (1-10 scale).
    """
    score = base_severity
    score *= get_consequence_multiplier(consequence)
    
    if str(clinical_sig).strip().lower() == 'pathogenic':
        score *= 1.1
    
    return min(10.0, max(1.0, round(score, 1)))


def severity_to_category(score):
    """Convert numeric severity to category."""
    if score >= 8:
        return 'Severe'
    elif score >= 6:
        return 'Moderate'
    elif score >= 4:
        return 'Mild'
    else:
        return 'Low'


# ============================================================================
# MAIN PIPELINE
# ============================================================================

def main():
    """
    Main execution pipeline for ALS synthetic patient generation.
    """
    
    print("=" * 80)
    print("ALS FEDERATED LEARNING SYNTHETIC PATIENT GENERATOR")
    print("Contextual Penetrance Model v2.0")
    print("=" * 80)
    
    print(f"\n[Configuration]")
    print(f"  • Patients per population: {N_PATIENTS_PER_POP:,}")
    print(f"  • Total patients: {N_PATIENTS_PER_POP * len(POPULATIONS):,}")
    print(f"  • Populations: {', '.join(POPULATIONS)}")
    print(f"  • Min carriers per variant: {MIN_CARRIERS_PER_VARIANT}")
    print(f"  • Output directory: {OUTPUT_DIR}")
    
    # Locate input file
    input_file = sys.argv[1] if len(sys.argv) > 1 else find_input_file()
    
    if not input_file or not os.path.exists(input_file):
        print("\n ERROR: Could not locate clinvar.cleaned.csv")
        print("Usage: python generate_als_synthetic_complete.py /path/to/clinvar.cleaned.csv")
        sys.exit(1)
    
    print(f"\n✓ Input file: {input_file}")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # ========================================================================
    # PHASE 1: DATA INGESTION & VARIANT STANDARDIZATION
    # ========================================================================
    print("\n" + "=" * 80)
    print("PHASE 1: DATA INGESTION & VARIANT STANDARDIZATION")
    print("=" * 80)
    
    print("\n[1.1] Loading ClinVar variant data...")
    df = pd.read_csv(input_file)
    print(f"       Raw input: {len(df)} variants")
    
    print("\n[1.2] Standardizing variant annotations...")
    df = clean_variant_data(df)
    
    # ========================================================================
    # PHASE 2: PHENOTYPE & gnomAD FREQUENCY ASSIGNMENT
    # ========================================================================
    print("\n" + "=" * 80)
    print("PHASE 2: PHENOTYPE & ALLELE FREQUENCY ESTIMATION")
    print("=" * 80)
    
    print("\n[2.1] Assigning gene-based phenotypes...")
    
    # Initialize phenotype columns
    phenotype_columns = ['condition', 'als_subtype', 'onset_type', 'onset_age', 
                         'base_progression', 'base_impact', 'phenotype_notes']
    for col in phenotype_columns:
        df[col] = ''
    df['base_severity'] = 0.0
    
    # Assign phenotypes based on gene
    for idx, row in df.iterrows():
        phenotype = get_gene_phenotype(row['gene_primary'])
        df.at[idx, 'condition'] = phenotype['condition']
        df.at[idx, 'als_subtype'] = phenotype['als_subtype']
        df.at[idx, 'onset_type'] = phenotype['onset_type']
        df.at[idx, 'onset_age'] = phenotype['onset_age']
        df.at[idx, 'base_progression'] = phenotype['base_progression']
        df.at[idx, 'base_severity'] = phenotype['base_severity']
        df.at[idx, 'base_impact'] = phenotype['base_impact']
        df.at[idx, 'phenotype_notes'] = phenotype['notes']
    
    print(f"       Assigned phenotypes to {len(df)} variants")
    print(f"       Unique subtypes: {df['als_subtype'].nunique()}")
    
    print("\n[2.2] Estimating population-specific allele frequencies...")
    
    # Initialize AF columns
    af_columns = ['gnomAD_AF', 'gnomAD_AF_AFR', 'gnomAD_AF_AMR', 
                  'gnomAD_AF_EAS', 'gnomAD_AF_EUR', 'gnomAD_AF_SAS']
    for col in af_columns:
        df[col] = 0.0
    
    # Estimate AFs
    for idx, row in df.iterrows():
        seed = hash(f"{row['variant_id']}_{row.get('consequence', '')}_{idx}") % (2**31)
        afs = estimate_gnomad_af(
            row.get('variant_type', ''),
            row.get('consequence', ''),
            row.get('clinical_sig', 'Pathogenic'),
            seed
        )
        for col, val in afs.items():
            df.at[idx, col] = val
    
    print(f"       Global AF: mean={df['gnomAD_AF'].mean():.6f}, range=[{df['gnomAD_AF'].min():.6f}, {df['gnomAD_AF'].max():.4f}]")
    
    # Calculate severity scores
    print("\n[2.3] Calculating variant severity scores...")
    df['severity_score'] = df.apply(
        lambda r: calculate_severity_score(r['base_severity'], r.get('consequence', ''), r.get('clinical_sig', '')),
        axis=1
    )
    df['severity_category'] = df['severity_score'].apply(severity_to_category)
    
    print(f"       Severity distribution: {df['severity_category'].value_counts().to_dict()}")
    
    # Save variant data
    variants_file = os.path.join(OUTPUT_DIR, "variants_complete.csv")
    df.to_csv(variants_file, index=False)
    print(f"\n       ✓ Saved: {variants_file}")
    
    # ========================================================================
    # PHASE 3: SYNTHETIC PATIENT GENERATION
    # ========================================================================
    print("\n" + "=" * 80)
    print("PHASE 3: SYNTHETIC PATIENT GENERATION (HWE Simulation)")
    print("=" * 80)
    
    np.random.seed(RANDOM_SEED)
    
    variant_ids = df['variant_id'].tolist()
    n_variants = len(df)
    
    # Pre-compute variant data for fast lookup
    variant_data = {}
    for _, row in df.iterrows():
        variant_data[row['variant_id']] = {
            'gene': row['gene_primary'],
            'condition': row['condition'],
            'als_subtype': row['als_subtype'],
            'onset_type': row['onset_type'],
            'base_progression': row['base_progression'],
            'base_impact': row['base_impact'],
            'base_severity': row['base_severity'],
            'severity_score': row['severity_score'],
            'severity_category': row['severity_category'],
            'consequence': row.get('consequence', ''),
            'gnomAD_AF': row['gnomAD_AF'],
        }
        # Add population-specific AFs
        for pop in POPULATIONS:
            variant_data[row['variant_id']][f'AF_{pop}'] = row[f'gnomAD_AF_{pop}']
    
    all_patients = []
    total_carriers_per_variant = np.zeros(n_variants, dtype=int)
    
    # Track contextual model examples for documentation
    contextual_examples = []
    
    print(f"\n[3.1] Generating patients with HWE-based genotypes...")
    
    for pop_idx, pop in enumerate(POPULATIONS):
        print(f"       {pop}: generating {N_PATIENTS_PER_POP:,} patients...", end=" ", flush=True)
        
        af_col = f'gnomAD_AF_{pop}'
        
        # Generate genotype matrix
        genotype_matrix = np.zeros((N_PATIENTS_PER_POP, n_variants), dtype=np.int8)
        
        for var_idx in range(n_variants):
            af = df.iloc[var_idx][af_col]
            if pd.isna(af) or af <= 0:
                af = 0.00001
            genotype_matrix[:, var_idx] = simulate_genotypes_hwe(af, N_PATIENTS_PER_POP)
        
        # STEP 4: Rare Variant Preservation (on last population)
        if pop_idx == len(POPULATIONS) - 1:
            current_totals = total_carriers_per_variant + np.sum(genotype_matrix > 0, axis=0)
            variants_needing_carriers = np.where(current_totals < MIN_CARRIERS_PER_VARIANT)[0]
            
            for var_idx in variants_needing_carriers:
                needed = MIN_CARRIERS_PER_VARIANT - int(current_totals[var_idx])
                non_carriers = np.where(genotype_matrix[:, var_idx] == 0)[0]
                if len(non_carriers) >= needed:
                    new_carriers = np.random.choice(non_carriers, size=needed, replace=False)
                    genotype_matrix[new_carriers, var_idx] = 1
        
        # Update carrier counts
        total_carriers_per_variant += np.sum(genotype_matrix > 0, axis=0)
        
        # Create patient records with Contextual Interaction Model
        for i in range(N_PATIENTS_PER_POP):
            # Find variants carried by this patient
            carried_variant_indices = np.where(genotype_matrix[i, :] > 0)[0]
            carried_variants = []
            
            for var_idx in carried_variant_indices:
                var_id = variant_ids[var_idx]
                var_info = variant_data[var_id].copy()
                var_info['genotype'] = int(genotype_matrix[i, var_idx])
                var_info['variant_id'] = var_id
                
                # === CONTEXTUAL INTERACTION MODEL ===
                # Calculate Relative Allelic Ratio
                pop_af = var_info[f'AF_{pop}']
                global_af = var_info['gnomAD_AF']
                allelic_ratio = calculate_relative_allelic_ratio(pop_af, global_af)
                
                # Assign Ancestral Modifier
                ancestral_mod, modifier_type = assign_ancestral_modifier(allelic_ratio)
                
                # Get consequence modifier
                consequence_mod = get_consequence_multiplier(var_info['consequence'])
                
                # Calculate Interaction Score
                interaction_score = calculate_interaction_score(
                    var_info['base_impact'],
                    ancestral_mod,
                    consequence_mod
                )
                
                # Determine Progression from Interaction Score
                contextual_progression = assign_progression_from_score(interaction_score)
                
                var_info['allelic_ratio'] = round(allelic_ratio, 3)
                var_info['ancestral_modifier'] = ancestral_mod
                var_info['modifier_type'] = modifier_type
                var_info['consequence_modifier'] = consequence_mod
                var_info['interaction_score'] = round(interaction_score, 3)
                var_info['contextual_progression'] = contextual_progression
                
                carried_variants.append(var_info)
                
                # Collect examples for documentation (first few interesting cases)
                if len(contextual_examples) < 20 and modifier_type != 'Neutral':
                    contextual_examples.append({
                        'patient_id': f"{pop}_{i:05d}",
                        'population': pop,
                        'variant_id': var_id,
                        'gene': var_info['gene'],
                        'pop_af': pop_af,
                        'global_af': global_af,
                        'allelic_ratio': allelic_ratio,
                        'modifier_type': modifier_type,
                        'ancestral_modifier': ancestral_mod,
                        'base_impact': var_info['base_impact'],
                        'interaction_score': interaction_score,
                        'progression': contextual_progression
                    })
            
            # Build patient record
            if carried_variants:
                # Sort by interaction score (most severe first)
                carried_variants.sort(key=lambda x: x['interaction_score'], reverse=True)
                primary = carried_variants[0]
                
                patient = {
                    'patient_id': f"{pop}_{i:05d}",
                    'superpopulation': pop,
                    'n_variants_carried': len(carried_variants),
                    # Primary variant info
                    'primary_gene': primary['gene'],
                    'primary_condition': primary['condition'],
                    'primary_subtype': primary['als_subtype'],
                    'onset_type': primary['onset_type'],
                    # Contextual Interaction Model outputs
                    'base_impact': primary['base_impact'],
                    'allelic_ratio': primary['allelic_ratio'],
                    'ancestral_modifier': primary['ancestral_modifier'],
                    'modifier_type': primary['modifier_type'],
                    'consequence_modifier': primary['consequence_modifier'],
                    'interaction_score': primary['interaction_score'],
                    'predicted_progression': primary['contextual_progression'],
                    # Severity
                    'severity_score': primary['severity_score'],
                    'severity_category': primary['severity_category'],
                    # All variants info
                    'all_genes': '; '.join(sorted(set(v['gene'] for v in carried_variants))),
                    'all_subtypes': '; '.join(sorted(set(v['als_subtype'] for v in carried_variants))),
                }
            else:
                patient = {
                    'patient_id': f"{pop}_{i:05d}",
                    'superpopulation': pop,
                    'n_variants_carried': 0,
                    'primary_gene': 'None',
                    'primary_condition': 'No pathogenic variant',
                    'primary_subtype': 'Unaffected',
                    'onset_type': 'N/A',
                    'base_impact': 0.0,
                    'allelic_ratio': 0.0,
                    'ancestral_modifier': 1.0,
                    'modifier_type': 'N/A',
                    'consequence_modifier': 1.0,
                    'interaction_score': 0.0,
                    'predicted_progression': 'N/A',
                    'severity_score': 0.0,
                    'severity_category': 'Unaffected',
                    'all_genes': '',
                    'all_subtypes': '',
                }
            
            # Add genotype columns
            for var_idx, var_id in enumerate(variant_ids):
                patient[f'geno_{var_id}'] = int(genotype_matrix[i, var_idx])
            
            all_patients.append(patient)
        
        n_carriers = (np.sum(genotype_matrix > 0, axis=0) > 0).sum()
        print(f"✓ ({np.sum(genotype_matrix > 0):,} carrier-variant pairs)")
    
    patients_df = pd.DataFrame(all_patients)
    
    # Update variant file with actual carrier counts
    df['actual_carriers'] = total_carriers_per_variant
    df.to_csv(variants_file, index=False)
    
    # ========================================================================
    # PHASE 4: OUTPUT GENERATION
    # ========================================================================
    print("\n" + "=" * 80)
    print("PHASE 4: OUTPUT GENERATION")
    print("=" * 80)
    
    # Save all patients
    print("\n[4.1] Saving patient datasets...")
    patients_file = os.path.join(OUTPUT_DIR, "patients_complete.csv")
    patients_df.to_csv(patients_file, index=False)
    print(f"       ✓ All patients: {patients_file}")
    
    # Save carriers only
    carriers_df = patients_df[patients_df['n_variants_carried'] > 0].copy()
    carriers_file = os.path.join(OUTPUT_DIR, "patients_carriers_only.csv")
    carriers_df.to_csv(carriers_file, index=False)
    print(f"       ✓ Carriers only: {carriers_file} ({len(carriers_df):,} patients)")
    
    # Save federated client files
    print("\n[4.2] Creating federated client datasets...")
    for pop in POPULATIONS:
        pop_df = patients_df[patients_df['superpopulation'] == pop]
        
        # All patients
        pop_file = os.path.join(OUTPUT_DIR, f"client_{pop}.csv")
        pop_df.to_csv(pop_file, index=False)
        
        # Carriers only
        pop_carriers = pop_df[pop_df['n_variants_carried'] > 0]
        pop_carriers_file = os.path.join(OUTPUT_DIR, f"client_{pop}_carriers.csv")
        pop_carriers.to_csv(pop_carriers_file, index=False)
        
        print(f"       ✓ {pop}: {len(pop_df):,} total, {len(pop_carriers):,} carriers")
    
    # ========================================================================
    # PHASE 5: VALIDATION & REPORTING
    # ========================================================================
    print("\n" + "=" * 80)
    print("PHASE 5: VALIDATION METRICS")
    print("=" * 80)
    
    total_patients = len(patients_df)
    total_carriers = len(carriers_df)
    
    print(f"\n[5.1] Patient Statistics:")
    print(f"       Total patients: {total_patients:,}")
    print(f"       Carriers (≥1 variant): {total_carriers:,} ({100*total_carriers/total_patients:.1f}%)")
    print(f"       Mean variants/carrier: {carriers_df['n_variants_carried'].mean():.2f}")
    
    print(f"\n[5.2] Variant Statistics:")
    print(f"       Total variants: {len(df):,}")
    print(f"       Variants with ≥1 carrier: {(df['actual_carriers'] > 0).sum()}")
    print(f"       Mean carriers/variant: {df['actual_carriers'].mean():.1f}")
    
    print(f"\n[5.3] AF vs Carrier Correlation (validation):")
    correlation = df[['gnomAD_AF', 'actual_carriers']].corr().iloc[0, 1]
    print(f"       Pearson r = {correlation:.3f}")
    print(f"       {'✓ Strong positive correlation confirms HWE simulation accuracy' if correlation > 0.7 else '⚠ Weaker correlation'}")
    
    print(f"\n[5.4] Subtype Distribution (carriers):")
    print(carriers_df['primary_subtype'].value_counts().head(10).to_string())
    
    print(f"\n[5.5] Progression Distribution (carriers):")
    print(carriers_df['predicted_progression'].value_counts().to_string())
    
    print(f"\n[5.6] Modifier Type Distribution (carriers):")
    print(carriers_df['modifier_type'].value_counts().to_string())
    
    # ========================================================================
    # GENERATE DOCUMENTATION
    # ========================================================================
    
    # Validation report
    validation_file = os.path.join(OUTPUT_DIR, "validation_report.txt")
    with open(validation_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("ALS SYNTHETIC PATIENT GENERATION - VALIDATION REPORT\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("CONFIGURATION\n")
        f.write("-" * 40 + "\n")
        f.write(f"Patients per population: {N_PATIENTS_PER_POP:,}\n")
        f.write(f"Total patients: {total_patients:,}\n")
        f.write(f"Populations: {', '.join(POPULATIONS)}\n")
        f.write(f"Random seed: {RANDOM_SEED}\n\n")
        
        f.write("PATIENT STATISTICS\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total patients: {total_patients:,}\n")
        f.write(f"Carriers: {total_carriers:,} ({100*total_carriers/total_patients:.1f}%)\n")
        f.write(f"Mean variants per carrier: {carriers_df['n_variants_carried'].mean():.2f}\n\n")
        
        f.write("VARIANT STATISTICS\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total variants: {len(df):,}\n")
        f.write(f"Variants with carriers: {(df['actual_carriers'] > 0).sum()}\n")
        f.write(f"Mean carriers per variant: {df['actual_carriers'].mean():.1f}\n")
        f.write(f"AF-Carrier correlation: {correlation:.3f}\n\n")
        
        f.write("SUBTYPE DISTRIBUTION\n")
        f.write("-" * 40 + "\n")
        f.write(carriers_df['primary_subtype'].value_counts().to_string() + "\n\n")
        
        f.write("PROGRESSION DISTRIBUTION\n")
        f.write("-" * 40 + "\n")
        f.write(carriers_df['predicted_progression'].value_counts().to_string() + "\n\n")
        
        f.write("MODIFIER TYPE DISTRIBUTION\n")
        f.write("-" * 40 + "\n")
        f.write(carriers_df['modifier_type'].value_counts().to_string() + "\n")
    
    print(f"\n       ✓ Validation report: {validation_file}")
    
    # Contextual model examples
    examples_file = os.path.join(OUTPUT_DIR, "contextual_model_examples.txt")
    with open(examples_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("CONTEXTUAL INTERACTION MODEL - WORKED EXAMPLES\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("This document demonstrates how population-specific allele frequencies\n")
        f.write("modulate disease severity through the Contextual Penetrance Model.\n\n")
        
        f.write("MODEL PARAMETERS:\n")
        f.write(f"  • Protective threshold: Ratio > {PROTECTIVE_THRESHOLD}\n")
        f.write(f"  • Aggravating threshold: Ratio < {AGGRAVATING_THRESHOLD}\n")
        f.write(f"  • Protective modifier: {PROTECTIVE_MODIFIER}\n")
        f.write(f"  • Aggravating modifier: {AGGRAVATING_MODIFIER}\n")
        f.write(f"  • Slow progression threshold: Score < {SLOW_PROGRESSION_THRESHOLD}\n")
        f.write(f"  • Fast progression threshold: Score > {FAST_PROGRESSION_THRESHOLD}\n\n")
        
        f.write("=" * 80 + "\n")
        f.write("EXAMPLES:\n")
        f.write("=" * 80 + "\n\n")
        
        for i, ex in enumerate(contextual_examples[:10], 1):
            f.write(f"EXAMPLE {i}: {ex['patient_id']}\n")
            f.write("-" * 40 + "\n")
            f.write(f"  Population: {ex['population']}\n")
            f.write(f"  Variant: {ex['variant_id']} ({ex['gene']})\n")
            f.write(f"  Population AF: {ex['pop_af']:.6f}\n")
            f.write(f"  Global AF: {ex['global_af']:.6f}\n")
            f.write(f"  Allelic Ratio: {ex['allelic_ratio']:.3f}\n")
            f.write(f"  Modifier Type: {ex['modifier_type']}\n")
            f.write(f"  Ancestral Modifier: {ex['ancestral_modifier']}\n")
            f.write(f"  Base Impact: {ex['base_impact']}\n")
            f.write(f"  Interaction Score: {ex['interaction_score']:.3f}\n")
            f.write(f"  Predicted Progression: {ex['progression']}\n\n")
            
            # Explanation
            if ex['modifier_type'] == 'Protective':
                f.write(f"  INTERPRETATION: Variant is {ex['allelic_ratio']:.1f}x MORE common in {ex['population']}\n")
                f.write(f"  than the global average. This suggests the population may have\n")
                f.write(f"  evolved genetic modifiers that buffer the pathogenic effect.\n")
                f.write(f"  Result: 20% REDUCED disease impact → {ex['progression']} progression\n\n")
            else:
                f.write(f"  INTERPRETATION: Variant is {1/ex['allelic_ratio']:.1f}x LESS common in {ex['population']}\n")
                f.write(f"  than the global average. This suggests stronger purifying selection,\n")
                f.write(f"  indicating the population is more sensitive to this variant.\n")
                f.write(f"  Result: 20% INCREASED disease impact → {ex['progression']} progression\n\n")
    
    print(f"       ✓ Contextual model examples: {examples_file}")
    
    # ========================================================================
    # FINAL SUMMARY
    # ========================================================================
    print("\n" + "=" * 80)
    print("✓ GENERATION COMPLETE")
    print("=" * 80)
    
    print(f"\n Output files in: {OUTPUT_DIR}/")
    print(f"   • variants_complete.csv")
    print(f"   • patients_complete.csv")
    print(f"   • patients_carriers_only.csv")
    for pop in POPULATIONS:
        print(f"   • client_{pop}.csv / client_{pop}_carriers.csv")
    print(f"   • validation_report.txt")
    print(f"   • contextual_model_examples.txt")
    
    return df, patients_df


if __name__ == "__main__":
    variants_df, patients_df = main()