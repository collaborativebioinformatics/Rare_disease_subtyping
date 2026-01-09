#!/usr/bin/env python3
"""
ALS Synthetic Patient Generator - Extended Version

Features:
- 15,000 synthetic patients (3,000 per population)
- Every variant guaranteed to have at least 1 carrier
- Number of carriers reflects allele frequency (higher AF = more carriers)

Usage: python generate_synthetic_patients.py
"""

import pandas as pd
import numpy as np
import sys
import os
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

# Number of patients PER POPULATION
# Total patients = N_PATIENTS_PER_POP × 5 populations
# Example: 3000 × 5 = 15,000 total patients
N_PATIENTS_PER_POP = 3000  

# Minimum carriers per variant (guarantees every variant appears in at least this many patients)
MIN_CARRIERS_PER_VARIANT = 1  # <-- CHANGE THIS (e.g., 5 for at least 5 carriers each)

# Populations (federated clients)
POPULATIONS = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']

# Output directory
OUTPUT_DIR = os.path.expanduser("~/Downloads/Synthetic_data")

# Random seed for reproducibility
RANDOM_SEED = 42

# ============================================================================
# FILE FINDER
# ============================================================================

def find_input_file():
    """Automatically find the input file in common locations."""
    possible_paths = [
        "clinvar.cleaned.csv",
        os.path.expanduser("~/Downloads/clinvar.cleaned.csv"),
        os.path.expanduser("~/Desktop/clinvar.cleaned.csv"),
        os.path.expanduser("~/Documents/clinvar.cleaned.csv"),
    ]
    
    for path in possible_paths:
        if os.path.exists(path):
            return path
    return None

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def estimate_gnomad_af(variant_type, consequence, clinical_sig, seed_val):
    """
    Estimate gnomAD allele frequencies based on variant characteristics.
    
    Returns population-specific AFs that will determine how many 
    patients carry each variant.
    """
    np.random.seed(seed_val % (2**31))
    
    # Base frequency - pathogenic variants are rare
    if clinical_sig == "Pathogenic":
        base_af = np.random.uniform(0.0001, 0.002)  # 0.01% - 0.2%
    else:  # Likely pathogenic
        base_af = np.random.uniform(0.0002, 0.003)  # 0.02% - 0.3%
    
    # Loss-of-function variants are typically rarer
    lof_consequences = ['nonsense', 'frameshift', 'splice donor', 'splice acceptor']
    if any(lof in str(consequence).lower() for lof in lof_consequences):
        base_af *= 0.5
    
    # Missense variants can be slightly more common
    if 'missense' in str(consequence).lower():
        base_af *= 1.5
    
    # Population-specific variation (founder effects, drift)
    pop_multipliers = {
        'AFR': np.random.uniform(0.8, 1.5),   # Highest diversity
        'AMR': np.random.uniform(0.7, 1.3),   # Admixed
        'EAS': np.random.uniform(0.5, 1.2),   # Some founder effects
        'EUR': np.random.uniform(0.6, 1.4),   # Founder effects (e.g., SOD1 in Scandinavia)
        'SAS': np.random.uniform(0.6, 1.3),   # Some consanguinity effects
    }
    
    afs = {}
    for pop, mult in pop_multipliers.items():
        af = base_af * mult * np.random.uniform(0.8, 1.2)
        af = max(0.00001, min(0.05, af))  # Bound between 0.001% and 5%
        afs[f'gnomAD_AF_{pop}'] = af
    
    afs['gnomAD_AF'] = np.mean(list(afs.values()))
    return afs


def simulate_genotypes_hwe(af, n_samples):
    """
    Simulate genotypes using Hardy-Weinberg Equilibrium.
    
    Higher AF → More carriers (this reflects AF in patient counts)
    
    Returns: array of genotypes (0=ref/ref, 1=het, 2=hom_alt)
    """
    if af <= 0 or af >= 1:
        af = 0.0001
    
    p = af          # Alt allele frequency
    q = 1 - p       # Ref allele frequency
    
    # HWE genotype frequencies
    prob_0 = q ** 2           # Homozygous reference
    prob_1 = 2 * p * q        # Heterozygous (most carriers will be het for rare variants)
    prob_2 = p ** 2           # Homozygous alternate (very rare)
    
    # Normalize
    total = prob_0 + prob_1 + prob_2
    probs = [prob_0/total, prob_1/total, prob_2/total]
    
    return np.random.choice([0, 1, 2], size=n_samples, p=probs)


def ensure_minimum_carriers(genotype_matrix, variant_afs, min_carriers=1):
    """
    Ensure every variant has at least `min_carriers` patients.
    
    For variants with 0 carriers after HWE simulation, 
    randomly assign carriers proportional to AF.
    """
    n_patients, n_variants = genotype_matrix.shape
    
    for var_idx in range(n_variants):
        current_carriers = np.sum(genotype_matrix[:, var_idx] > 0)
        
        if current_carriers < min_carriers:
            # Need to add more carriers
            needed = min_carriers - current_carriers
            
            # Find patients who don't carry this variant
            non_carriers = np.where(genotype_matrix[:, var_idx] == 0)[0]
            
            if len(non_carriers) >= needed:
                # Randomly select patients to become carriers
                new_carriers = np.random.choice(non_carriers, size=needed, replace=False)
                
                # Assign heterozygous genotype (most realistic for rare variants)
                genotype_matrix[new_carriers, var_idx] = 1
    
    return genotype_matrix


def clean_variant_data(df):
    """Clean and standardize the variant data."""
    df = df.copy()
    df = df.drop(columns=[c for c in df.columns if 'Unnamed' in c or c == ''], errors='ignore')
    
    def clean_rsid(rsid, idx):
        if pd.isna(rsid) or rsid == ':' or rsid == '':
            return f"var_{idx}"
        return str(rsid).strip()
    
    df['variant_id'] = [clean_rsid(rsid, idx) for idx, rsid in enumerate(df['rsID'])]
    
    def clean_position(pos):
        if pd.isna(pos) or pos == '':
            return None
        pos_str = str(pos)
        if ' - ' in pos_str:
            return int(pos_str.split(' - ')[0])
        try:
            return int(float(pos_str))
        except:
            return None
    
    df['position_clean'] = df['position'].apply(clean_position)
    df['gene_primary'] = df['gene'].apply(lambda x: str(x).split('|')[0] if pd.notna(x) else 'Unknown')
    df['chrom_clean'] = df['chromosome'].apply(lambda x: str(x).replace('chr', '') if pd.notna(x) else '')
    
    df = df[df['position_clean'].notna()].copy()
    df = df[df['gene_primary'] != ''].copy()
    
    return df


# ============================================================================
# MAIN FUNCTION
# ============================================================================

def main():
    print("=" * 70)
    print("ALS SYNTHETIC PATIENT GENERATOR - EXTENDED VERSION")
    print("=" * 70)
    print(f"\nSettings:")
    print(f"  • Patients per population: {N_PATIENTS_PER_POP:,}")
    print(f"  • Total patients: {N_PATIENTS_PER_POP * len(POPULATIONS):,}")
    print(f"  • Populations: {POPULATIONS}")
    print(f"  • Min carriers per variant: {MIN_CARRIERS_PER_VARIANT}")
    
    # Find input file
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    else:
        input_file = find_input_file()
    
    if input_file is None or not os.path.exists(input_file):
        print("\n ERROR: Could not find clinvar_cleaned.csv")
        print("\nPlease run with the path:")
        print("  python generate_synthetic_patients.py /path/to/clinvar_cleaned.csv")
        sys.exit(1)
    
    print(f"\n✓ Found input file: {input_file}")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Load data
    print(f"\n[1/6] Loading ClinVar variants...")
    df = pd.read_csv(input_file)
    print(f"       Loaded {len(df)} variants")
    
    # Clean data
    print("\n[2/6] Cleaning variant data...")
    df = clean_variant_data(df)
    print(f"       After cleaning: {len(df)} variants")
    
    # Add gnomAD frequencies
    print("\n[3/6] Estimating gnomAD allele frequencies...")
    af_columns = ['gnomAD_AF', 'gnomAD_AF_AFR', 'gnomAD_AF_AMR', 
                  'gnomAD_AF_EAS', 'gnomAD_AF_EUR', 'gnomAD_AF_SAS']
    
    for col in af_columns:
        df[col] = 0.0
    
    for idx, row in df.iterrows():
        seed_val = hash(f"{row.get('variant_type', '')}_{row.get('consequence', '')}_{idx}")
        afs = estimate_gnomad_af(
            str(row.get('variant_type', '')),
            str(row.get('consequence', '')),
            str(row.get('clinical_sig', 'Pathogenic')),
            seed_val
        )
        for col, val in afs.items():
            df.at[idx, col] = val
    
    print(f"       Mean AF: {df['gnomAD_AF'].mean():.4f} ({df['gnomAD_AF'].mean()*100:.2f}%)")
    print(f"       AF range: {df['gnomAD_AF'].min():.6f} - {df['gnomAD_AF'].max():.4f}")
    
    # Calculate expected carriers per variant
    df['expected_carriers'] = df['gnomAD_AF'].apply(
        lambda af: int(2 * af * (1-af) * N_PATIENTS_PER_POP * len(POPULATIONS))
    )
    
    # Save variants
    variants_file = os.path.join(OUTPUT_DIR, "variants_with_gnomad.csv")
    df.to_csv(variants_file, index=False)
    print(f"       ✓ Saved: {variants_file}")
    
    # Generate synthetic patients
    print(f"\n[4/6] Generating {N_PATIENTS_PER_POP * len(POPULATIONS):,} synthetic patients...")
    np.random.seed(RANDOM_SEED)
    
    all_patients = []
    variant_ids = df['variant_id'].tolist()
    n_variants = len(df)
    
    # Track carriers per variant across all populations
    total_carriers_per_variant = np.zeros(n_variants, dtype=int)
    
    for pop_idx, pop in enumerate(POPULATIONS):
        print(f"       Generating {pop} ({N_PATIENTS_PER_POP:,} patients)...", end=" ", flush=True)
        
        af_col = f'gnomAD_AF_{pop}'
        
        # Generate genotype matrix using HWE
        genotype_matrix = np.zeros((N_PATIENTS_PER_POP, n_variants), dtype=np.int8)
        
        for var_idx, (_, variant) in enumerate(df.iterrows()):
            af = variant.get(af_col, 0.0001)
            if pd.isna(af) or af <= 0:
                af = 0.0001
            genotype_matrix[:, var_idx] = simulate_genotypes_hwe(af, N_PATIENTS_PER_POP)
        
        # Ensure minimum carriers (only for last population to check totals)
        if pop_idx == len(POPULATIONS) - 1:
            # Check which variants still need carriers
            current_total = total_carriers_per_variant + np.sum(genotype_matrix > 0, axis=0)
            variants_needing_carriers = np.where(current_total < MIN_CARRIERS_PER_VARIANT)[0]
            
            for var_idx in variants_needing_carriers:
                needed = MIN_CARRIERS_PER_VARIANT - current_total[var_idx]
                non_carriers = np.where(genotype_matrix[:, var_idx] == 0)[0]
                if len(non_carriers) >= needed:
                    new_carriers = np.random.choice(non_carriers, size=int(needed), replace=False)
                    genotype_matrix[new_carriers, var_idx] = 1
        
        # Update total carriers
        total_carriers_per_variant += np.sum(genotype_matrix > 0, axis=0)
        
        # Create patient records
        for i in range(N_PATIENTS_PER_POP):
            patient = {
                'patient_id': f"{pop}_{i:05d}",
                'superpopulation': pop,
                'n_variants_carried': int(np.sum(genotype_matrix[i, :] > 0)),
            }
            
            # Add genotypes for each variant
            for var_idx, var_id in enumerate(variant_ids):
                patient[f"geno_{var_id}"] = int(genotype_matrix[i, var_idx])
            
            all_patients.append(patient)
        
        carriers_this_pop = np.sum(genotype_matrix > 0)
        print(f"✓ ({carriers_this_pop:,} variant-carrier pairs)")
    
    patients_df = pd.DataFrame(all_patients)
    
    # Add carrier counts to variants dataframe
    df['actual_carriers'] = total_carriers_per_variant
    df.to_csv(variants_file, index=False)  # Update with actual carriers
    
    # [5/6] Saving patient data
    print(f"\n[5/6] Saving patient data...")
    patients_file = os.path.join(OUTPUT_DIR, "synthetic_patients.csv")
    patients_df.to_csv(patients_file, index=False)
    print(f"       ✓ Saved Original: {patients_file}")

    # ADD: Save filtered version (n_variants_carried > 0)
    filtered_patients_file = os.path.join(OUTPUT_DIR, "synthetic_patients_filtered.csv")
    filtered_df = patients_df[patients_df['n_variants_carried'] > 0]
    filtered_df.to_csv(filtered_patients_file, index=False)
    print(f"       ✓ Saved Filtered: {filtered_patients_file} ({len(filtered_df):,} patients)")

    # [6/6] Creating federated   datasets
    print(f"\n[6/6] Creating federated   datasets...")
    for pop in POPULATIONS:
    # Keep original   files as they were
     _df = patients_df[patients_df['superpopulation'] == pop].copy()
     _file = os.path.join(OUTPUT_DIR, f" _{pop}.csv")
     _df.to_csv( _file, index=False)
        
    # OPTIONAL: Create filtered versions of   files too
     _filtered_file = os.path.join(OUTPUT_DIR, f" _{pop}_filtered.csv")
     _df[ _df['n_variants_carried'] > 0].to_csv( _filtered_file, index=False)
        
    n_carriers = ( _df['n_variants_carried'] > 0).sum()
    print(f"     ✓ {pop}: {len( _df):,} total, {n_carriers:,} carriers saved to filtered file")
    
    
    # ========================================================================
    # SUMMARY STATISTICS
    # ========================================================================
    print("\n" + "=" * 70)
    print("✓ GENERATION COMPLETE!")
    print("=" * 70)
    
    total_patients = len(patients_df)
    total_carriers = (patients_df['n_variants_carried'] > 0).sum()
    
    print(f"\n PATIENT STATISTICS:")
    print(f"   Total patients: {total_patients:,}")
    print(f"   Patients carrying ≥1 variant: {total_carriers:,} ({100*total_carriers/total_patients:.1f}%)")
    print(f"   Mean variants per patient: {patients_df['n_variants_carried'].mean():.2f}")
    print(f"   Max variants per patient: {patients_df['n_variants_carried'].max()}")
    
    print(f"\n VARIANT STATISTICS:")
    print(f"   Total variants: {len(df):,}")
    print(f"   Variants with ≥1 carrier: {(df['actual_carriers'] > 0).sum()}")
    print(f"   Min carriers per variant: {df['actual_carriers'].min()}")
    print(f"   Max carriers per variant: {df['actual_carriers'].max()}")
    print(f"   Mean carriers per variant: {df['actual_carriers'].mean():.1f}")
    
    print(f"\n AF vs CARRIER CORRELATION:")
    # Show that higher AF = more carriers
    print("   (Higher AF should have more carriers)")
    af_carrier_corr = df[['gnomAD_AF', 'actual_carriers']].corr().iloc[0, 1]
    print(f"   Correlation(AF, carriers): {af_carrier_corr:.3f}")
    
    # Show examples
    print(f"\n   Top 5 variants by carrier count:")
    top_carriers = df.nlargest(5, 'actual_carriers')[['variant_id', 'gene_primary', 'gnomAD_AF', 'actual_carriers']]
    for _, row in top_carriers.iterrows():
        print(f"     {row['variant_id']}: {row['gene_primary']}, AF={row['gnomAD_AF']:.4f}, carriers={row['actual_carriers']}")
    
    print(f"\n   Bottom 5 variants by carrier count:")
    bottom_carriers = df.nsmallest(5, 'actual_carriers')[['variant_id', 'gene_primary', 'gnomAD_AF', 'actual_carriers']]
    for _, row in bottom_carriers.iterrows():
        print(f"     {row['variant_id']}: {row['gene_primary']}, AF={row['gnomAD_AF']:.4f}, carriers={row['actual_carriers']}")
    
    print(f"\n OUTPUT FILES in '{OUTPUT_DIR}/':")
    print(f"   • variants_with_gnomad.csv  (variants + AFs + carrier counts)")
    print(f"   • synthetic_patients.csv    (all {total_patients:,} patients)")
    for pop in POPULATIONS:
        print(f"   • client_{pop}.csv")
    
    print("\n GENE DISTRIBUTION (top 10):")
    print(df['gene_primary'].value_counts().head(10).to_string())
    
    return df, patients_df


if __name__ == "__main__":
    variants_df, patients_df = main()
