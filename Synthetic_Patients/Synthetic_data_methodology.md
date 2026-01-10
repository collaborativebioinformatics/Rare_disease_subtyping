# Technical Methodology: Synthetic Biobank Generation and Interaction Modeling

This document provides a comprehensive technical elaboration of the **RAIDers** (RAre disease AI and Radar) pipeline. It details the algorithmic transformation of raw genomic annotations into a high-fidelity synthetic patient cohort of 15,000 individuals partitioned across five ancestral nodes.

---

### 1. Data Ingestion and Variant Standardization

The initial phase focus on converting raw, heterogeneous variant data from `clinvar.cleaned.csv` into a structured format suitable for population-scale simulation.

* **Metadata Integration**: The script preserves essential clinical metadata such as `clinical_sig` (Pathogenicity) and `consequence` (e.g., missense, nonsense).
* **Biological Anchors**: These attributes serve as the fundamental biological anchors for subsequent simulations, ensuring every patient record is linked to a documented ALS mutation.
* **Standardization**: High-fidelity cleaning involves creating unique variant IDs, normalizing genomic coordinates, and extracting primary gene symbols to ensure accurate mapping to the gene-to-phenotype knowledge base.

---

### 2. The Rationale for gnomAD Allele Frequency Estimation

The engine implements `estimate_gnomad_af()` to generate population-specific Allele Frequencies () for five superpopulations: **AFR, AMR, EAS, EUR, and SAS**. This step is foundational for establishing a "Genomic Ground Truth" for the federated model.

#### A. Biological Logic Modeling (Impact and Rarity)

The calculation is driven by evolutionary principles:

* **Pathogenicity Constraints**: Variants labeled as "Pathogenic" are assigned lower base frequencies (0.005% − 0.15%), reflecting evolutionary pressure against harmful mutations.
* **Functional Penalties**: "Loss of Function" (LoF) consequences undergo a 0.5 × frequency penalty to simulate higher levels of purifying selection. Missense variants are assigned a 1.2 × multiplier as they are often better tolerated.

**Implementation Snippet:**

```python
# A. Pathogenicity-based base frequency
if str(clinical_sig).strip() == "Pathogenic":
    base_af = np.random.uniform(0.00005, 0.0015)
else:
    base_af = np.random.uniform(0.0001, 0.002)

# B. Functional impact penalty for LoF variants
lof_consequences = ['nonsense', 'frameshift', 'splice donor', 'splice acceptor', 'start lost']
if any(lof in str(consequence).lower() for lof in lof_consequences):
    base_af *= 0.5  # 50% reduction for LoF

```

#### B. Enabling the "Ancestral Modifier" Interaction

By generating variation between populations (e.g., 0.7−1.6 for AFR vs. 0.4−1.3 for EAS), the pipeline creates the necessary signal for the AI to learn **Contextualized Penetrance**.

---

### 3. Stochastic Genotype Simulation: The HWE Model

The script translates estimated population frequencies into individual-level genetic data using the **Hardy-Weinberg Equilibrium (HWE)** principle.

**Mathematical Framework:**
The probability distribution for individual genotypes $G \in \{0, 1, 2 } \$ given the alternate allele frequency $p$and reference frequency $q = 1 - p$:

The probability distribution for individual genotypes $G \in \{0, 1, 2\}$ given allele frequency $p$ and $q = 1 - p$:

<img width="501" height="136" alt="Screenshot 2026-01-09 at 1 36 58 pm" src="https://github.com/user-attachments/assets/38cb7c0a-0bd9-4819-9c91-63ec3ea30d38" />​

For each of the **3,000 samples per population**, the script performs a weighted random choice based on these normalized probabilities to assign a specific genotype for every variant.

---

### 4. Rare Variant Preservation Logic

In rare disease modeling, stochastic sampling often results in "lost" variants where zero carriers are generated.

* **Carrier Assurance**: The `ensure_minimum_carriers` function identifies variants resulting in zero carriers after HWE simulation.
  
* **Manual Injection**: To guarantee analytical viability, the script selectively assigns heterozygous genotypes () to a minimum number of patients (defined by `MIN_CARRIERS_PER_VARIANT`). This ensures every pathogenic variant is available for downstream subtyping.

---

### 5. Contextual Interaction Model: Calculating Penetrance

The core innovation is the mapping of interaction scores to clinical labels based on ancestral context.

#### A. Relative Allelic Ratio ()

This ratio determines contextual sensitivity relative to the global mean:

$$R = \frac{AF_{\text{Population}}}{AF_{\text{Global}}}$$

#### B. Ancestral Modifier ()

These ratios determine the clinical modifier assigned to a patient:

* **PROTECTIVE (M=0.8)**: Assigned when R > 1.5 (High frequency suggests evolved tolerance).

* **AGGRAVATING (M=1.2)**: Assigned when R < 0.5 (Extreme rarity suggests high sensitivity).

* **NEUTRAL (M=1.0)**: Assigned when 0.5 ≤ R ≤ 1.5.

#### C. Final Interaction Score ()

The score synthesizes the baseline impact with ancestral and molecular modifiers:

S = ( I × M × C) + ϵ

Where:

I is the Mutation Anchor (Base Impact).

M is the Ancestral Modifier.

C is the Consequence Multiplier.

ϵ ∼ N(0,0.05) is the stochastic noise factor.

Where  represents stochastic noise.

**Phenotype Mapping Thresholds:**
| Final Score ($S$) | Clinical Label |
| :--- | :--- |
| $S > 0.85$ | Fast Progression |
| $0.60 < S \le 0.85$ | Slow Progression |
| $S \le 0.60$ | Asymptomatic / Low Penetrance |


---

### 6. Severity Modeling and Categorization

The framework calculates a numeric severity score by adjusting a gene's baseline severity against its molecular consequence and clinical significance.

**6.1 Mathematical Framework**

The calculation follows a multiplicative logic to represent the compounding effect of genetic disruption:

$$
S_{\text{sev}} = \min \left( 10.0,\ \max \left( 1.0,\( B \times C \times P \right)) \right)
$$

Where:

- **B (Base Severity):** Curated baseline severity assigned per gene  
  - Example: FUS = 9.0, SOD1 = 7.0

- **C (Consequence Multiplier):** Weight based on molecular consequence  
  - Frameshift: 1.4  
  - Nonsense: 1.3  
  - Missense: 1.0  

- **P (Pathogenicity Bonus):**  
  - P = 1.1 if `clinical_sig = "Pathogenic"`  
  - P = 1.0 otherwise

**Implementation Snippet:**

```python
def calculate_severity_score(base_severity, consequence, clinical_sig):
    score = base_severity
    score *= get_consequence_multiplier(consequence)
    if str(clinical_sig).strip().lower() == 'pathogenic':
        score *= 1.1
    return min(10.0, max(1.0, round(score, 1)))

```

**6.2. Clinical Categorization Logic**

To facilitate easier interpretation for downstream federated subtyping, the continuous S sev is discretized into four distinct clinical categories.


| Score Range ($S_{\text{sev}}$) | Clinical Category | Rationale |
|-------------------------------|------------------|-----------|
| $S_{\text{sev}} \geq 8.0$ | Severe | High-impact mutations in aggressive gene backgrounds (e.g., FUS nonsense). |
| $6.0 \leq S_{\text{sev}} < 8.0$ | Moderate | Typical presentation of classical ALS genes (e.g., SOD1 or TARDBP). |
| $4.0 \leq S_{\text{sev}} < 6.0$ | Mild | Juvenile-onset or slower progression genes (e.g., ALS2 or SETX). |
| $S_{\text{sev}} < 4.0$ | Low | Risk modifiers or variants with very high ancestral tolerance. |
---

### 7. Validation and Statistical Integrity

* **AF vs. Carrier Correlation**: A Pearson correlation coefficient is calculated between the estimated  and actual carrier counts to validate simulation accuracy.
* **Filtered Exporting**: The pipeline generates `patients_carriers_only.csv`, providing a focused environment for high-dimensional subtyping analysis.
