# RAIDers (RAre disease AI and Radar)

RAIDers is a federated computational framework designed to resolve the phenotypic heterogeneity of Amyotrophic Lateral Sclerosis (ALS) while maintaining data sovereignty. By synthesizing global genomic annotations with simulated patient cohorts, this framework establishes a scalable architecture for rare disease subtyping.

<img width="1536" height="1024" alt="F4D4A0A4-4A7B-4710-AA6F-9776780AEBF4" src="https://github.com/user-attachments/assets/c2144ff0-0934-4b54-9f8e-785707df3014" />

---

## Project Overview

Research in rare disease genomics is primarily hindered by extreme data scarcity and institutional data silos mandated by privacy regulations. While authoritative repositories provide critical intelligence on pathogenic variants, these resources are seldom integrated into a unified analytical space.

RAIDers addresses this fragmentation by consolidating disparate genomic signals into a high-fidelity feature matrix. This architecture serves as a "Genomic Flight Simulator," validating a federated subtyping pipeline on synthetic data to demonstrate readiness for integration with controlled-access biobank datasets.

---

## 1. Pipeline Execution (Placeholder)

To initialize the pipeline from synthetic cohort generation through federated clustering:

```bash
# Clone the repository
git clone https://github.com/project/RAIDers.git
cd RAIDers

# Run the primary simulation and analysis (Filename TBD)
python main_pipeline.py

```

*Note: A curated version of `clinvar.cleaned.csv` must be present in the local directory for execution.*

---

## 2. Synthetic Cohort Generation

To overcome the "mathematical invisibility" of rare variants in standard population samples, RAIDers employs a digital mutagenesis strategy. This allows for the generation of a balanced, statistically significant cohort of 15,000 patients partitioned into five ancestral nodes.

### 2.1 Genomic Anchors: ClinVar

The pipeline utilizes `clinvar.cleaned.csv` to identify approximately 450 pathogenic ALS variants. These records provide the biological ground truth for the simulation, including:

* **Gene Association:** (e.g., *SOD1, TARDBP, C9orf72*)
* **Clinical Significance:** Standardized pathogenicity classifications.
* **Molecular Consequence:** Variant-level impact (missense, nonsense, frameshift).

### 2.2 The Interaction Model: gnomAD AF Integration

Rather than using static lookups, we simulate variable penetrance by treating the ancestral background as a clinical modifier. We utilize **gnomAD** as our Genomic Reference Frame. Specifically, we adopt gnomAD’s **superpopulation divisions** and **relative allelic ratios** to construct 'Ancestral Modifiers,' allowing us to simulate how different genomic backgrounds influence clinical expression. We then utilize **gnomAD Allele Frequency (AF)** logic to determine population-specific "tolerance" to pathogenic variants.

**Rationale for AF Estimation:**
Empirical gnomAD frequencies for rare ALS variants are often  or zero in specific subpopulations. Direct application would result in a sparse matrix with insufficient carrier counts for machine learning. We estimate and amplify these frequencies (targeting 0.01% – 0.2%) to ensure analytical viability while maintaining biological realism through Selection penalties (e.g., a 50% AF reduction for Loss-of-Function mutations).

### 2.3 Phenotype Severity Assignment Logic

Clinical labels (e.g., Fast vs. Slow Progression) are derived from the interaction between a variant’s baseline impact and its ancestral modifier (the AF Ratio).

```python
def assign_contextual_phenotype(variant_row, population_id):
    # 1. Mutation Impact (Anchor derived from ClinVar)
    base_impact = 0.8 if "Pathogenic" in variant_row['clinical_sig'] else 0.5
    
    # 2. Ancestral Modifier (gnomAD AF Ratio)
    # High AF ratio implies population tolerance (Protective Modifier)
    # Low AF ratio implies population sensitivity (Aggravating Modifier)
    af_ratio = variant_row[f'gnomAD_AF_{population_id}'] / variant_row['gnomAD_AF']
    modifier = 0.8 if af_ratio > 1.5 else (1.2 if af_ratio < 0.5 else 1.0)
    
    # 3. Probabilistic Interaction (with 5-10% stochastic noise)
    interaction_score = (base_impact * modifier) + np.random.normal(0, 0.05)
    
    return "Fast Progression" if interaction_score > 0.85 else "Slow Progression"

```

### 2.4 Age of Onset Assignment Logic (Placeholder for William)

---

## 3. Federated Analysis: Subtype Discovery

The framework simulates five institutional silos partitioned by superpopulation (AFR, AMR, EAS, EUR, SAS).

### 3.1 Federated Learning Across Simulated Hospitals (Placeholder for Arnav)

Molecular subtypes are discovered through a decentralized K-Means algorithm:

* **Local Iteration:** Clients compute cluster assignments and centroids based on local synthetic cohorts.
* **Global Aggregation:** Centroids are sent to a central server for federated averaging.
* **Broadcast:** Updated global centroids are returned to clients; the process repeats until convergence (change < 0.001).

### 3.2 Analytical Metrics (Placeholder for Arnav)

* **Silhouette Score:** Evaluates the cohesion and separation of discovered molecular clusters.
* **Clustering Stability:** Validated via bootstrap resampling (100 iterations).
* **Biological Validation (Placeholder):** Analysis of cluster-driving features and pathway enrichment.

---

## 4. Scientific Objectives

* **Subtype Discovery:** Identifying coherent molecular signatures across diverse ancestral backgrounds.
* **Feature Validity:** Confirming that integrated annotations contain sufficient signal to distinguish ALS-associated genes.
* **Federated Feasibility:** Demonstrating that analytical fidelity is maintained when data is physically separated across institutional nodes.
* **Ground Truth Validation:** Assessing whether discovered subtypes align with the simulated interaction rules.

---

## Data Sources

* **ClinVar:** Pathogenic variant curation.
* **gnomAD:** Population-level allele frequencies.
* **OMIM/Orphanet:** Clinical gene–disease associations.
* **STRING:** Protein–protein interaction network topology.

# Contributers 

- Aastha Shah
- Arnav Kharbanda
- Bill Paseman
- Chantera Lazard
- Jialan Ma
- Kushal Koirala
- Kyulin Kim
- Nikita Rajesh
- Pu Kao
- Shreya Nandakumar
- Vibha Acharya
- William Lu






