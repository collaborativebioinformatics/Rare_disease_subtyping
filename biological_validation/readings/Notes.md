## ALS Molecular Pathway Definitions (Pathway placeholders)

| Pathway | Description | Therapeutic Relevance | Example Therapeutics | Associated Genes |
|:--------|:------------|:----------------------|:---------------------|:-----------------|
| **RNA Metabolism** | Dysfunction in RNA processing, splicing, and transport; TDP-43/FUS pathology | Antisense oligonucleotides (ASOs), splicing modulators | Tofersen (SOD1 ASO), Jacifusen (FUS ASO), BIIB078 (C9orf72) | TARDBP, FUS, SETX, HNRNPA1, MATR3, ANG, ATXN2, TAF15, EWSR1 |
| **Proteostasis** | Impaired protein folding, aggregation, and degradation via ubiquitin-proteasome system | Heat shock protein inducers, proteasome modulators | Arimoclomol, CuATSM | SOD1, VCP, UBQLN2 |
| **Autophagy-Lysosomal** | Defective autophagy initiation, cargo recognition, or lysosomal clearance | Autophagy inducers, lysosomal modulators, TFEB activators | Rapamycin analogs, Trehalose, Lithium | TBK1, OPTN, SQSTM1, C9orf72, GRN |
| **Axonal Transport** | Disrupted motor protein function and cargo trafficking along axons | Microtubule stabilizers, motor protein modulators | Noscapine derivatives, HDAC6 inhibitors | DCTN1, KIF5A, SPAST |
| **Cytoskeleton** | Destabilized actin/microtubule networks affecting neuronal structure | Cytoskeletal stabilizers, actin modulators | Fasudil (Rho kinase inhibitor) | PFN1, TUBA4A, NEFH |
| **Mitochondrial** | Impaired mitochondrial dynamics, bioenergetics, or oxidative stress response | Antioxidants, mitochondrial biogenesis enhancers | Edaravone (approved), PTC857, Omaveloxolone | CHCHD10, SOD1 |
| **Vesicle Trafficking** | Defective endosomal sorting, ER-Golgi transport, or membrane dynamics | Endosomal trafficking modulators | Emerging targets | ALS2, CHMP2B, FIG4, VAPB, ANXA11 |
| **DNA Damage Response** | Impaired DNA repair mechanisms and cell cycle checkpoint control | PARP inhibitors, DNA repair enhancers | Research stage | NEK1, C21orf2 |
| **Excitotoxicity** | Glutamate-mediated neuronal toxicity and calcium dysregulation | Glutamate modulators, calcium channel blockers | Riluzole (approved), Memantine | DAO |

## Gene-to-Pathway Mapping

| Gene | Primary Pathway | Secondary Pathway | Protein Function | Key Reference |
|:-----|:----------------|:------------------|:-----------------|:--------------|
| TARDBP | RNA Metabolism | - | TDP-43; RNA binding and splicing regulation | Neumann et al., 2006 |
| FUS | RNA Metabolism | - | FUS; RNA binding and transcription regulation | Kwiatkowski et al., 2009 |
| SETX | RNA Metabolism | - | Senataxin; RNA helicase, transcription termination | Chen et al., 2004 |
| HNRNPA1 | RNA Metabolism | - | hnRNP A1; pre-mRNA processing | Kim et al., 2013 |
| MATR3 | RNA Metabolism | - | Matrin 3; RNA processing and retention | Johnson et al., 2014 |
| ANG | RNA Metabolism | - | Angiogenin; tRNA processing, stress response | Greenway et al., 2006 |
| ATXN2 | RNA Metabolism | - | Ataxin-2; RNA metabolism, stress granules | Elden et al., 2010 |
| SOD1 | Proteostasis | Mitochondrial | Superoxide dismutase 1; oxidative stress defense | Rosen et al., 1993 |
| VCP | Proteostasis | Autophagy | Valosin-containing protein; ER protein extraction | Johnson et al., 2010 |
| UBQLN2 | Proteostasis | - | Ubiquilin 2; proteasome shuttling factor | Deng et al., 2011 |
| TBK1 | Autophagy | - | TANK-binding kinase 1; autophagy initiation | Freischmidt et al., 2015 |
| OPTN | Autophagy | - | Optineurin; autophagy receptor | Maruyama et al., 2010 |
| SQSTM1 | Autophagy | - | p62/Sequestosome-1; autophagy adaptor | Fecto et al., 2011 |
| C9orf72 | Autophagy | RNA Metabolism | C9orf72; autophagy regulation, RNA foci | DeJesus-Hernandez et al., 2011 |
| GRN | Autophagy | Inflammation | Progranulin; lysosomal function | Baker et al., 2006 |
| DCTN1 | Axonal Transport | - | Dynactin subunit 1; retrograde transport | Puls et al., 2003 |
| KIF5A | Axonal Transport | - | Kinesin heavy chain; anterograde transport | Nicolas et al., 2018 |
| SPAST | Axonal Transport | - | Spastin; microtubule severing | Meyer et al., 2005 |
| PFN1 | Cytoskeleton | - | Profilin 1; actin polymerization | Wu et al., 2012 |
| TUBA4A | Cytoskeleton | - | Tubulin alpha-4A; microtubule structure | Smith et al., 2014 |
| NEFH | Cytoskeleton | - | Neurofilament heavy chain; axonal caliber | Figlewicz et al., 1994 |
| CHCHD10 | Mitochondrial | - | CHCHD10; mitochondrial cristae maintenance | Bannwarth et al., 2014 |
| ALS2 | Vesicle Trafficking | - | Alsin; endosomal trafficking, Rab5 GEF | Yang et al., 2001 |
| CHMP2B | Vesicle Trafficking | - | CHMP2B; ESCRT-III, endosomal sorting | Skibinski et al., 2005 |
| FIG4 | Vesicle Trafficking | - | FIG4; phosphoinositide phosphatase | Chow et al., 2009 |
| VAPB | Vesicle Trafficking | - | VAPB; ER-Golgi trafficking, UPR | Nishimura et al., 2004 |
| ANXA11 | Vesicle Trafficking | - | Annexin A11; RNA granule transport | Smith et al., 2017 |
| NEK1 | DNA Damage Response | - | NEK1; DNA damage checkpoint, cilia | Kenna et al., 2016 |
| C21orf2 | DNA Damage Response | - | C21orf2; ciliogenesis, DNA damage | van Rheenen et al., 2016 |
| DAO | Excitotoxicity | - | D-amino acid oxidase; D-serine metabolism | Mitchell et al., 2010 |

## Pathway Feature Encoding for Clustering

| Pathway Code | Pathway Name | One-Hot Column Name | Genes (n) |
|:-------------|:-------------|:--------------------|:----------|
| P1 | RNA Metabolism | `pathway_RNA_Metabolism` | 8 |
| P2 | Proteostasis | `pathway_Proteostasis` | 3 |
| P3 | Autophagy-Lysosomal | `pathway_Autophagy` | 5 |
| P4 | Axonal Transport | `pathway_Axonal_Transport` | 3 |
| P5 | Cytoskeleton | `pathway_Cytoskeleton` | 3 |
| P6 | Mitochondrial | `pathway_Mitochondrial` | 2 |
| P7 | Vesicle Trafficking | `pathway_Vesicle_Trafficking` | 5 |
| P8 | DNA Damage Response | `pathway_DNA_Damage` | 2 |
| P9 | Excitotoxicity | `pathway_Excitotoxicity` | 1 |

## Multi-Pathway Gene Assignments

'C9orf72': ['Autophagy', 'RNA_Metabolism']  # Both pathways

# One-Hot Encoding for Pathway Features

## Overview

One-hot encoding converts categorical pathway labels into binary numerical columns that machine learning algorithms (like k-means clustering) can process. This document explains the rationale, implementation, and expected outcomes for pathway-based subtyping.

---

## The Problem: Why One-Hot Encoding?

Clustering algorithms operate on **numerical vectors**, not categorical strings. Raw pathway labels cannot be used directly:
```python
# ❌ This doesn't work - algorithms can't interpret strings
patient = {
    'patient_id': 'AFR_00001',
    'primary_pathway': 'RNA_Metabolism'  # String - not usable for clustering
}
```

---

## The Solution: Binary Pathway Columns

Create a separate column for each pathway, with values of `1` (pathway affected) or `0` (pathway not affected):
```python
# ✅ This works - each pathway is a numerical feature
patient = {
    'patient_id': 'AFR_00001',
    'pathway_RNA_Metabolism': 1,
    'pathway_Autophagy': 0,
    'pathway_Proteostasis': 0,
    'pathway_Mitochondrial': 0,
    # ... other pathways
}
```

---

## Visual Example

### Before One-Hot Encoding

| patient_id | gene | primary_pathway |
|:-----------|:-----|:----------------|
| AFR_00001 | TARDBP | RNA_Metabolism |
| AFR_00002 | TBK1 | Autophagy |
| AFR_00003 | FUS | RNA_Metabolism |
| AFR_00004 | SOD1 | Proteostasis |
| AFR_00005 | C9orf72 | Autophagy |

### After One-Hot Encoding

| patient_id | gene | pathway_RNA_Metabolism | pathway_Autophagy | pathway_Proteostasis | pathway_Mitochondrial |
|:-----------|:-----|:----------------------:|:-----------------:|:--------------------:|:---------------------:|
| AFR_00001 | TARDBP | 1 | 0 | 0 | 0 |
| AFR_00002 | TBK1 | 0 | 1 | 0 | 0 |
| AFR_00003 | FUS | 1 | 0 | 0 | 0 |
| AFR_00004 | SOD1 | 0 | 0 | 1 | 1 |
| AFR_00005 | C9orf72 | 1 | 1 | 0 | 0 |

> **Note:** SOD1 and C9orf72 have multiple `1` values because they affect multiple pathways.

---

## Handling Multi-Pathway Genes

Some genes affect multiple cellular pathways. One-hot encoding handles this elegantly by setting multiple columns to `1`:

### Multi-Pathway Gene Examples

| Gene | Pathways Affected | pathway_RNA | pathway_Autophagy | pathway_Proteostasis | pathway_Mito |
|:-----|:------------------|:-----------:|:-----------------:|:--------------------:|:------------:|
| TARDBP | RNA Metabolism only | 1 | 0 | 0 | 0 |
| TBK1 | Autophagy only | 0 | 1 | 0 | 0 |
| C9orf72 | RNA + Autophagy | 1 | 1 | 0 | 0 |
| SOD1 | Proteostasis + Mitochondrial | 0 | 0 | 1 | 1 |
| VCP | Proteostasis + Autophagy | 0 | 1 | 1 | 0 |

---

## Implementation

### Step 1: Define Pathway Constants
```python
# =============================================================================
# PATHWAY DEFINITIONS
# =============================================================================

PATHWAYS = [
    'RNA_Metabolism',
    'Proteostasis',
    'Autophagy',
    'Axonal_Transport',
    'Cytoskeleton',
    'Mitochondrial',
    'Vesicle_Trafficking',
    'DNA_Damage',
    'Excitotoxicity',
    'Unknown'
]
```

### Step 2: Create Gene-to-Pathway Mapping
```python
# =============================================================================
# GENE TO PATHWAY MAPPING
# =============================================================================

GENE_TO_PATHWAY = {
    # RNA Metabolism
    'TARDBP': ['RNA_Metabolism'],
    'FUS': ['RNA_Metabolism'],
    'SETX': ['RNA_Metabolism'],
    'HNRNPA1': ['RNA_Metabolism'],
    'MATR3': ['RNA_Metabolism'],
    'ANG': ['RNA_Metabolism'],
    'ATXN2': ['RNA_Metabolism'],
    
    # Proteostasis (some with secondary pathways)
    'SOD1': ['Proteostasis', 'Mitochondrial'],
    'VCP': ['Proteostasis', 'Autophagy'],
    'UBQLN2': ['Proteostasis'],
    
    # Autophagy-Lysosomal
    'TBK1': ['Autophagy'],
    'OPTN': ['Autophagy'],
    'SQSTM1': ['Autophagy'],
    'C9orf72': ['Autophagy', 'RNA_Metabolism'],
    'GRN': ['Autophagy'],
    
    # Axonal Transport
    'DCTN1': ['Axonal_Transport'],
    'KIF5A': ['Axonal_Transport'],
    'SPAST': ['Axonal_Transport'],
    
    # Cytoskeleton
    'PFN1': ['Cytoskeleton'],
    'TUBA4A': ['Cytoskeleton'],
    'NEFH': ['Cytoskeleton'],
    
    # Mitochondrial
    'CHCHD10': ['Mitochondrial'],
    
    # Vesicle Trafficking
    'ALS2': ['Vesicle_Trafficking'],
    'CHMP2B': ['Vesicle_Trafficking'],
    'FIG4': ['Vesicle_Trafficking'],
    'VAPB': ['Vesicle_Trafficking'],
    'ANXA11': ['Vesicle_Trafficking'],
    
    # DNA Damage Response
    'NEK1': ['DNA_Damage'],
    'C21orf2': ['DNA_Damage'],
    
    # Excitotoxicity
    'DAO': ['Excitotoxicity'],
}
```

### Step 3: One-Hot Encoding Function
```python
def get_pathway_one_hot(gene_name):
    """
    Convert gene name to one-hot encoded pathway dictionary.
    
    Parameters
    ----------
    gene_name : str
        Gene symbol (e.g., 'SOD1', 'C9orf72')
    
    Returns
    -------
    dict
        Dictionary with pathway columns as keys and 0/1 as values
    
    Examples
    --------
    >>> get_pathway_one_hot('TARDBP')
    {
        'pathway_RNA_Metabolism': 1,
        'pathway_Proteostasis': 0,
        'pathway_Autophagy': 0,
        'pathway_Axonal_Transport': 0,
        'pathway_Cytoskeleton': 0,
        'pathway_Mitochondrial': 0,
        'pathway_Vesicle_Trafficking': 0,
        'pathway_DNA_Damage': 0,
        'pathway_Excitotoxicity': 0,
        'pathway_Unknown': 0
    }
    
    >>> get_pathway_one_hot('C9orf72')
    {
        'pathway_RNA_Metabolism': 1,    # Multi-pathway gene
        'pathway_Proteostasis': 0,
        'pathway_Autophagy': 1,          # Multi-pathway gene
        'pathway_Axonal_Transport': 0,
        'pathway_Cytoskeleton': 0,
        'pathway_Mitochondrial': 0,
        'pathway_Vesicle_Trafficking': 0,
        'pathway_DNA_Damage': 0,
        'pathway_Excitotoxicity': 0,
        'pathway_Unknown': 0
    }
    """
    # Initialize all pathways to 0
    one_hot = {f'pathway_{pathway}': 0 for pathway in PATHWAYS}
    
    # Get pathways for this gene (default to 'Unknown' if not found)
    gene_clean = str(gene_name).upper().strip()
    gene_pathways = GENE_TO_PATHWAY.get(gene_clean, ['Unknown'])
    
    # Set matching pathways to 1
    for pathway in gene_pathways:
        column_name = f'pathway_{pathway}'
        if column_name in one_hot:
            one_hot[column_name] = 1
    
    return one_hot
```

### Step 4: Get Primary Pathway Label
```python
def get_primary_pathway(gene_name):
    """
    Get the primary (first-listed) pathway for a gene.
    
    Parameters
    ----------
    gene_name : str
        Gene symbol (e.g., 'SOD1', 'C9orf72')
    
    Returns
    -------
    str
        Primary pathway name
    
    Examples
    --------
    >>> get_primary_pathway('SOD1')
    'Proteostasis'
    
    >>> get_primary_pathway('C9orf72')
    'Autophagy'
    
    >>> get_primary_pathway('UNKNOWN_GENE')
    'Unknown'
    """
    gene_clean = str(gene_name).upper().strip()
    pathways = GENE_TO_PATHWAY.get(gene_clean, ['Unknown'])
    return pathways[0]
```

### Step 5: Integration with Patient Record Creation
```python
# Inside the patient generation loop:

# Create base patient record
patient = {
    'patient_id': f"{pop}_{i:05d}",
    'superpopulation': pop,
    'primary_gene': primary['gene'],
    'primary_condition': primary['condition'],
    'primary_subtype': primary['als_subtype'],
    
    # Severity (kept as outcome variable, NOT used in clustering)
    'severity_score': primary['severity_score'],
    'severity_category': primary['severity_category'],
    
    # Interaction model outputs
    'interaction_score': primary['interaction_score'],
    'predicted_progression': primary['contextual_progression'],
}

# Add primary pathway label
patient['primary_pathway'] = get_primary_pathway(primary['gene'])

# Add one-hot encoded pathway features (USED in clustering)
pathway_features = get_pathway_one_hot(primary['gene'])
patient.update(pathway_features)

# Final patient record now includes:
# {
#     'patient_id': 'AFR_00001',
#     'primary_gene': 'C9orf72',
#     'primary_pathway': 'Autophagy',
#     'severity_score': 8.0,                    # NOT for clustering
#     'pathway_RNA_Metabolism': 1,              # FOR clustering
#     'pathway_Autophagy': 1,                   # FOR clustering
#     'pathway_Proteostasis': 0,                # FOR clustering
#     'pathway_Mitochondrial': 0,               # FOR clustering
#     ...
# }
```

---

## Clustering Feature Selection

### Before (Severity-Driven Clustering)
```python
# Old approach - severity dominates
clustering_features = [
    'severity_score',        # Highly correlated with gene
    'interaction_score',     # Derived from severity
    'geno_var1', 'geno_var2', ...  # Sparse genotype columns
]
```

### After (Pathway-Driven Clustering)
```python
# New approach - pathway is explicit
clustering_features = [
    # Pathway features (primary signal)
    'pathway_RNA_Metabolism',
    'pathway_Proteostasis',
    'pathway_Autophagy',
    'pathway_Axonal_Transport',
    'pathway_Cytoskeleton',
    'pathway_Mitochondrial',
    'pathway_Vesicle_Trafficking',
    'pathway_DNA_Damage',
    'pathway_Excitotoxicity',
    
    # Ancestral modifier (contextual signal)
    'ancestral_modifier',
    
    # Genotype columns (variant signal)
    'geno_var1', 'geno_var2', ...
]

# Severity is EXCLUDED from clustering but kept for analysis
outcome_variables = [
    'severity_score',
    'severity_category',
    'predicted_progression',
]
```

---

## Expected Clustering Outcomes

### Before: Severity-Stratified Clusters

| Cluster | Dominant Feature | Severity | Pathway Mix |
|:--------|:-----------------|:---------|:------------|
| C0 | Low severity | 4.0 - 5.0 | Mixed (SETX, ALS2, ...) |
| C1 | Moderate severity | 6.0 - 7.0 | Mixed (SOD1, OPTN, ...) |
| C2 | High severity | 8.0 - 9.0 | Mixed (FUS, TBK1, ...) |

> **Problem:** Patients with same pathway but different severity are split across clusters.

### After: Pathway-Based Clusters

| Cluster | Dominant Pathway | Severity Range | Key Genes |
|:--------|:-----------------|:---------------|:----------|
| C0 | RNA Metabolism | 4.0 - 9.0 | TARDBP, FUS, SETX |
| C1 | Autophagy | 5.0 - 8.0 | TBK1, OPTN, C9orf72 |
| C2 | Proteostasis | 6.0 - 9.0 | SOD1, VCP, UBQLN2 |
| C3 | Vesicle Trafficking | 4.0 - 7.0 | ALS2, CHMP2B, FIG4 |

> **Benefit:** Severity varies **within** each pathway cluster — true molecular subtypes.

---

## Therapeutic Relevance

The pathway-based clustering directly maps to therapeutic strategies:

| Cluster | Pathway | Therapeutic Approach | Example Drugs |
|:--------|:--------|:---------------------|:--------------|
| C0 | RNA Metabolism | ASOs, splicing modulators | Tofersen, Jacifusen |
| C1 | Autophagy | Autophagy inducers | Rapamycin, Trehalose |
| C2 | Proteostasis | HSP inducers, proteasome modulators | Arimoclomol, CuATSM |
| C3 | Vesicle Trafficking | Trafficking modulators | Emerging targets |

This is the core value proposition: **cluster assignment predicts therapeutic strategy**, regardless of severity.

---

## Complete Working Example
```python
#!/usr/bin/env python3
"""
One-Hot Encoding Example for ALS Pathway Features
"""

# =============================================================================
# CONSTANTS
# =============================================================================

PATHWAYS = [
    'RNA_Metabolism',
    'Proteostasis',
    'Autophagy',
    'Axonal_Transport',
    'Cytoskeleton',
    'Mitochondrial',
    'Vesicle_Trafficking',
    'DNA_Damage',
    'Excitotoxicity',
    'Unknown'
]

GENE_TO_PATHWAY = {
    'TARDBP': ['RNA_Metabolism'],
    'FUS': ['RNA_Metabolism'],
    'SOD1': ['Proteostasis', 'Mitochondrial'],
    'C9orf72': ['Autophagy', 'RNA_Metabolism'],
    'TBK1': ['Autophagy'],
    'ALS2': ['Vesicle_Trafficking'],
    'CHCHD10': ['Mitochondrial'],
    'NEK1': ['DNA_Damage'],
}

# =============================================================================
# FUNCTIONS
# =============================================================================

def get_pathway_one_hot(gene_name):
    """Convert gene to one-hot encoded pathway dict."""
    one_hot = {f'pathway_{p}': 0 for p in PATHWAYS}
    gene_clean = str(gene_name).upper().strip()
    gene_pathways = GENE_TO_PATHWAY.get(gene_clean, ['Unknown'])
    
    for pathway in gene_pathways:
        col = f'pathway_{pathway}'
        if col in one_hot:
            one_hot[col] = 1
    
    return one_hot


def get_primary_pathway(gene_name):
    """Get primary pathway for a gene."""
    gene_clean = str(gene_name).upper().strip()
    return GENE_TO_PATHWAY.get(gene_clean, ['Unknown'])[0]


# =============================================================================
# DEMONSTRATION
# =============================================================================

if __name__ == "__main__":
    test_genes = ['TARDBP', 'SOD1', 'C9orf72', 'TBK1', 'UNKNOWN_GENE']
    
    print("Gene-to-Pathway One-Hot Encoding Examples")
    print("=" * 60)
    
    for gene in test_genes:
        print(f"\nGene: {gene}")
        print(f"  Primary Pathway: {get_primary_pathway(gene)}")
        print(f"  One-Hot Encoding:")
        
        one_hot = get_pathway_one_hot(gene)
        active_pathways = [k for k, v in one_hot.items() if v == 1]
        
        for pathway in active_pathways:
            print(f"    {pathway}: 1")
```

### Example Output
```
Gene-to-Pathway One-Hot Encoding Examples
============================================================

Gene: TARDBP
  Primary Pathway: RNA_Metabolism
  One-Hot Encoding:
    pathway_RNA_Metabolism: 1

Gene: SOD1
  Primary Pathway: Proteostasis
  One-Hot Encoding:
    pathway_Proteostasis: 1
    pathway_Mitochondrial: 1

Gene: C9orf72
  Primary Pathway: Autophagy
  One-Hot Encoding:
    pathway_Autophagy: 1
    pathway_RNA_Metabolism: 1

Gene: TBK1
  Primary Pathway: Autophagy
  One-Hot Encoding:
    pathway_Autophagy: 1

Gene: UNKNOWN_GENE
  Primary Pathway: Unknown
  One-Hot Encoding:
    pathway_Unknown: 1
```

---

## Summary

| Aspect | Before | After |
|:-------|:-------|:------|
| Pathway representation | Implicit (via gene → severity) | Explicit (one-hot columns) |
| Clustering driver | Severity score | Pathway features |
| Multi-pathway genes | Not handled | Multiple columns = 1 |
| Therapeutic relevance | Indirect | Direct mapping |
| Subtype definition | Severity-stratified | Molecularly-grounded |

# Then in patient record:
patient['pathway_Autophagy'] = 1
patient['pathway_RNA_Metabolism'] = 1  # Both = 1
```
- Biologically accurate
- Allows clustering to find "dual-pathway" patients

ALSoD (ALS Online Database)-Curated gene list with annotations- https://alsod.ac.uk/
OMIM  Gene-disease relationships, pathways  -  https://omim.org/
GeneCards   Pathway annotations per gene     https://www.genecards.org/
KEGG Pathways Formal pathway definitions https://www.kegg.jp/
Gene Ontology (GO) Biological process annotations http://geneontology.org/ 
UniProt   Protein function, pathway involvement    https://www.uniprot.org/
