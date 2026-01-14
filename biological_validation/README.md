## Overview 

This folder contains biological validation analyses of federated k-means clustering results on synthetic ALS patient data. Our analysis reveals **severity-stratified molecular subtypes** or patient clusters that combine disease severity with distinct pathway dysfunction patterns. 

## Key Findings
We found four clinically-informed molecular clusters. 
| Cluster | Size | Severity | Progression | Dominant Pathways | Key Gene(s) |
| :--: | :--: | :--: | :--: | :--: | :--: |
| C0 | 287 | Mild (5.0) | Slow (74%) | RNA metabolism, Proteostasis | SETX, ALS2 |
| C4 | 845 | Moderate (6.3) | Moderate (91%) | Multipathway | ALS2, OPTN |
| C1 | 3,287 | Moderately Severe (7.3) | Moderate (93%) | Multipathway (SOD1-dominant) | SOD1 |
| C2 | 1,999 | Severe (9.3) | Fast (96.5%) | Autophagy-Transport collapse | TBK1, FUS |
| C3 | 8,582 | Control (0.0) | - | Healthy Controls | - |

## Major Discoveries 
1. **Severity Gradient**: Clear progression from mild (C0, 5.0) → moderate (C4, 6.3) → moderately-severe (C1, 7.3) → severe (C2, 9.3)
2. **Progression Velocity**: Different pathways predicts different progression rates
3. **Molecular Subtypes Within Severity**: Clusters 1 and 4 share the same severity category, but have different distinct molecular profiles.
4. **Variant Heterogeneity as Biomarker**:
   - Low std (C0: 0.03) = homogeneous variants -> predictable outcomes
   - High std (C1: 0.10) = heterogeneous variants -> variable phenotypes -> needs sub-clustering
  
