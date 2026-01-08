# RAIDers



---

# 1. Project Overview
Rare disease genomics faces two major challenges:
1. Data scarcity
2. Data silos across biobanks due to privacy constraints.

Rare disease genomics research is limited by small sample sizes, heterogeneous annotations, and fragmented data across institutions. Although many public databases curate valuable information on disease-causing variants, these resources are rarely analysed together in a unified and scalable way.

This project focuses on Amyotrophic Lateral Sclerosis (ALS) and aims to integrate real, publicly available genomic and biological annotations from multiple databases into a single feature space for exploratory analysis

<img width="2456" height="1842" alt="RAIDer flow-3" src="https://github.com/user-attachments/assets/2c14a579-d1b7-4e6a-9009-99d4ca71ec2b" />



# 2. Scientific Goal


The primary scientific goal of this project is to determine whether molecular subtypes of ALS-associated genes and variants can be identified using integrated real-world annotations, and whether these subtypes reflect biological signal rather than methodological or population-driven artifacts.

Specifically, we aim to answer the following questions:
1. Subtype Discovery:
   
Can unsupervised clustering of integrated variant-level features reveal coherent molecular subtypes among ALS-associated genes?

2. Feature Validity:

Do the extracted features contain sufficient biological information to accurately distinguish ALS genes in a supervised classification setting?

3. Methodological Consistency:
   
Are the subtypes discovered via unsupervised learning supported by supervised model performance and feature importance analysis?

4. Federated Feasibility

Can clustering and validation be performed in a federated manner that mirrors institutional data separation while preserving analytical fidelity?

By addressing these questions, the project evaluates whether federated, annotation-driven analysis can support meaningful rare disease stratification using only real, publicly available data, providing a foundation for future extension to controlled-access biobank datasets.


# 3. Data Sources 

Each participating site has access to a subset of ALS-related genes and variants and integrates information from:
- ClinVar – pathogenic / likely pathogenic ALS variants
- gnomAD – population allele frequencies
- OMIM – gene–disease associations
- STRING – protein–protein interaction networks


# Local Data Processing 

All preprocessing is performed locally at each client before any information is shared.

1. Text Embedding Generation
- Free-text fields (e.g. variant descriptions, disease annotations, phenotype summaries)
- Embedded locally using PubMedBERT
- Output: dense semantic embedding vectors

2. Numerical Feature Normalization

Examples include:
- Population allele frequencies
- Network-derived scores
- Variant-level quantitative annotations

Standard normalization or scaling is applied locally.

3. Categorical Feature Encoding

Examples include:
	•	Variant type
	•	Gene identifier
	•	Inheritance pattern

Encoded using:
- One-hot encoding or
- Label encoding (depending on feature cardinality)

4. Feature Vector Construction


For each record:
- Text columns are embedded using PubMedBERT
- Numerical features are normalized (e.g. allele frequencies, network scores)
- Categorical features are encoded (one-hot or label encoding)
- All features are concatenated into a single fixed-length feature vector

The resulting feature matrix is used for downstream analysis.



# 2. Unsupervised Analysis: Subtype Discovery

Goal:

To discover molecular subtypes of ALS-associated genes and variants via clustering.

# 2.1 Vertical Federated Clustering Implementation

Federated learning is introduced at the algorithmic level, not at raw data ingestion.

Client Definition
- Each client represents an institution
- Each client is a public database

Federated K-Means
1. Each client computes:
- Local cluster assignments
- Local centroids based on its subset of data
2. Clients send locally computed embeddings to the server
3. The server aggregates centroids via federated averaging
4. Updated global centroids are broadcast back to clients

This process is repeated until convergence.

Convergence criterion:
Change in centroid position < 0.001


# 2.2 Unsupervised Validation
- Internal metrics
- Silhouette score
- Davies–Bouldin index
- Biological validation
- Gene enrichment analysis
- Pathway and functional coherence
- Stability analysis
- Bootstrap resampling (100 iterations)
- Cluster consistency across resamples


# 3. Supervised Analysis: Feature Validation

Goal

To validate that extracted features are biologically informative, rather than clustering artifacts.

# 3.1 Classification Task
- Target: Gene identity
(SOD1, C9orf72, FUS, TARDBP, others)
- Input: Same feature set used in clustering
- Train/test split: 80/20, stratified by gene


# 3.2 Models
- XGBoost classifier
- TabTransformer (comparison model)


# 3.3 Supervised Validation
- Metrics
- Accuracy
- Weighted F1-score
- Confusion matrix
- Feature importance
- SHAP value analysis
- Cross-validation
- 5-fold stratified cross-validation


# 4. Integration & Cross-Validation

To connect unsupervised discovery with supervised validation:
- Compare cluster assignments with gene prediction performance
- Assess whether gene pairs with high supervised separability form distinct clusters
- Use supervised feature importance to interpret cluster-driving features
- Evaluate consistency between molecular subtypes and gene-level predictability



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






