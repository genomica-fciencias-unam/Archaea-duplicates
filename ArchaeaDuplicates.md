```python
# Archaeal Duplication
# Stage 1: Sequence Acquisition and Curation

#This section details the first phase of the Archaeal Duplication pipeline. It covers the complete sequence acquisition workflow, raw assembly standardization, de novo protein prediction, and metadata integration.

```


```python
# Step 1 Genome Acquisition via RefSeq (rsync)

#We retrieved 1,829 complete archaeal genomes from the NCBI RefSeq database via `rsync`.
```


```bash
%%bash
# Download complete archaeal genomic FNA files from NCBI RefSeq
rsync -avz rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/*/reference/*/*_genomic.fna.gz .

# Download corresponding CDS sequences
rsync -avz rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/*/reference/*/*_cds_from_genomic.fna.gz .
```


```python
# Step 2 Sequence Header Standardization (`header_fasta.pl`)

#To ensure traceability across the pipeline, we reformatted all FASTA headers using a custom Perl script (`bincomunista/header_fasta.pl`) to embed the Genome Assembly Identifier (GCF/GCA) directly into each sequence header.
```


```bash
%%bash
#!/bin/bash
# Description: Uncompress, apply header_fasta.pl, and recompress FASTA files

input_dir="/home/isramv/Archaea_Duplicados"
perl_script="bincomunista/header_fasta.pl"

for file in $input_dir/*.fna.gz; do
    if [[ -f "$file" ]]; then
        base=$(basename "$file" .fna.gz)
        echo "Processing $base..."

        # Decompress to temporary file
        zcat "$file" > "$input_dir/$base.fna"

        # Apply custom Perl script to embed Genome ID in headers
        perl "$perl_script" seq "$input_dir/$base.fna"

        # Re-compress output file
        gzip -f "$input_dir/$base.fna"
    fi
done

echo "Header processing completed."
```


```python
# Step 3 De novo Protein Prediction (Prodigal)

#To bypass heterogeneous annotation sources and standardize gene calling across all genomes, we performed de novo protein prediction using Prodigal v2.6.3
```


```bash
%%bash
#!/bin/bash
# Description: Batch Prodigal prediction on standardized Archaeal genomes

INPUT_DIR="/home/isramv/Archaea_Duplicados/Genomas_ref_Archaea/GenomasCompletos"
OUTPUT_DIR="${INPUT_DIR}/prodigal_results"
mkdir -p "$OUTPUT_DIR"

for fasta in "$INPUT_DIR"/*.fna; do
    if [[ -f "$fasta" ]]; then
        base=$(basename "$fasta" .fna)
        echo "Processing $base with Prodigal..."
        
        prodigal -i "$fasta" \
                 -a "${OUTPUT_DIR}/${base}_proteins.faa" \
                 -d "${OUTPUT_DIR}/${base}_genes.fna" \
                 -o "${OUTPUT_DIR}/${base}_genes.gbk" \
                 -p single \
                 -g 11 \
                 -q
        
        if [[ $? -eq 0 ]]; then
            echo " Completed $base: proteins in ${base}_proteins.faa"
        else
            echo "✗ Error in $base"
        fi
    fi
done

echo "Prodigal prediction completed. Results in: $OUTPUT_DIR"
```


```python
# Step 4 Genome Statistics Generation (EMBOSS infoseq)

#We generated genome-level statistics, including total sequence length (in base pairs) and protein counts, using EMBOSS `infoseq` and custom Bash utilities.
```


```bash
%%bash
#!/bin/bash
# Description: Compute genome length per assembly using EMBOSS infoseq

OUTPUT_FILE="longitudes_totales_$(date +%Y%m%d_%H%M%S).txt"
> "$OUTPUT_FILE"

echo "Generating assembly length table..."

for archivo in *_genomic.fna; do
    if [[ -f "$archivo" ]]; then
        # Extract GCF Accession Identifier
        id_archivo=$(basename "$archivo" | sed 's/\..*//')
        
        # Calculate total bases via infoseq
        total_bases=$(infoseq -only -length "$archivo" | awk '{sum += $NF} END {print sum}')
        
        echo -e "$id_archivo\t$total_bases" >> "$OUTPUT_FILE"
    fi
done

echo "File generated: $OUTPUT_FILE"
head -n 5 "$OUTPUT_FILE"
```


```python
# Step 5 Metadata Integration & Taxonomic Curation in R (`taxize` / `taxizedb`)

#We standardized taxonomic metadata using the NCBI Taxonomy database. We integrated genome statistics (`infoseq` output) with taxonomy in R using the `taxize` and `taxizedb` packages. 

#We systematically retrieved full lineages for each GCF identifier (Phylum to Genus) and applied recursive assignment logic for incomplete nomenclature or *incertae sedis* to prevent data loss.
```


```python
# R Script: Taxonomic curation
library(taxize)
library(taxizedb)
library(dplyr)
library(readr)

# 1. Load genome assembly statistics
infoseq_data <- read_tsv("longitudes_totales_20260301_100000.txt", 
                         col_names = c("GCF_ID", "Genome_Length"))

ncbi_meta <- read_tsv("Tabla_NCBI_Taxa.tsv")

# Merge length statistics with initial NCBI table
genome_metadata <- left_join(infoseq_data, ncbi_meta, by = c("GCF_ID" = "Assembly_Accession"))

# 2. Retrieve full taxonomic lineages using taxizedb
cat("Retrieving full NCBI lineages for GCF identifiers...\n")

# Extract unique TaxIDs or species names
taxids <- genome_metadata$TaxID

# Fetch full classification lineage
lineages <- taxizedb::classification(taxids, db = "ncbi")

# 3. Recursive logic parser for missing or incertae sedis ranks
parse_lineage <- function(lineage_df) {
  if (is.null(lineage_df) || is.na(lineage_df)) {
    return(c(phylum=NA, class=NA, order=NA, family=NA, genus=NA))
  }
  
  ranks <- c("phylum", "class", "order", "family", "genus")
  parsed <- setNames(rep(NA_character_, length(ranks)), ranks)
  
  last_valid <- "Unclassified_Archaea"
  
  for (r in ranks) {
    match_row <- lineage_df %>% filter(rank == r)
    
    if (nrow(match_row) > 0 && !is.na(match_row$name) && match_row$name != "" && !grepl("incertae sedis", match_row$name, ignore.case = TRUE)) {
      parsed[r] <- match_row$name
      last_valid <- match_row$name
    } else {
      # Fallback to the immediately higher valid rank to preserve classification hierarchy
      parsed[r] <- paste0("Unclassified_", last_valid)
    }
  }
  return(parsed)
}

# Apply lineage parsing across dataset
tax_matrix <- do.call(rbind, lapply(lineages, parse_lineage)) %>% as.data.frame()

# Combine with final master metadata frame
final_metadata <- cbind(genome_metadata, tax_matrix)

# 4. Save standardized master taxonomy table
write_tsv(final_metadata, "Tabla_InfoseqCompleta.tsv")
```


```python
# Section Over
```


```python
# Archaeal Duplication
# Section 2 Quantification of gene duplication and clustering

#This section details the identification of gene family expansions and the quantification of gene duplication across Archaeal genomes.
```


```python
# Step 1: Global Clustering (MMseqs2)

#We aggregate all predicted proteomes (`*_proteins.faa`) into a single master FASTA file and execute MMseqs2 across sequence identity cutoffs of **40%**, **70%**, and **90%** with a **70%** coverage filter.
```


```bash
%%bash
#!/bin/bash
# Description: Global clustering (MMseqs2)

PROT_DIR="/home/isramv/Archaea_Duplicados/Genomas_ref_Archaea/GenomasCompletos/prodigal_results"
GLOBAL_DIR="./global_clustering_mmseqs2"
MASTER_FASTA="${GLOBAL_DIR}/all_archaeal_proteins.faa"
TMP_DIR="${GLOBAL_DIR}/tmp"

mkdir -p "$GLOBAL_DIR" "$TMP_DIR"

# 1. Pool all predicted proteomes into a single aggregate dataset
echo "Concatenating all predicted proteomes into a master FASTA..."
cat "$PROT_DIR"/*_proteins.faa > "$MASTER_FASTA"

# Build MMseqs2 database
mmseqs createdb "$MASTER_FASTA" "${GLOBAL_DIR}/archaeal_global_db"

# 2. Iterate across sequence identity thresholds (40%, 70%, 90%)
THRESHOLDS=(0.40 0.70 0.90)

for id in "${THRESHOLDS[@]}"; do
    echo "=================================================="
    echo "Running Global MMseqs2 Clustering at ${id} min identity..."
    echo "=================================================="
    
    OUT_CLUSTER="${GLOBAL_DIR}/global_cluster_id${id}"
    OUT_TSV="${GLOBAL_DIR}/global_clusters_id${id}.tsv"
    
    # Cluster with 70% coverage filter (-c 0.7) to ensure structural homology
    mmseqs cluster "${GLOBAL_DIR}/archaeal_global_db" "$OUT_CLUSTER" "$TMP_DIR" \
        --min-seq-id "$id" \
        -c 0.7 \
        --cov-mode 0 \
        --threads 16
    
    # Generate flat TSV (Representative ID \t Member ID)
    mmseqs createtsv "${GLOBAL_DIR}/archaeal_global_db" "${GLOBAL_DIR}/archaeal_global_db" "$OUT_CLUSTER" "$OUT_TSV"
    
    echo " Global clustering for threshold ${id} finished -> ${OUT_TSV}"
done

rm -rf "$TMP_DIR"
echo "Global clustering complete."
```


```python
# Step 2 Individual Clustering (MMseqs2)

#To quantify intra-genome paralogy without the noise of inter-genomic orthology, we execute an independent all-versus-all self-clustering for each of the 1,829 archaeal proteomes.
```


```bash
%%bash
#!/bin/bash
# Description: Batch intra-genomic self-clustering for each proteome individually

PROT_DIR="/home/isramv/Archaea_Duplicados/Genomas_ref_Archaea/GenomasCompletos/prodigal_results"
INTRA_DIR="./intra_genomic_clustering"
TMP_DIR="${INTRA_DIR}/tmp"

mkdir -p "$INTRA_DIR" "$TMP_DIR"

THRESHOLDS=(0.40 0.70 0.90)

for faa in "$PROT_DIR"/*_proteins.faa; do
    if [[ -f "$faa" ]]; then
        gcf_id=$(basename "$faa" | sed 's/_proteins\.faa//')
        echo "Processing intra-genomic self-clustering for assembly: $gcf_id"
        
        genome_db="${TMP_DIR}/${gcf_id}_db"
        mmseqs createdb "$faa" "$genome_db" &>/dev/null
        
        for id in "${THRESHOLDS[@]}"; do
            out_cluster="${TMP_DIR}/${gcf_id}_cluster_id${id}"
            out_tsv="${INTRA_DIR}/${gcf_id}_self_cluster_id${id}.tsv"
            
            mmseqs cluster "$genome_db" "$out_cluster" "$TMP_DIR" \
                --min-seq-id "$id" \
                -c 0.7 \
                --cov-mode 0 \
                --threads 4 &>/dev/null
            
            mmseqs createtsv "$genome_db" "$genome_db" "$out_cluster" "$out_tsv" &>/dev/null
        done
        
        # Clean up database files for this genome
        rm -f ${TMP_DIR}/${gcf_id}_db* ${TMP_DIR}/${gcf_id}_cluster*
    fi
done

rm -rf "$TMP_DIR"
echo "Intra-genomic clustering completed for all genomes."
```


```python
# Step 3 Genomic Metrics(Python)

#We parse the clustering outputs to calculate three primary metrics per genome
```


```python
import os
import glob
import pandas as pd

def calculate_paralogy_metrics(cluster_tsv):
    """
    Parses MMseqs2 clustering TSV (Col1: Representative, Col2: Member)
    and computes duplication metrics per genome.
    """
    df = pd.read_csv(cluster_tsv, sep="\t", header=None, names=["Representative", "Member"])
    
    # Calculate cluster sizes (N)
    cluster_sizes = df.groupby("Representative").size().reset_index(name="Cluster_Size")
    
    # Merge size back to member level
    df = df.merge(cluster_sizes, on="Representative")
    
    # Extract Genome Accession (GCF) from Member IDs
    df["Genome_ID"] = df["Member"].apply(lambda x: "_".join(x.split("_")[:2]))
    
    metrics = []
    for gcf, group in df.groupby("Genome_ID"):
        total_proteome = len(group)
        
        # Singletons: clusters of size 1
        singletons = len(group[group["Cluster_Size"] == 1])
        singleton_frequency = singletons / total_proteome if total_proteome > 0 else 0
        
        # Duplicated proteins: members belonging to clusters of size > 1
        duplicated_proteins = len(group[group["Cluster_Size"] > 1])
        duplicated_proteome_fraction = (duplicated_proteins / total_proteome) * 100 if total_proteome > 0 else 0
        
        # Total Duplicate Count: Sum of (N - 1) per cluster family in the genome
        family_dups = group.groupby("Representative")["Cluster_Size"].first() - 1
        total_duplicate_count = family_dups[family_dups > 0].sum()
        
        metrics.append({
            "Genome_ID": gcf,
            "Total_Proteome": total_proteome,
            "Total_Duplicate_Count": total_duplicate_count,
            "Singletons": singletons,
            "Singleton_Frequency": round(singleton_frequency, 4),
            "Duplicated_Proteins": duplicated_proteins,
            "Duplicated_Proteome_Fraction_Pct": round(duplicated_proteome_fraction, 2)
        })
        
    return pd.DataFrame(metrics)

# Run metric computation across global identity thresholds
global_tsvs = sorted(glob.glob("./global_clustering_mmseqs2/global_clusters_id*.tsv"))

for tsv in global_tsvs:
    id_cutoff = tsv.split("_id")[-1].replace(".tsv", "")
    print(f"Calculating duplication metrics for Global Threshold: {id_cutoff}")
    
    metrics_df = calculate_paralogy_metrics(tsv)
    out_path = f"./global_clustering_mmseqs2/duplication_metrics_global_id{id_cutoff}.tsv"
    metrics_df.to_csv(out_path, sep="\t", index=False)
    print(f"  -> Saved: {out_path}")
```


```python
# Section Over
```


```python
# Archaeal Duplication
# Section 3 Abundance matrix construction and filtering
#This section details the transformation of raw clustering outputs into high-dimensional count matrices for each sequence identity threshold (40%, 70%, and 90%). 
```


```python
# Step 1 High-Dimensional Count Matrix Construction (`Python`)

#We parse the global MMseqs2 cluster output TSVs (40%, 70%, and 90% identity) to extract protein-to-genome mappings and pivot the data into wide-format abundance matrices.
```


```python
import os
import glob
import pandas as pd
import numpy as np

def build_raw_abundance_matrix(cluster_tsv):
    """
    Transforms flat MMseqs2 TSV into a Cluster x Genome count matrix.
    """
    print(f"Reading raw cluster file: {cluster_tsv}")
    # Read TSV (Col 1: Cluster Representative, Col 2: Member ID)
    df = pd.read_csv(cluster_tsv, sep="\t", header=None, names=["Cluster_ID", "Member_ID"])
    
    # Extract Genome Assembly Identifier (GCF/GCA prefix)
    df["Genome_ID"] = df["Member_ID"].apply(lambda x: "_".join(x.split("_")[:2]))
    
    # Group by Cluster_ID and Genome_ID to calculate copy numbers (paralogs per family per genome)
    print("  -> Pivot-table transformation: counting paralogs per cluster per genome...")
    count_df = df.groupby(["Cluster_ID", "Genome_ID"]).size().reset_index(name="Copy_Number")
    
    # Pivot to wide format matrix (Rows: Clusters, Cols: Genomes)
    matrix = count_df.pivot(index="Cluster_ID", columns="Genome_ID", values="Copy_Number").fillna(0).astype(int)
    
    return matrix

# Output directory for matrices
matrix_dir = "./abundance_matrices"
os.makedirs(matrix_dir, exist_ok=True)

# Process all identity thresholds
global_tsvs = sorted(glob.glob("./global_clustering_mmseqs2/global_clusters_id*.tsv"))

for tsv in global_tsvs:
    id_cutoff = tsv.split("_id")[-1].replace(".tsv", "")
    print(f"\n==================================================")
    print(f"Building Raw Abundance Matrix for Threshold ID: {id_cutoff}")
    print(f"==================================================")
    
    raw_matrix = build_raw_abundance_matrix(tsv)
    
    out_raw = f"{matrix_dir}/raw_abundance_matrix_id{id_cutoff}.tsv"
    raw_matrix.to_csv(out_raw, sep="\t")
    
    print(f" Raw Matrix generated: {raw_matrix.shape[0]} Clusters x {raw_matrix.shape[1]} Genomes")
    print(f"  -> Saved to: {out_raw}")
```


```python
# Step 2: Filtering Protocol for Duplication

#To focus downstream statistical modeling on gene duplication dynamics and paralogous expansion, we filter out non-duplicated lineage-specific singletons.
```


```python
def filter_abundance_matrix(matrix):
    """
    Applies filtering protocol to prioritize clusters with evidence of gene duplication
    or inter-genomic conservation.
    """
    initial_clusters = matrix.shape[0]
    
    # 1. Identify clusters with duplication (Copy Number >= 2 in at least 1 genome)
    max_copy_per_cluster = matrix.max(axis=1)
    duplicated_mask = max_copy_per_cluster >= 2
    
    # 2. Identify clusters with inter-genomic conservation (Present in >= 2 genomes)
    genome_occupancy = (matrix > 0).sum(axis=1)
    conserved_mask = genome_occupancy >= 2
    
    # Combined filtering mask: Keep if it has duplications OR inter-genomic conservation
    keep_mask = duplicated_mask | conserved_mask
    
    filtered_matrix = matrix[keep_mask]
    
    print(f"  -> Initial Clusters: {initial_clusters}")
    print(f"  -> Clusters with intra-genome duplication (max copy >= 2): {duplicated_mask.sum()}")
    print(f"  -> Clusters with inter-genome conservation (present in >= 2 genomes): {conserved_mask.sum()}")
    print(f"  -> Final Filtered Clusters Retained: {filtered_matrix.shape[0]} (Removed {initial_clusters - filtered_matrix.shape[0]} uninformative singletons)")
    
    return filtered_matrix

# Execute filtering across all raw generated matrices
raw_matrix_files = sorted(glob.glob(f"{matrix_dir}/raw_abundance_matrix_id*.tsv"))

for matrix_file in raw_matrix_files:
    id_cutoff = matrix_file.split("_id")[-1].replace(".tsv", "")
    print(f"\n--------------------------------------------------")
    print(f"Filtering Abundance Matrix for Threshold ID: {id_cutoff}")
    print(f"--------------------------------------------------")
    
    # Load raw matrix
    raw_mat = pd.read_csv(matrix_file, sep="\t", index_col=0)
    
    # Apply filtering protocol
    filtered_mat = filter_abundance_matrix(raw_mat)
    
    # Export filtered matrix
    out_filtered = f"{matrix_dir}/filtered_abundance_matrix_id{id_cutoff}.tsv"
    filtered_mat.to_csv(out_filtered, sep="\t")
    print(f" Filtered Matrix exported to: {out_filtered}")
```


```python
# Step 3: Matrix Validation and Quality Check (R)

#We load the filtered abundance matrices into R to verify dimensions and generate diagnostic summary statistics prior to Z-score normalization.
```


```python
# R Script: Matrix validation and structural inspection
library(dplyr)
library(readr)
library(tibble)

matrix_dir <- "./abundance_matrices"
thresholds <- c("0.40", "0.70", "0.90")

for (id in thresholds) {
  file_path <- file.path(matrix_dir, paste0("filtered_abundance_matrix_id", id, ".tsv"))
  
  if (file.exists(file_path)) {
    cat("\n==================================================\n")
    cat("Validating Filtered Matrix ID:", id, "\n")
    cat("==================================================\n")
    
    # Load matrix
    mat <- read_tsv(file_path, show_col_types = FALSE) %>% column_to_rownames(var = "Cluster_ID")
    
    cat("Dimensions:", nrow(mat), "Clusters x", ncol(mat), "Genomes\n")
    cat("Total Gene Copies in Matrix:", sum(mat), "\n")
    cat("Max Copy Number observed in a single cell:", max(mat), "\n")
    cat("Sparsity (% zeros):", round(sum(mat == 0) / (nrow(mat) * ncol(mat)) * 100, 2), "%\n")
  }
}
```


```python
# Section Over
```


```python
# Archaeal Duplication 
# Section 4 Statistical modeling and outlier detection

#This section details the linear regression modeling between total duplicated proteins and total proteome size across the 1,829 Archaeal genomes. 
```


```python
# Step 1 Linear Regression & Standardized Residual Calculation (R)

#We load the duplicated proteome metrics alongside genome metadata, fit linear models for each identity threshold (40%, 70%, 90%), calculate leverage-adjusted standardized residuals using `rstandard()`, and identify outlier lineages.
```


```python
# R Script: Linear regression and leverage-adjusted Z-score modeling
library(tidyverse)
library(broom)

# Load master metadata and duplication metrics
metadata <- read_tsv("Tabla_InfoseqCompleta.tsv", show_col_types = FALSE)
metrics_dir <- "./global_clustering_mmseqs2"
thresholds <- c("0.40", "0.70", "0.90")

# Results container
all_residuals_list <- list()

for (id in thresholds) {
  cat("\n==================================================\n")
  cat("Modeling Duplication Dynamics for Threshold ID:", id, "\n")
  cat("==================================================\n")
  
  metrics_file <- file.path(metrics_dir, paste0("duplication_metrics_global_id", id, ".tsv"))
  
  if (!file.exists(metrics_file)) {
    warning(paste("File not found:", metrics_file))
    next
  }
  
  # Load metrics for current threshold
  metrics_df <- read_tsv(metrics_file, show_col_types = FALSE)
  
  # Merge with taxonomy metadata
  model_data <- inner_join(metrics_df, metadata, by = c("Genome_ID" = "GCF_ID"))
  
  # Fit linear regression: Total Duplicated Proteins ~ Total Proteome Size
  lm_model <- lm(Duplicated_Proteins ~ Total_Proteome, data = model_data)
  
  cat("Linear Model Summary (ID", id, "):\n")
  print(summary(lm_model))
  
  # Calculate leverage-adjusted standardized residuals via rstandard()
  # Standardized residuals follow N(0, 1) and function directly as Z-scores
  model_data$Standardized_Residual_Z <- rstandard(lm_model)
  model_data$Raw_Residual <- residuals(lm_model)
  model_data$Leverage <- hatvalues(lm_model)
  model_data$Identity_Threshold <- id
  
  # Flag statistically significant outliers (|Z| > 2.0)
  model_data <- model_data %>% 
    mutate(Is_Outlier = abs(Standardized_Residual_Z) > 2.0,
           Outlier_Type = case_when(
             Standardized_Residual_Z > 2.0 ~ "Expansion_Outlier",
             Standardized_Residual_Z < -2.0 ~ "Depletion_Outlier",
             TRUE ~ "Background"
           ))
  
  num_expansion_outliers <- sum(model_data$Outlier_Type == "Expansion_Outlier")
  cat("  -> Identified", num_expansion_outliers, "significant expansion outliers (Z > +2.0)\n")
  
  all_residuals_list[[id]] <- model_data
}

# Combine multi-scale results
multi_scale_df <- bind_rows(all_residuals_list)
write_tsv(multi_scale_df, "regression_multi_scale_residuals.tsv")
```


```python
# Step 2 Longitudinal Trajectory Classification (Python)

#We parse the multi-scale $Z$-scores across the 40%, 70%, and 90% identity thresholds
```


```python
import pandas as pd
import numpy as np

# Load multi-scale residual output from R
df = pd.read_csv("regression_multi_scale_residuals.tsv", sep="\t")

# Pivot matrix to get Z-scores per genome across thresholds (Columns: 0.40, 0.70, 0.90)
z_matrix = df.pivot(index="Genome_ID", columns="Identity_Threshold", values="Standardized_Residual_Z")
z_matrix.columns = [f"Z_{col}" for col in z_matrix.columns]

def classify_trajectory(row):
    """
    Classifies evolutionary trajectories based on multi-scale Z-score patterns.
    """
    z40 = row["Z_0.40"]
    z70 = row["Z_0.70"]
    z90 = row["Z_0.90"]
    
    # 1. Permanent Duplication: Consistent Z > 2.0 across all scales
    if z40 >= 2.0 and z70 >= 2.0 and z90 >= 2.0:
        return "Permanent Duplication"
    
    # 2. Emerging Expansions: Increasing Z towards 90% (Z90 > Z40 and Z90 >= 2.0)
    elif z90 >= 2.0 and z90 > z40:
        return "Emerging Expansions"
    
    # 3. Historical Divergences: High at 40% (Z40 >= 2.0) that collapse at 70%/90% (Z90 < 1.0 or Z70/Z90 << Z40)
    elif z40 >= 2.0 and (z90 < 1.0 or (z40 - z90) >= 1.5):
        return "Historical Divergence"
    
    # 4. Cooling or Diverging Lineages: Steady monotonic decline in residuals (Z40 > Z70 > Z90)
    elif z40 > z70 and z70 > z90 and (z40 - z90) >= 0.5:
        return "Cooling Lineage"
    
    else:
        return "Stable/Background"

# Apply classification logic
z_matrix["Evolutionary_Trajectory"] = z_matrix.apply(classify_trajectory, axis=1)

# Merge back taxonomy metadata for reporting
meta_cols = ["Genome_ID", "phylum", "class", "order", "family", "genus"]
metadata_unique = df[meta_cols].drop_duplicates(subset=["Genome_ID"])
final_trajectories = z_matrix.reset_index().merge(metadata_unique, on="Genome_ID")

# Output path
final_trajectories.to_csv("archaeal_duplication_trajectories.tsv", sep="\t", index=False)

# Summary counts per trajectory
print("\n==================================================")
print("EVOLUTIONARY TRAJECTORY CLASSIFICATION SUMMARY")
print("==================================================")
summary_counts = final_trajectories["Evolutionary_Trajectory"].value_counts()
print(summary_counts.to_string())
```


```python
# Step 3 Regression and Trajectory Diagnostic Plots (R)

#Generate diagnostic visualizations of the linear models, Z-score distributions, and multi-scale trajectory trends.
```


```python
# R Script: Visualization of regression models and trajectory shifts
library(ggplot2)
library(dplyr)
library(readr)

trajectories_df <- read_tsv("archaeal_duplication_trajectories.tsv", show_col_types = FALSE)
residuals_df <- read_tsv("regression_multi_scale_residuals.tsv", show_col_types = FALSE)

# 1. Plot Linear Regression (Total Duplicated Proteins vs. Total Proteome Size)
p_reg <- ggplot(residuals_df, aes(x = Total_Proteome, y = Duplicated_Proteins)) +
  geom_point(aes(color = Outlier_Type), alpha = 0.6, size = 1.8) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  facet_wrap(~ Identity_Threshold, labeller = label_both) +
  scale_color_manual(values = c("Expansion_Outlier" = "#D55E00", 
                                "Depletion_Outlier" = "#0072B2", 
                                "Background" = "gray60")) +
  theme_minimal() +
  labs(title = "Duplication Scaling Relative to Proteome Size",
       x = "Total Proteome Size (Proteins)",
       y = "Total Duplicated Proteins",
       color = "Status")

ggsave("regression_scaling_plots.png", plot = p_reg, width = 10, height = 4, dpi = 300)

# 2. Plot Multi-Scale Trajectory Shift Profiles
trajectory_long <- trajectories_df %>%
  filter(Evolutionary_Trajectory != "Stable/Background") %>%
  pivot_longer(cols = starts_with("Z_"), names_to = "Threshold", values_to = "Z_Score") %>%
  mutate(Threshold = gsub("Z_", "", Threshold))

p_traj <- ggplot(trajectory_long, aes(x = Threshold, y = Z_Score, group = Genome_ID, color = Evolutionary_Trajectory)) +
  geom_line(alpha = 0.4) +
  geom_hline(yintercept = 2.0, linetype = "dashed", color = "red") +
  facet_wrap(~ Evolutionary_Trajectory) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  labs(title = "Evolutionary Trajectories of Paralog-Enriched Lineages",
       x = "Sequence Identity Threshold",
       y = "Standardized Residual (Z-score)",
       color = "Trajectory")

ggsave("trajectory_shift_profiles.png", plot = p_traj, width = 9, height = 6, dpi = 300)
```


```python
# Section Over
```


```python
# Archaeal Duplication 
# Section 5 Functional annotation

#This section details the multi-layered functional annotation of representative protein sequences using **eggNOG-mapper v2.1.9** against the eggNOG v5.0 database. It covers similarity searches via **DIAMOND v2**, mobilome filtering via a curated keyword lexicon, and high-level KEGG pathway/module mapping using the R package **KEGGREST**.
```


```python
# Step 1: Representative Protein Sequence Extraction & eggNOG-mapper Execution

#We extract representative protein sequences for all filtered clusters and run **eggNOG-mapper v2.1.9** using **DIAMOND v2** in `--sens` (sensitive) mode to ensure remote homolog detection.
```


```bash
%%bash
#!/bin/bash
# Description: Extract representative sequences and perform eggNOG-mapper functional annotation

CLUSTER_DIR="./global_clustering_mmseqs2"
OUTPUT_DIR="./functional_annotation"
EGGNOG_DB="/db/eggnog_v5.0" # Adjust to your local eggNOG database path

mkdir -p "$OUTPUT_DIR"

# 1. Extract cluster representative protein FASTA using MMseqs2
echo "Extracting representative protein sequences for annotation..."
mmseqs createsubdb "${CLUSTER_DIR}/global_cluster_id0.40" "${CLUSTER_DIR}/archaeal_global_db" "${OUTPUT_DIR}/rep_seqs_db"
mmseqs convert2fasta "${OUTPUT_DIR}/rep_seqs_db" "${OUTPUT_DIR}/representative_proteins_id0.40.faa"

# 2. Run eggNOG-mapper v2.1.9 with DIAMOND v2 in sensitive mode
echo "Running eggNOG-mapper v2.1.9 (DIAMOND sensitive mode)..."
emapper.py \
    -i "${OUTPUT_DIR}/representative_proteins_id0.40.faa" \
    --output_dir "$OUTPUT_DIR" \
    -o "eggnog_results_id0.40" \
    -m diamond \
    --sens \
    --cpu 16 \
    --data_dir "$EGGNOG_DB" \
    --override

echo "eggNOG-mapper completed. Output generated in $OUTPUT_DIR"
```


```python
# Step 2 Mobilome Lexicon Filtering and Functional Categorization (Python)

#To isolate core informational and metabolic duplications from mobile genetic element (MGE) activity, we implement a regex-based keyword classification protocol. Proteins matching mobilome terms are categorized under **"Mobile Elements"**.
```


```python
import pandas as pd
import re

# File paths
emapper_tsv = "./functional_annotation/eggnog_results_id0.40.emapper.annotations"
output_curated = "./functional_annotation/functional_annotations_curated.tsv"

# Read eggNOG-mapper output (skipping comment lines starting with #)
print("Loading eggNOG-mapper output...")
df = pd.read_csv(emapper_tsv, sep="\t", comment="#", header=None)

# Standard eggNOG-mapper v2.1.9 header columns
emapper_cols = [
    "query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl",
    "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko",
    "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE",
    "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs"
]
df.columns = emapper_cols[:len(df.columns)]

# Define mobilome regex lexicon
mge_pattern = re.compile(
    r"\b(transposase|integrase|IS[0-9]+|DDE|HTH|resolvase|recombinase|insertion sequence|invader|prophage|conjugative|plasmid)\b",
    re.IGNORECASE
)

def classify_function(row):
    description = str(row["Description"])
    cog_cat = str(row["COG_category"])
    
    # 1. Check for Mobile Genetic Element (MGE) keywords in Description or COG 'X'
    if mge_pattern.search(description) or "X" in cog_cat:
        return "Mobile Elements"
    
    # 2. Return primary COG category or Unclassified
    elif cog_cat != "-" and pd.notna(cog_cat):
        return cog_cat
    else:
        return "Unclassified"

print("Applying mobilome keyword-based classification protocol...")
df["Functional_Category"] = df.apply(classify_function, axis=1)

# Summary of MGE vs Core categorization
mge_count = (df["Functional_Category"] == "Mobile Elements").sum()
total_annotated = len(df)
print(f"  -> Total Representatives Annotated: {total_annotated}")
print(f"  -> Mobilome/MGE Cluster Families Identified: {mge_count} ({round(mge_count/total_annotated * 100, 2)}%)")

# Save curated functional annotations
df.to_csv(output_curated, sep="\t", index=False)
print(f"Curated functional dataset saved to {output_curated}")
```


```python
# Step 3 Higher-Level KEGG Mapping via `KEGGREST` (R)

#We integrate the curated functional annotations into an R pipeline using the `KEGGREST` package to map KEGG Orthology (KO) identifiers to higher-level biological structures, including metabolic pathways and functional modules.
```


```python
# R Script KEGGREST integration for higher-level pathway and module mapping
library(KEGGREST)
library(tidyverse)

# 1. Load curated functional annotation dataset
annot_path <- "./functional_annotation/functional_annotations_curated.tsv"
annot_df <- read_tsv(annot_path, show_col_types = FALSE)

# Filter for valid KEGG KO entries
ko_df <- annot_df %>%
  filter(!is.na(KEGG_ko) & KEGG_ko != "-") %>%
  separate_rows(KEGG_ko, sep = ",") %>%
  mutate(KEGG_ko = gsub("ko:", "", KEGG_ko))

unique_kos <- unique(ko_df$KEGG_ko)
cat("Found", length(unique_kos), "unique KO identifiers for KEGGREST mapping.\n")

# 2. Query KEGG REST API in batches to retrieve Pathway and Module names
get_kegg_mappings <- function(ko_list) {
  pathway_map <- list()
  module_map <- list()
  
  # Process in chunks of 50 to respect API rate limits
  chunks <- split(ko_list, ceiling(seq_along(ko_list) / 50))
  
  for (i in seq_along(chunks)) {
    chunk <- chunks[[i]]
    cat("Querying KEGGREST batch", i, "of", length(chunks), "...\n")
    
    tryCatch({
      # Retrieve KO entry details
      query_res <- keggGet(chunk)
      
      for (entry in query_res) {
        current_ko <- entry$ENTRY
        
        # Extract Pathways
        if (!is.null(entry$PATHWAY)) {
          pathway_map[[current_ko]] <- paste(names(entry$PATHWAY), entry$PATHWAY, collapse = "; ")
        } else {
          pathway_map[[current_ko]] <- NA
        }
        
        # Extract Modules
        if (!is.null(entry$MODULE)) {
          module_map[[current_ko]] <- paste(names(entry$MODULE), entry$MODULE, collapse = "; ")
        } else {
          module_map[[current_ko]] <- NA
        }
      }
    }, error = function(e) {
      warning(paste("Error querying batch", i, ":", e$message))
    })
    
    Sys.sleep(0.2) # API rate limit control
  }
  
  # Build mapping data frame
  mapping_df <- tibble(
    KEGG_ko = names(pathway_map),
    KEGG_Pathway_Name = unlist(pathway_map),
    KEGG_Module_Name = unlist(module_map)
  )
  return(mapping_df)
}

# Run KEGGREST Mapping
kegg_hierarchy_df <- get_kegg_mappings(unique_kos)

# 3. Merge higher-level KEGG metadata back to main functional annotation frame
final_functional_dataset <- left_join(ko_df, kegg_hierarchy_df, by = "KEGG_ko")

# Save final integrated table
write_tsv(final_functional_dataset, "./functional_annotation/integrated_kegg_functional_dataset.tsv")
```


```python
# Section Over
```


```python
# Archaeal Duplication
# Section 6 Transposase database construction and distribution analysis

#This section details the construction of a specialized transposase reference database and the assessment of mobilome distribution across prokaryotic domains. Seed sequences from TnCentral/ISFinder (10,657 proteins) are clustered into 2,086 family representatives, followed by exploratory domain comparisons (iTOL, Bray-Curtis PCoA, PERMANOVA) and a comprehensive similarity search against NCBI RefSeq (81.6 million sequences) using DIAMOND v2.0.15.
```


```python
# Step 1 Transposase Seed Clustering (MMseqs2)

#We cluster the 10,657 seed transposase sequences retrieved from TnCentral/ISFinder using MMseqs2 at **40% sequence identity** and **70% coverage** (`-c 0.7`) to extract 2,086 representative sequences.
```


```bash
%%bash
#!/bin/bash
# Description: Cluster TnCentral seed transposases to 2,086 family representatives

RAW_TN_FASTA="./transposase_db/TnCentral_raw_seeds.fasta"
OUT_DIR="./transposase_db"
TMP_DIR="${OUT_DIR}/tmp"

mkdir -p "$OUT_DIR" "$TMP_DIR"

echo "Creating MMseqs2 database from raw TnCentral seeds..."
mmseqs createdb "$RAW_TN_FASTA" "${OUT_DIR}/tn_seeds_db"

echo "Clustering transposases at 40% identity and 70% query coverage..."
mmseqs cluster "${OUT_DIR}/tn_seeds_db" "${OUT_DIR}/tn_cluster_id0.40" "$TMP_DIR" \
    --min-seq-id 0.40 \
    -c 0.7 \
    --cov-mode 0 \
    --threads 16

# Extract representative sequences (2,086 families)
mmseqs createsubdb "${OUT_DIR}/tn_cluster_id0.40" "${OUT_DIR}/tn_seeds_db" "${OUT_DIR}/tn_rep_db"
mmseqs convert2fasta "${OUT_DIR}/tn_rep_db" "${OUT_DIR}/tn_representatives_2086.faa"

REP_COUNT=$(grep -c "^>" "${OUT_DIR}/tn_representatives_2086.faa")
echo "Extraction complete: ${REP_COUNT} transposase family representatives generated."

rm -rf "$TMP_DIR"
```


```python
# Step 2 Exploratory Domain Comparison, Bray-Curtis PCoA & PERMANOVA (R)

#We compare transposase family relative abundances across Archaea and a subset of Bacteria, perform Principal Coordinate Analysis (PCoA) on Bray-Curtis distances, and test for statistical significance using PERMANOVA (999 permutations).
```


```python
# R Script: Bray-Curtis PCoA and PERMANOVA for Transposase Distribution
library(vegan)
library(tidyverse)

# 1. Load Transposase Abundance Matrix (Rows: Genomes, Cols: Transposase Families)
abundance_matrix <- read_tsv("./transposase_db/domain_tnp_abundance_matrix.tsv", show_col_types = FALSE)

# Extract metadata (Domain classification) and count matrix
sample_info <- abundance_matrix %>% select(Genome_ID, Domain)
counts_mat <- abundance_matrix %>% select(-Genome_ID, -Domain) %>% as.matrix()
rownames(counts_mat) <- sample_info$Genome_ID

# 2. Transform raw counts to Relative Abundances
rel_abundance <- decostand(counts_mat, method = "total")

# 3. Calculate Bray-Curtis Dissimilarity Distance Matrix
cat("Computing Bray-Curtis dissimilarity matrix...\n")
bray_dist <- vegdist(rel_abundance, method = "bray")

# 4. Perform Principal Coordinate Analysis (PCoA / CmdScale)
pcoa_res <- cmdscale(bray_dist, k = 2, eig = TRUE)
pcoa_df <- data.frame(
  Genome_ID = rownames(pcoa_res$points),
  PCoA1 = pcoa_res$points[, 1],
  PCoA2 = pcoa_res$points[, 2]
) %>% left_join(sample_info, by = "Genome_ID")

# Calculate Variance Explained
var_exp <- round((pcoa_res$eig[1:2] / sum(pcoa_res$eig)) * 100, 2)

# Save PCoA coordinates for iTOL / R plotting
write_tsv(pcoa_df, "./transposase_db/pcoa_bray_curtis_coords.tsv")

# 5. PERMANOVA Test (999 Permutations) via adonis2
cat("\n==================================================\n")
cat("PERMANOVA TEST (Domain Structure Significance)\n")
cat("==================================================\n")
permanova_res <- adonis2(bray_dist ~ Domain, data = sample_info, permutations = 999)
print(permanova_res)

# Export PERMANOVA summary
write_tsv(as.data.frame(permanova_res), "./transposase_db/permanova_results.tsv")
```


```python
# Step 3 Global Similarity Search Against NCBI RefSeq (DIAMOND v2.0.15)

#We perform similarity searches using the 2,086 representative transposase sequences against the custom-indexed complete NCBI RefSeq database (81,617,152 sequences; 31.2 billion residues) using **DIAMOND v2.0.15** in `blastp` mode, applying an **$E$-value threshold of $10^{-10}$** and a **minimum query and subject coverage of 70%**.
```


```bash
%%bash
#!/bin/bash
# Description: High-throughput DIAMOND blastp search against complete RefSeq

QUERY_FAA="./transposase_db/tn_representatives_2086.faa"
REFSEQ_DB="/db/refseq_2025/refseq_protein.dmnd" # Custom DIAMOND index (81,617,152 seqs)
OUT_TSV="./transposase_db/diamond_refseq_tnp_hits.tsv"

echo "Executing DIAMOND blastp search against RefSeq (81.6M sequences)..."

diamond blastp \
    --query "$QUERY_FAA" \
    --db "$REFSEQ_DB" \
    --out "$OUT_TSV" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp \
    --evalue 1e-10 \
    --query-cover 70 \
    --threads 32 \
    --max-target-seqs 500

echo " DIAMOND search completed. Results saved to $OUT_TSV"
```


```python
# Step 4 Prevalence Calculation & Dispersal Categorization (Python)

#We compute transposase family prevalence as the percentage of unique genera within a lineage harboring specific families. To isolate global dispersal signals from background noise, families detected in > 10 unique genera are analyzed individually, while those detected in fewer than 10 unique genera are aggregated into the "Low Abundance" category.
```


```python
import pandas as pd

# File paths
diamond_hits_path = "./transposase_db/diamond_refseq_tnp_hits.tsv"
taxonomy_metadata_path = "Tabla_InfoseqCompleta.tsv"
out_prevalence_path = "./transposase_db/transposase_family_prevalence_summary.tsv"

print("Loading DIAMOND RefSeq search results...")
hits_cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
              "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp"]

hits_df = pd.read_csv(diamond_hits_path, sep="\t", names=hits_cols)

# Extract Genome/Subject accession and merge with genus taxonomy metadata
hits_df["Subject_GCF"] = hits_df["sseqid"].apply(lambda x: "_".join(x.split("_")[:2]))

metadata = pd.read_csv(taxonomy_metadata_path, sep="\t")
mapped_df = hits_df.merge(metadata[["GCF_ID", "genus", "phylum"]], left_on="Subject_GCF", right_on="GCF_ID", how="inner")

print("Calculating prevalence across unique genera per transposase family...")

# Count unique genera per transposase family (qseqid)
family_genus_counts = mapped_df.groupby("qseqid")["genus"].nunique().reset_index(name="Unique_Genera_Count")

# Total unique genera in dataset
total_unique_genera = metadata["genus"].nunique()
family_genus_counts["Prevalence_Pct"] = (family_genus_counts["Unique_Genera_Count"] / total_unique_genera) * 100

# Categorization protocol: > 10 unique genera vs Low Abundance (<= 10 unique genera)
def categorize_dispersal(count):
    if count > 10:
        return "Widespread Family (>10 Genera)"
    else:
        return "Low Abundance Category (<=10 Genera)"

family_genus_counts["Dispersal_Category"] = family_genus_counts["Unique_Genera_Count"].apply(categorize_dispersal)

# Aggregate Low Abundance stats
widespread_count = (family_genus_counts["Dispersal_Category"] == "Widespread Family (>10 Genera)").sum()
low_abundance_count = (family_genus_counts["Dispersal_Category"] == "Low Abundance Category (<=10 Genera)").sum()

print("\n==================================================")
print("TRANSPOSASE DISPERSAL CATEGORIZATION SUMMARY")
print("==================================================")
print(f"Total Transposase Families Analyzed: {len(family_genus_counts)}")
print(f"  -> Widespread Families (> 10 Genera): {widespread_count}")
print(f"  -> Low Abundance Category (<= 10 Genera): {low_abundance_count}")

# Save detailed prevalence table
family_genus_counts.to_csv(out_prevalence_path, sep="\t", index=False)
print(f" Summary saved to: {out_prevalence_path}")
```


```python
# Section Over
```
