# Python → Seurat snippets

Purpose: provide "side-by-side" snippets so you can replicate the intent of the Python/Scanpy code in R/Seurat.
Please consider that I am just selecting snippets, that means that there's more code pre-post this snippets
These R snippets are *not* fully adapted to your objects; you should adjust assay/column names and gene patterns as needed.

## Global assumptions (minimal)
- Seurat QC fields live in `obj@meta.data`.
- The R code below assumes these meta columns exist: `nCount_RNA`, `nFeature_RNA`, `percent.mt`, `percent.ribo`.
- If your object uses different column names (e.g. SCT assay), map them accordingly.
- `percent.ribo` must be computed by you (definition depends on species/gene naming).

---

## 1) Hard-bounds QC filter (mito, ribo, UMI, genes, novelty)

### Python (Scanpy)

```python
def filter_adata_qc_hardbounds(
    adata,
    mt_max: float = 15.0,
    ribo_min: float = 5.0,
    ribo_max: float = 50.0,
    umi_min: int = 400,
    umi_max: int = 50000,
    genes_min: int = 300,
    genes_max: int = 10000,
    novelty_min: float = 0.8,
    inplace: bool = False
):
    """
    Apply the agreed hard-bounds filter:
      - pct_counts_mt > mt_max → remove
      - pct_counts_ribo > ribo_max or < ribo_min → remove
      - total_counts < umi_min or > umi_max → remove
      - n_genes_by_counts < genes_min or > genes_max → remove
      - log10GenesPerUMI < novelty_min → remove
    Erythroid % is computed but not used for filtering.
    Returns filtered AnnData unless inplace=True.
    """
    req = ["pct_counts_mt","pct_counts_ribo","total_counts","n_genes_by_counts"]
    for k in req:
        if k not in adata.obs.columns:
            raise ValueError(f"Missing QC metric '{k}' in adata.obs. Run compute_basic_qc() first.")

    obs = adata.obs
    if "log10GenesPerUMI" not in obs.columns:
        obs["log10GenesPerUMI"] = np.log10(obs["n_genes_by_counts"] + 1) / np.log10(obs["total_counts"] + 1)

    obs["high_ribo"] = obs["pct_counts_ribo"] > ribo_max
    obs["low_ribo"] = obs["pct_counts_ribo"] < ribo_min
    obs["low_novelty"] = obs["log10GenesPerUMI"] < novelty_min

    keep = (
        (obs["pct_counts_mt"] <= mt_max) &
        (obs["pct_counts_ribo"] <= ribo_max) &
        (obs["pct_counts_ribo"] >= ribo_min) &
        (obs["total_counts"] >= umi_min) &
        (obs["total_counts"] <= umi_max) &
        (obs["n_genes_by_counts"] >= genes_min) &
        (obs["n_genes_by_counts"] <= genes_max) &
        (obs["log10GenesPerUMI"] >= novelty_min)
    )

    filtered = adata[keep].copy()
    print(f"[QC] Removed {adata.n_obs - filtered.n_obs} / {adata.n_obs} cells; kept {filtered.n_obs}.")
    print(f"[QC flags] high_ribo: {int(obs['high_ribo'].sum())}, low_ribo: {int(obs['low_ribo'].sum())}, low_novelty: {int(obs['low_novelty'].sum())}")

    if inplace:
        adata._inplace_subset_obs(keep)
        return None
    return filtered
```

### R(Seurat)

```R
filter_seurat_qc_hardbounds <- function(
  obj,
  mt_max = 15.0,
  ribo_min = 5.0,
  ribo_max = 50.0,
  umi_min = 400,
  umi_max = 50000,
  genes_min = 300,
  genes_max = 10000,
  novelty_min = 0.8,
  inplace = FALSE
) {
  md <- obj@meta.data

  req <- c("percent.mt", "percent.ribo", "nCount_RNA", "nFeature_RNA")
  miss <- setdiff(req, colnames(md))
  if (length(miss) > 0) {
    stop(paste0(
      "Missing QC meta columns: ",
      paste(miss, collapse = ", "),
      ". Compute/map them first."
    ))
  }

  if (!("log10GenesPerUMI" %in% colnames(md))) {
    md$log10GenesPerUMI <- log10(md$nFeature_RNA + 1) / log10(md$nCount_RNA + 1)
  }

  md$high_ribo   <- md$percent.ribo > ribo_max
  md$low_ribo    <- md$percent.ribo < ribo_min
  md$low_novelty <- md$log10GenesPerUMI < novelty_min

  keep <- (
    (md$percent.mt <= mt_max) &
    (md$percent.ribo <= ribo_max) &
    (md$percent.ribo >= ribo_min) &
    (md$nCount_RNA >= umi_min) &
    (md$nCount_RNA <= umi_max) &
    (md$nFeature_RNA >= genes_min) &
    (md$nFeature_RNA <= genes_max) &
    (md$log10GenesPerUMI >= novelty_min)
  )

  removed <- sum(!keep)
  kept <- sum(keep)
  message(sprintf("[QC] Removed %d / %d cells; kept %d.", removed, nrow(md), kept))
  message(sprintf(
    "[QC flags] high_ribo: %d, low_ribo: %d, low_novelty: %d",
    sum(md$high_ribo), sum(md$low_ribo), sum(md$low_novelty)
  ))

  # persist computed columns/flags
  obj@meta.data <- md

  if (inplace) {
    obj <- subset(obj, cells = rownames(md)[keep])
    return(obj)
  } else {
    obj_f <- subset(obj, cells = rownames(md)[keep])
    return(obj_f)
  }
}
```
### What you must adapt (by design)
If your Seurat object uses different QC column names, map them (e.g. percent.mt vs pct_mt).
Check for the correct assay (e.g. SCT), depending on your downstream!! 
Ensure percent.ribo is computed with the gene pattern that matches... in case you'll do it.

## 2) HVG selection (from raw counts) + HVG exclusion (confounder sets: PCs/Ig/Ribo/MT/Xist/TCR)
### Python (Scanpy): HVG selection

```python
def select_highly_variable_genes(adata, batch_key="Sample", n_top=3000):
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=n_top,
        layer="counts",
        batch_key=batch_key,
        span=0.6,
        subset=True,
    )
    return adata
```

### Python (Scanpy): exclude gene sets from the HVG space, then subset & save (for scVI)
```python
from pathlib import Path
# ---- Categories / gene sets to exclude from the HVG space (includes PCs) ----
PCS = {"Xbp1", "Jchain", "Scd1", "Irf4", "Prdm1"}

# Genes we do NOT want to influence the HVG space -> embedding/clustering/scVI
mask_exclude = (
    ad.var_names.isin(PCS) |                             # PCs (custom set)
    ad.var_names.str.startswith(("Igk", "Igl", "Igh")) |  # Ig loci
    ad.var_names.str.startswith(("Rps", "Rpl")) |         # Ribosomal
    ad.var_names.str.startswith(("mt-")) |         # Mito
    (ad.var_names == "Xist") |                            # Xist
    ad.var_names.str.startswith(("Trav", "Trbv", "Trcv")) # TCR V segments
)

# ---- Safety: ensure HVGs already computed ----
if "highly_variable" not in ad.var.columns:
    raise ValueError("Missing ad.var['highly_variable'] (compute HVGs before this block)")

# ---- Report: how many excluded genes were HVG before excluding ----
n_exclude_total = int(mask_exclude.sum())
n_exclude_in_hvg = int((mask_exclude & ad.var["highly_variable"]).sum())
print(f"Genes to exclude (total): {n_exclude_total}")
print(f"Genes to exclude that were marked as HVG: {n_exclude_in_hvg}")

# ---- Mark as NOT-HVG (does not recalculate HVGs; only flips the flag) ----
ad = ad.copy()  # avoid mutating the object in memory (recommended)
ad.var.loc[mask_exclude, "highly_variable"] = False

# ---- Subset: keep ONLY HVGs after exclusion ----
ad_hvg = ad[:, ad.var["highly_variable"]].copy()
print(f"After exclusion: final HVG genes = {ad_hvg.n_vars}")

# ---- Save for scVI (reduced object to 'clean' HVGs) ----
out_path = Path(
    "combined_singlets_HVG_excluded.h5ad"
)
ad_hvg.write(out_path)
print(f"[SAVE] {out_path}")
```

### R (Seurat): conceptual equivalence (minimal)

- HVGs live in VariableFeatures(obj) (a character vector of feature names).
- Exclusion is a set difference: hvgs_clean <- setdiff(VariableFeatures(obj), excluded_genes).


```R
# ---- Define the exclusion set (adapt patterns/species as needed) ----
PCS <- c("Xbp1","Jchain","Scd1","Irf4","Prdm1")

# Current HVG set (Seurat stores this as a character vector)
hvgs <- VariableFeatures(obj)

# Genes to exclude from the HVG space
excluded_genes <- unique(c(
  PCS,
  grep("^Igk|^Igl|^Igh", rownames(obj), value = TRUE),   
  grep("^Rps|^Rpl", rownames(obj), value = TRUE),        
  grep("^mt-", rownames(obj), value = TRUE),             
  "Xist",
  grep("^Trav|^Trbv|^Trcv", rownames(obj), value = TRUE) 
))

# ---- Report (optional but handy) ----
cat("Excluded genes (total):", length(excluded_genes), "\n")
cat("Excluded genes that were HVG:", sum(hvgs %in% excluded_genes), "\n")

# ---- Clean HVGs (this is the key action) ----
hvgs_clean <- setdiff(hvgs, excluded_genes)

# Store back into the object (now downstream steps use this HVG list)
VariableFeatures(obj) <- hvgs_clean

cat("Final HVGs after exclusion:", length(VariableFeatures(obj)), "\n")

```