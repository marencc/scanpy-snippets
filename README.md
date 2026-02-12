# Scanpy → Seurat snippets

This repo contains a small, practical collection of "side-by-side" snippets showing how a few Scanpy/Python steps map to Seurat/R.

## What’s inside
- `python_to_seurat_snippets.md`: the main document (Python code + suggested Seurat equivalents).

## Scope
- The Seurat side is intentionally minimal: it captures the *intent* and the common Seurat idioms, but it is not fully adapted to any specific project object structure.
- You should adjust assay/column names and gene-pattern definitions as needed.

## Notes
- HVGs in Seurat are stored as a vector via `VariableFeatures(obj)`.
- If you need to physically subset to HVGs, you can do it, but it’s usually enough to keep a cleaned HVG vector.

## Contact
Open an issue or message the maintainer if you want another Python snippet or further details!
