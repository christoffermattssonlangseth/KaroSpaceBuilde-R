# KaroSpaceR – Implementation Plan

Goal: Achieve feature parity between the R export pipeline and the Python pipeline while ensuring the viewer supports gene expression, metadata defaults, spatial graphs, and scalable payload sizes.

---

## Phase 1: Gene Expression Parity

### Tasks
- [x] Add a proper assay/layer resolver for each supported input type.
- [x] For Seurat, default to normalized expression in `SCT` when available.
- [x] Implement practical assay fallback with slot/layer fallback for expression extraction.
- [x] Populate `available_genes` in the payload.
- [x] Populate `genes_meta` in the payload.
- [x] Export per-section gene matrices.
- [x] Add CLI/example-script option `--assay`.
- [x] Add CLI/example-script option `--genes`.

### Acceptance Criteria
- `Heart_A90_karospacer.html` can switch between cluster coloring and true gene-expression coloring.

---

## Phase 2: Better Metadata and Viewer Defaults

### Tasks
- [x] Improve auto-detection of `groupby`.
- [x] Improve auto-detection of `initial_color`.
- [x] Improve auto-detection of `additional_colors`.
- [x] Improve detection of section-level metadata.
- [x] Add optional `outline_by` support.
- [x] Add cleaner metadata filtering for single-section datasets.
- [x] Add cleaner metadata filtering for multi-section datasets.
- [x] Normalize Seurat inputs to the standard payload structure.
- [x] Normalize SingleCellExperiment inputs to the standard payload structure.
- [x] Normalize SpatialExperiment inputs to the standard payload structure.
- [x] Normalize plain-list inputs to the standard payload structure.

### Acceptance Criteria
- The same object exports with sensible defaults and minimal manual flags.

### Remaining Polish
- [x] Tighten default color selection for partial external-metadata merges.
- [x] Add clearer export-time reporting of annotation coverage for partially annotated datasets.

---

## Phase 3: Spatial Graph and Neighbor Features

### Tasks
- [x] Add optional neighbor graph export.
- [x] For Seurat, detect existing graph objects.
- [x] Prefer existing graphs when available.
- [x] Add option to derive graphs from spatial coordinates when missing.
- [x] Populate `has_neighbors` in the payload.
- [x] Export per-section edge lists.
- [x] Export `neighbor_stats`.

### Acceptance Criteria
- Neighbor overlays work in the viewer.
- Neighbor summary panels function in the R export.

---

## Phase 4: Marker Genes and Interaction Features

### Tasks
- [ ] Implement marker-gene computation for categorical groupings.
- [ ] Add marker export to viewer payload.
- [ ] Add optional interaction/contact-style marker computation.
- [ ] Ensure interaction markers depend on the graph layer.
- [ ] Keep marker computations optional due to computational cost.

### Acceptance Criteria
- Viewer marker panels show real results instead of empty states.

---

## Phase 5: Performance and Payload Size

### Tasks
- [ ] Implement sparse gene encoding for expression matrices.
- [ ] Add optional array packing for large payloads.
- [ ] Add base64 encoding for large section payloads when needed.
- [ ] Implement a `lightweight export` mode.
- [ ] Reduce HTML payload size where possible.

### Acceptance Criteria
- Large datasets export successfully.
- HTML files remain usable and do not become excessively large.

---

## Status

- Phase 1: ✅ Complete  
- Phase 2: ✅ Complete  
- Phase 3: ✅ Complete  
- Phase 4: ⬜ Not started  
- Phase 5: ⬜ Not started
