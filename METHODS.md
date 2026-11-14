# Methods and interpretation

## Detection and labels

The final development fit uses 37 crops and 440 accepted Type 1 reference marks from Mouse 1 sections 1, 3, 7, 12, and 13. Whole-section cross-validation excludes each of those sections in turn. Test sections 5 and 9 remain excluded from the final fit. Crops are training examples; the biological sample size is one mouse.

Three Cellpose 3 detectors were fine-tuned with seeds 0, 1, and 2. The morphology filter averages three ExtraTrees classifiers; the appearance filter requires agreement from three random-feature SVM pipelines, with feature seeds 31, 47, and 89. The random feature encoders are frozen. The classifier cutoffs are in [methods/classifier_cutoffs.json](methods/classifier_cutoffs.json).

Accepted somata include recognizable dim or partially visible cell bodies. Type 2 locations have zero supervised training weight and are excluded near uncertain marks; they are neither positive nor negative labels. Type 3 marks provide explicit examples of processes and other non-cells. Unmarked locations in the additional annotation packets are not blanket negative labels. Training targets include historical supervision and approximate object masks; point labels alone do not specify exact cell boundaries.

Inference uses 512-pixel contexts around 256-pixel cores, a 1-pixel Gaussian smoothing sigma, a 12 µm Cellpose diameter, cell-probability threshold 0, flow threshold 0.4, and normalization tiles of 200 pixels. Physical pixel size is 1.449 µm. The completed application covered 8,415 tiles across 16 sections.

After filtering, two of three detectors must agree within 8 µm. Consensus uses complete-link groups, with at most one member per detector, and retains an actual predicted center as the medoid. Cross-tile deduplication uses a strict distance below 8 µm and prioritizes context margin; close detections from the same tile are retained. See [spatial.py](spatial.py) for these routines. This release does not include trained weights, classifier feature extraction, or the full detector training/inference code.

## Comparisons and dense windows

The coordinate files preserve physical-image detections and flags for membership in the existing atlas and historical tissue envelope. Comparisons retain points satisfying both flags. Regional densities use the actual retained area of each label, with hemispheres combined. `analyze.py` checks densities against the archived areas; reconstructing those areas requires the private label images and tissue-envelope procedure.

Method agreement is maximum-cardinality one-to-one matching within 9 µm, not a greedy nearest-neighbor count. The 9,064 common matches measure agreement only. Neither method is treated as truth in the full-slice comparison.

Dense review windows use a 32-pixel grid of candidate centers inside the shared comparison area. Centers are ranked by detections within 250 µm, with coordinate tie-breaking. Centers are accepted greedily at least 500 µm apart, provided there are at least 20 detections. Selection is independent for each method; both methods are counted in each selected window. There are 150 model-selected and 210 alternative-selected windows. These are exploratory two-dimensional windows, not anatomical borders, distinct three-dimensional nodes, or a reason to discard isolated detections.

The gallery's 320-pixel local windows come from the archived disagreement review. Selection considered both positive and negative count differences, at least 75% common-area coverage, and at least 20 detections from either method. These examples intentionally emphasize differences. Full-slice and local panels use fixed contrast within each comparison; local contrast differs from the full-slice display.

## Atlas provenance

The source maps are the existing GUI-refined aligned label TIFFs named `ROI1_green1_labels.tif` through `ROI1_green16_labels.tif`. They contain Allen structure IDs. No new registration is performed here. Per-file hashes in [data/slices.csv](data/slices.csv) identify the exact source maps.

The saved official Allen structure graph supplies names, acronyms, and ancestry. Every dense-window member retains its exact assigned ID; display-only pooling merges actual cortical layers into their immediate area. CA1, CA2, and CA3 remain separate hippocampal fields. Colors are custom and encode region identity, not density or confidence. Mixed windows retain their full composition; 106 of the 150 model windows span multiple displayed areas.

Correct ID lookup and ontology names do not independently validate anatomical registration. The maps do not establish precise bregma coordinates or anatomical hemisphere. Registration error can affect regional interpretation, especially near boundaries.

## Validation limits

Development cross-validation on all 440 references yielded 371 matched marks, 69 misses, and 35 unmatched detections for the candidate recipe, versus 353, 87, and 23 for the previous baseline. F1 improved from 0.865 to 0.877, while precision decreased. The original conservative gate failed, including increased unmatched detections and acceptance of explicit negative locations (2 to 3). Selection for cluster-oriented analysis did not change that result.

The neural-only change with the old filters achieved F1 0.878, comparable to the combined candidate. The classifier-only refit achieved F1 0.861. These ablations do not establish an added benefit from refitting the classifiers. The full CV results, including detector-pair controls and subgroup results, are retained in [results/development_cv.json](results/development_cv.json). Development was adaptive, so those results are not an independent final test.

The six retained test crops contain 127 marks from sections 5 and 9, with no Type 2 marks. The current full model matches 115, misses 12, and produces 13 unmatched detections at 9 µm. Tests at 6 and 12 µm and an interior-only analysis are included in [results/validation_summary.csv](results/validation_summary.csv). These test annotations had been evaluated earlier; the current run is a repeat check on the same mouse. Replacement of an unjudgeable original crop limits generalization to ambiguous tissue. Hand-count agreement is conditional on the annotation policy and retained crops.

There is no independent-animal test, randomized whole-tissue accuracy audit, or validated biological interpretation of every detection. Dense-window visualization is useful for review but does not prove model accuracy, functional connectivity, or a treatment effect.

## Release provenance

The working model is the full-development V91 fit of the V86 recipe. Whole-slice results come from `mouse1_full_v91_v96`, repeat-test results from `gfp_full_test_v92`, and corrected atlas reports from `gfp_hotspot_allen_v101`. Source names are retained for traceability, not as separate models to choose among.

[methods/source_manifest.json](methods/source_manifest.json) records byte-identical copied artifacts and their hashes. Alternative coordinate CSVs were exported from the original saved alternative points with atlas and tissue flags assigned by the existing comparison procedure. Gallery PNGs were rendered for this release. The archived environment file was recorded during the original inference workflow; the lightweight dependencies in `requirements.txt` are for this repository's postprocessing only.
