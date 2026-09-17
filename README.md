# Reproduction of rNMP Mapping in the Human Mitochondrial Genome

## Overview

This repository presents an independent computational reproduction of the ribonucleotide monophosphate (rNMP) mapping analysis reported by Xu, Yang *et al.* in *Nucleic Acids Research* (2024). The study examined the distribution, strand bias, enriched zones, sequence context, and nucleotide composition of rNMPs embedded in human mitochondrial DNA (hmtDNA).

The analysis was rebuilt from the deposited ribose-seq reads and evaluated against the values reported in the publication and its supplementary tables. Across all 32 libraries, the reconstructed workflow recovers **23,073,140 rNMPs compared with 23,089,962 reported in the paper**, corresponding to a cohort-level deviation of **-0.07%** after restriction-site subtraction.

The repository reproduces six major analytical panels from the publication and includes additional validation analyses intended to distinguish numerical agreement, methodological sensitivity, and unresolved discrepancies.

## Reference study

> Xu P, Yang T, Kundnani DL, *et al.* **Light-strand bias and enriched zones of embedded ribonucleotides are associated with DNA replication and transcription in the human-mitochondrial genome.** *Nucleic Acids Research*. 2024;52(3):1207-1225. DOI: 10.1093/nar/gkad1204

## Main reproduction results

The principal results are summarized below.

| Analysis | Reproduction outcome | Validation source |
|---|---|---|
| Figure 1A, strand bias by cell type | Maximum deviation 4.49 percentage points; median deviation 0.63 percentage points; Pearson r = 0.996 | Digitized published panel |
| Figure 1B, per-nucleotide rNMP frequency | FS185: 143,269 vs 143,927 published; FS259: 4,632,030 vs 4,665,409 published | Figure labels and Table S7 |
| Figure 1C, rNMP-enriched zones | 13 of 13 common REZs recovered; 10 of 13 match the exact published library count | Table S3 |
| Figure 2, hotspots and sequence context | 84.3% of 1,536 published hotspot counts reproduced as identical integers; r = 0.9948 | Table S4 |
| Figure 5, composition and strand-bias contribution | All dominant rNMP-type trends described in the paper are reproduced | Table S2 and Methods |
| Figure 6A, dinucleotide context | 87 of 112 enrichment calls agree overall; 45 of 48 agree where N > 4 | Table S5 |

## Scientific objective

Ribonucleotides are among the most common non-canonical nucleotides incorporated into DNA. In human mitochondrial DNA, their distribution can provide information about replication, transcription, and nucleotide incorporation processes. The original study reported three major features of the mitochondrial rNMP landscape:

1. strand-specific differences in rNMP abundance,
2. recurrent rNMP-enriched zones across the mitochondrial genome,
3. nucleotide and sequence-context preferences at high-frequency incorporation sites.

This repository tests whether these findings can be recovered directly from the deposited ribose-seq data using a transparent and independently reconstructed analysis workflow.

## Data sources

| Source | Accession or file | Description |
|---|---|---|
| Ribose-seq data | PRJNA941970 | 32 libraries spanning 10 cell-type or experimental groups |
| Mitochondrial reference | GRCh38 chrM | 16,569 nt mitochondrial genome, with an additional rotated reference for circular recovery |
| Supplementary tables | Study supplementary files | Tables S1-S5 and S7 parsed to tabular files in `data/external/` |
| Restriction-site filter | `xph9876/ArtificialRiboseDetection` | Authors' code used as the reference implementation for restriction-site background subtraction |

The deposited sequencing reads appear to be already adapter-demultiplexed and UMI-extracted. For that reason, barcode and pattern extraction are not re-applied in the reconstructed Ribose-Map workflow, and PCR duplicate collapse cannot be performed from the deposited FASTQ files.

## Computational workflow

The reconstructed pipeline follows the sequence below:

```text
ENA FASTQ retrieval and md5 verification
        -> Trim Galore quality and length filtering
        -> Bowtie 2 alignment to GRCh38 chrM
        -> Ribose-Map alignment and coordinate modules
        -> reference-mismatch filtering
        -> restriction-site subtraction
        -> circular genome recovery
        -> reproduction of published summary analyses
```

The principal analysis code is contained in `src/rnmp/`, while the complete reproduction is driven from:

```text
notebooks/rnmp_hmtdna_reproduction.ipynb
```

Formulas used for normalized frequency, enrichment factor, strand bias, and strand-bias contribution are implemented in `src/rnmp/metrics.py` and tested against the published definitions.

---

# Results

## 1. Pilot validation of the reconstructed pipeline

Before processing the full cohort, the workflow was evaluated on pilot library FS310. This library was used to confirm read processing, mapping, rNMP recovery, nucleotide composition, and strand assignment.

![Pilot pipeline validation](figures/figA1_pilot_pipeline_funnel.png)

The pipeline recovers **4,798 rNMPs compared with 4,765 reported in the publication**, a difference of approximately **+0.7%**.

![Pilot composition comparison](figures/figA2_pilot_composition.png)

The reproduced rNMP composition differs from the published values by no more than **1.1 percentage points** across rAMP, rCMP, rGMP, and rUMP.

![Pilot strand split](figures/figA3_pilot_strand_split.png)

The pilot library contains 45.2% light-strand rNMPs and 54.8% heavy-strand rNMPs. This panel is included as a pipeline diagnostic only. It is based on a single library and is not intended to substitute for the replicate-based analysis used in the paper's Figure 1A.

## 2. Restriction-site subtraction

A restriction-site background filter was reconstructed from the authors' published GPL-v3 code. The implementation was validated directly against the authors' `re_build()` logic on GRCh38 chrM and produces identical masked position sets for all 18 enzymes used in the study.

![Effect of restriction-site subtraction](figures/fig_re_filter_effect.png)

Across the full cohort, restriction-site subtraction reduces the aggregate deviation from approximately **+0.71% before subtraction to -0.07% after subtraction**. Two hESC-H9 libraries that were generated without restriction enzymes remain substantial positive outliers, indicating that their excess cannot be explained by restriction-site background.

## 3. Strand bias across cell types, reproduced Figure 1A

![Reproduced Figure 1A](figures/fig1A_strand_bias_by_celltype.png)

The reconstructed strand-bias profile closely matches the published panel, with a **maximum deviation of 4.49 percentage points**, a **median deviation of 0.63 percentage points**, and **Pearson r = 0.996** across the ten cell-type or experimental groups.

Most groups show a light-strand excess. In contrast, both liver biopsy groups show the opposite pattern, with a larger fraction of rNMPs on the heavy strand. This behavior is visible in both the reproduced analysis and the original published figure.

## 4. Per-nucleotide rNMP frequency, reproduced Figure 1B

![Reproduced Figure 1B](figures/fig1B_rnmp_frequency_circular.png)

The circular frequency plots reproduce the positional structure reported for representative libraries FS185 and FS259. Spike heights are scaled using the rule described in the publication, where a full-height spike corresponds to half of the maximum count observed on that strand.

The reproduced library totals are:

- **FS185:** 143,269 rNMPs, compared with 143,927 published
- **FS259:** 4,632,030 rNMPs, compared with 4,665,409 published

Because the paper does not provide a complete table of per-position counts for these panels, this comparison is primarily a visual and total-count reproduction rather than an exact position-by-position numerical validation.

## 5. rNMP-enriched zones, reproduced Figure 1C

![Reproduced Figure 1C](figures/fig1C_rez.png)

All **13 of 13 common rNMP-enriched zones (REZs)** listed in the publication are recovered. Ten of the thirteen zones match the published number of positive libraries exactly.

One notable case is REZ5H, which is recovered at the published count only after circular genome recovery is applied. This supports the need to account for mitochondrial genome circularity when evaluating enrichment near the reference boundary.

## 6. Hotspots and local sequence context, reproduced Figure 2

![Reproduced Figure 2](figures/fig2_hotspots.png)

The reproduced hotspot analysis agrees closely with Supplementary Table S4. Of the 1,536 published library-by-position values considered, **84.3% of the reconstructed underlying counts are identical integers**, with **r = 0.9948** across the full comparison.

The sequence-context analysis also reproduces the main nucleotide preferences reported in the paper. In particular, the dominant rNMP identities for the examined cell-type subsets agree with the published claims, including the HEK293T and RNASEH2A knockout comparisons.

## 7. rNMP composition and strand-bias contribution, reproduced Figure 5

![Reproduced Figure 5](figures/fig5_composition.png)

The reconstructed composition analysis reproduces the major cell-type-specific nucleotide preferences:

- rAMP is dominant in the liver biopsy groups,
- rCMP is prominent in hESC-H9 and whole-blood libraries,
- rGMP is elevated in the RNASEH2A knockout groups,
- rUMP remains a minor component across all groups.

The strand-bias contribution heatmap shows that the overall bias in a given cell type is frequently driven by one predominant rNMP class rather than by a uniform shift across all four nucleotide types.

## 8. Dinucleotide context, reproduced Figure 6A

![Reproduced Figure 6A](figures/fig6A_dinucleotide.png)

The reconstructed dinucleotide analysis recovers the major sequence-context patterns described in the publication. C-C is strongly enriched for rCMP in many libraries, while T-A, A-G, and C-U are recurrent patterns for rAMP, rGMP, and rUMP, respectively.

Across the 112 cell-type-by-pattern significance comparisons represented in Table S5, **87 calls agree** with the publication. Agreement increases to **45 of 48 comparisons** for groups with more than four libraries.

---

# Additional findings from the reproduction

The analyses below are not simple restatements of the original paper. They emerged while reconciling the deposited data, supplementary tables, and reconstructed pipeline.

## 1. Statistical significance at small sample sizes

Some reported significance levels appear difficult to reconcile with the stated non-parametric tests and the available group sizes. The clearest examples occur in groups with N = 4 or N = 5, where the minimum attainable exact p-values are constrained by the small number of observations.

This repository therefore treats the associated significance annotations cautiously and separates numerical reproduction of the measurements from reproduction of the reported p-values.

Relevant outputs:

```text
results/tables/fig6A_nr_significance_vs_paper.csv
docs/dataset_notes.md
```

## 2. Restriction-site filtering is narrower than whole-site masking

The authors' restriction-site correction does not simply remove every nucleotide within a recognition sequence. Instead, it acts on specific strand-aware positions derived from each matched site and applies an additional condition related to non-templated dA tails for one class of events.

The reconstructed implementation matches the authors' code exactly at the level of masked genomic positions and materially improves agreement with the published cohort totals.

Relevant implementation:

```text
src/rnmp/re_sites.py
tests/test_re_sites.py
```

## 3. Circular recovery affects REZ analysis and total-count agreement differently

Circular recovery is necessary to reproduce REZ5H at its published frequency across libraries. However, including all recovered junction-spanning reads causes some per-library totals to move away from the values reported in Supplementary Table S2.

This suggests that the published REZ analysis and the published total-count table may not have been generated from exactly the same processed count set. The repository therefore retains both non-circular and circularly recovered coordinate files and uses each where it is analytically appropriate.

## 4. One SRA sample annotation is inconsistent with the supplementary tables

Run `SRR23726642` (`FS303`) is deposited in SRA as HEK293T-WT, whereas the study's supplementary tables place it in the hESC-H9 group. The group sizes reported in Figure 1A are consistent with the supplementary-table assignment.

For this reason, `config/samples_xu2024.tsv` uses the study tables as the source of cell-type labels.

## 5. The reported 70% common-hotspot rCMP proportion depends on the denominator

The paper states that approximately 70% of common hotspots are rCMP. Supplementary Table S4 contains 48 listed positions, of which 32 are rCMP, corresponding to 66.7%. A closely related subset of 46 shared positions gives 69.6%, which rounds to 70%.

The reproduction therefore supports the qualitative conclusion that rCMP is the dominant common hotspot identity, while indicating that the exact percentage depends on the set of positions used as the denominator.

---

# Limitations

Several points remain unresolved and are intentionally reported as limitations rather than forced into agreement with the publication.

- Four hESC-H9 libraries remain more than 10% above their published total after all reconstructed filters. Two of these libraries were fragmented without restriction enzymes, so the excess cannot be attributed to restriction-site background.
- Two published hotspot calls are not supported by any aligned read at the relevant position in the corresponding reproduced library. The cause therefore appears to occur upstream of the reconstructed filtering steps.
- PCR duplicate removal cannot be repeated because the deposited reads are already UMI-extracted.
- The reproduced Figure 1C is rendered as a raster representation of positive 200-nt bins rather than as the original ring-style visualization. The underlying REZ calls, not the graphical style, are the validated quantity.
- The strand-bias contribution color scale is not fully reconcilable with the scale printed in the publication, although the implemented contribution formula is internally consistent.
- Figure 6B and Figures 3, 4, and 7 were not reproduced in this repository.

---

# Reproducibility

## Quick start

The repository can regenerate all committed figures and validation tables without downloading the original FASTQ files.

```bash
conda env create -f envs/rnmp-analysis.yaml
conda activate rnmp-analysis
python -m pytest -q
python scripts/execute_notebook.py notebooks/rnmp_hmtdna_reproduction.ipynb
```

Expected outcome:

```text
71 tests passed
32 notebook code cells executed
10 inline figures generated
6 rendered tables generated
```

## Rebuilding the count data from FASTQ files

A second environment contains the aligners and is only required when rebuilding coordinates from the raw sequencing reads.

```bash
conda env create -f envs/rnmp-pipeline.yaml
conda activate rnmp-pipeline

bash scripts/00_setup_reference.sh
bash scripts/run_library.sh SRR23726637 FS310

# Full 32-library cohort
bash scripts/run_all_libraries.sh
python scripts/export_coordinates.py
python scripts/filter_re_sites.py
python scripts/04_circularize.py
```

The full rebuild downloads approximately 747 MB of sequencing data and regenerates the per-library alignment and coordinate intermediates used in the downstream analysis.

## Repository structure

```text
envs/          conda environments for analysis and alignment
config/        sample metadata and pipeline parameters
data/external/ parsed supplementary tables and external reference inputs
data/reference/ mitochondrial reference files and rotated reference
src/rnmp/      analysis, metrics, filtering, and plotting code
scripts/       download, alignment, filtering, export, and execution utilities
notebooks/     main end-to-end reproduction notebook
figures/       rendered figures used in this README
results/       per-position coordinates and validated result tables
tests/         unit and reproduction tests
docs/          detailed provenance notes and unresolved discrepancies
```

## Reproducibility boundary

The repository commits the processed coordinate files and summary tables required to reproduce the figures without rerunning the expensive alignment stage. Large intermediate SAM and BED files are intentionally excluded and can be regenerated from the raw sequencing reads using the scripts above.

The committed reproduction includes:

```text
results/coordinates/*.csv.gz
results/tables/*.csv
data/reference/chrM.fa
pilot-library alignment and coordinate BED files
```

This separation allows the analytical results to be regenerated quickly while preserving a complete path back to the deposited raw data.

---

# Interpretation

Overall, the deposited data support the major biological patterns reported in the original study. The reconstructed workflow reproduces the cohort-level rNMP count, strand-bias structure, common enriched zones, hotspot composition, nucleotide preferences, and several sequence-context effects with close numerical agreement.

At the same time, the reproduction identifies a small number of methodological and reporting issues that are important for interpreting exact numerical claims. These include sensitivity to restriction-site subtraction, circular-genome handling, small-sample significance testing, one inconsistent sample annotation, and a few remaining differences in library-level counts.

The purpose of the repository is therefore not only to recreate the published figures, but also to document which results reproduce exactly, which reproduce approximately, and which depend on analytical choices that are not fully specified in the published record.

## Licence

The analysis code in this repository is released under the MIT Licence.

`src/rnmp/re_sites.py` and `data/external/res_all.list` are distributed under GPL v3 because they are derived from the authors' `ArtificialRiboseDetection` implementation. The original paper is published under CC BY.
