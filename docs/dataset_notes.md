# Dataset notes and deviations from the original study

## 1. Xu, Yang et al. 2024 -- PRJNA941970 (reproduction target)

49 runs, ~747 MB total. Resolved live from SRA; see `config/samples_xu2024.tsv`.

| role | n | LibraryStrategy | use |
|---|---|---|---|
| `ribose_seq` | 32 | OTHER | the 32 rNMP libraries analysed in the paper |
| `dna_seq_strandbias_control` | 13 | WGS | Methods, "Compare light and heavy strand sequencing bias" |
| `fragmentation_background` | 4 | WGS | Methods, "Background coverage calculation" (CD4+T; dsf, RE1, RE2, RE3) |

The 32 ribose-seq runs map onto the paper's 10 cell categories. Both single-end
(13 runs: CD4+T, hESC-H9, HCT116, WB-GTP) and paired-end (19 runs: liver
biopsies, HEK293T) libraries are present, which matters twice over: Ribose-Map's
`coordinate.sh` keeps only read 1 of a pair (`samtools view -f67`), and its
`alignment.sh` takes different branches depending on whether `read2` is set in
the config. `02_trim.sh` and `03_ribosemap.sh` therefore both switch on the
`layout` column rather than globbing FASTQ names.

### Known gaps in the public metadata

Two per-library parameters the pipeline needs are **not** in SRA metadata and
not in the main text:

1. **`barcode` and `pattern`** (the UMI layout consumed by `umi_tools extract`
   in Ribose-Map's Alignment module). The published libraries use Adapter.L1-L8
   with Adapter.S; the per-library assignment is in Supplementary Table S1/S2.
2. **Fragmentation method per library** (dsDNA Fragmentase, RE1, RE2, RE3, or
   Fragmentase+RE) -- Supplementary Table S2. This is required for the
   restriction-enzyme background-removal filter, since the sites to mask depend
   on which enzyme set cut that library.

Both are `TODO` in `config/samples_xu2024.tsv`. Fill them from Supplementary
Table S1/S2 or from the authors' deposited configs
(Zenodo 10.5281/zenodo.10211459) before running `03_ribosemap.sh`.

### Steps that are custom, not Ribose-Map

- **Circular recovery.** Ribose-Map treats chrM as linear, so rNMPs spanning
  position 1/16569 are lost. The paper aligns a second time to recover
  control-region rNMPs and merges the two coordinate sets. Implemented in
  `scripts/04_circularize.py` against a rotated chrM reference.
- **RE-site background removal** (Zenodo 10.5281/zenodo.8121711).
- **Reference-mismatch removal** of rNMPs whose called base disagrees with chrM.

## 2. Berglund et al. 2017 -- PRJNA354396 / GSE90054 (extension)

**This is Alk-HydEn-seq, not ribose-seq.** 44 runs, ~20 GB; see
`config/samples_berglund2017.tsv`. Design:

- 3 HeLa replicates (A, B, C) + 8 fibroblast lines: `FB1` (control) and 7
  patient lines carrying mtDNA-depletion-syndrome mutations in *MPV17*, *TK2*
  and *DGUOK*.
- Every sample appears as a **KOH** (alkali-cleaved, rNMP signal) and a matched
  **KCl** (mock, background) library -- 22 KOH / 22 KCl.
- Each of those in two fragmentation modes: `random` and `HincII`.

`background_run` in the sample sheet pairs each KOH library with its matched KCl
library on `(cell_line, fragmentation)`; all 22 pair cleanly.

### Consequences for "run the same pipeline"

The trimming/alignment half transfers directly. The rNMP-calling half does not:

1. **Different coordinate convention.** Ribose-Map supports this data via
   `technique='Alk-HydEn-seq'`, which takes the rNMP as the base 5'-adjacent to
   the read start on the same strand. Ribose-seq uses the opposite-strand
   convention. Getting `technique` wrong silently produces coordinates offset by
   one base and flipped in strand -- which would fake a strand-bias result.
2. **rNMP calling requires background subtraction.** In ribose-seq every
   sequenced fragment end *is* an rNMP. In Alk-HydEn-seq a 5' end is only
   evidence of an rNMP above the matched-KCl rate, so KOH counts must be
   normalised against the paired KCl library per position. Running
   `technique='Alk-HydEn-seq'` on KOH alone yields all free 5' ends, not rNMPs.
3. **No UMIs** -> leave `pattern`/`barcode` unset so Alignment skips
   `umi_tools`; there is no PCR-duplicate collapse available.
4. **Origin-proximal 5'-end pile-up.** 5'-end methods pile up at mtDNA
   replication origins and at the HincII sites, which overlaps exactly where the
   paper reports control-region REZs. Any REZ called in the Berglund data needs
   the KCl track shown alongside it.

### The comparison is not a like-for-like validation

Xu et al. already discuss this dataset (their ref. 12) without reprocessing it,
and report that Berglund et al. found elevated light-strand rNMPs in normal
fibroblasts but **no major strand bias in HeLa**. So a strand-bias disagreement
between the two datasets is an expected outcome, and the interesting question is
whether it is technique bias or biology. Treat the following as required
controls rather than optional extras:

- KCl-only strand ratio (pure technique/background bias, no rNMP signal)
- `random` vs `HincII` fragmentation as an internal replicate of the same DNA
- the paper's own DNA-seq strand-bias control (PRJNA941970 WGS runs) as the
  ribose-seq-side equivalent

A third option, if a cleaner methodological contrast is wanted: Xu et al.
reprocessed emRiboSeq data from GSE64521, but only for *S. cerevisiae* mtDNA
(their Supplementary Figure S1B), never for human. Human emRiboSeq is therefore
also unclaimed ground.

## 3. The EF normalisation: the paper uses TWO conventions

**The paper applies two different definitions of `R_total,L` in two different
analyses, and says so in neither.** Each is established here against the paper's
own published numbers, and neither convention reproduces both.

| analysis | `R_total,L` | `bin_enrichment(normalize_by=...)` | evidence |
|---|---|---|---|
| single-nucleotide hotspot EF (Suppl. Table S4, Fig 2) | pooled, both strands | `"genome"` | all 1,536 published EF values invert to exact integer counts |
| 200-nt REZ bins (Suppl. Table S3, Fig 1C) | per-strand | `"strand"` | 12 of 13 published common REZs recovered; 0 of 13 under pooling |

### Why pooled is right for the hotspots

At a single nucleotide `Length(G) = 1`, so the EF formula reduces to

    EF = ( R_pos / R_total,L ) x 16569

Table S4 prints that EF for 48 common hotspot positions in each of the 32
libraries, and Table S2 prints each library's total rNMP count. Inverting every
one of those 1,536 values with the **pooled** total recovers an exact integer
rNMP count (max residual 2.8e-14, i.e. floating-point noise). A per-strand total
does not. Locked by `tests/test_published_values.py`.

### Why per-strand is right for the REZ bins

Applying the pooled convention to the 200-nt bins yields **zero** common REZs:
no bin on either strand reaches the ">80% of libraries" threshold (the maximum
attained is 78.1%). This is not a marginal failure, and the cause is mechanical
rather than numerical.

Under pooling, a light-strand bin is measured against a denominator dominated by
whichever strand the library favours. Eight of the 32 libraries are the liver
biopsies (DLTB, TLTB), which are strongly *heavy*-strand biased -- only 23-27% of
their rNMPs are on the light strand (see `results/tables/fig1A_strand_bias_by_celltype.csv`).
Averaged over the eight light-strand zones that Table S3 lists, those libraries
clear EF > 1 in just 1.5-1.75 of them, versus 6.2 for CD4+T and 8.0 for
WB-GTP-control. Any light-strand zone is therefore vetoed by the liver libraries
before the 80% criterion can be met.

Switching to per-strand `R_total,L` recovers **12 of the 13** published common
REZs at their exact bin indices, with per-library counts matching Table S3's own
`n_libraries` within 0-5 -- several exactly (27/27, 32/32, 26/26, 31/31, 32/32).
See `results/tables/fig1C_rez_vs_paper.csv`.

### The one miss, and the two extras

- **REZ5H** (16400-16569, heavy) is reported in 26/32 libraries; we find 21. This
  is the truncated 169-nt final bin at the end of the coordinate system, and
  circular recovery is not yet applied, so rNMPs on reads spanning position 1 are
  lost. Its light-strand partner **REZ8L occupies the same bin and matches exactly**
  (30/32 in both), so the deficit is strand-asymmetric -- the signature of losing
  origin-spanning reads, not of a threshold artefact.
- We additionally call light-strand bin 81 (87.5%) and heavy-strand bin 65
  (81.25%). Both sit just above the cutoff, and at each of those bin indices the
  paper calls the *opposite* strand, so these are near-threshold disagreements
  consistent with the missing restriction-site filter (§4).

### What "enriched in all cell categories" means

The common-REZ definition requires enrichment "in all cell categories and
genotypes, and more than 80% of the libraries". `common_rez(category_rule=...)`
implements both readings, and the paper's is forced rather than chosen: under
`"all"` (every library of every category enriched) a common REZ would need 100% of
libraries, making the separate >80% criterion vacuous -- yet Table S3 lists zones
at 26/32 libraries (81.25%). `"any"` (at least one library per category) is
therefore the paper's rule, and is the default. Running `"all"` against this
cohort yields zero common REZs.

### Consequence for reporting

Because the two conventions place the EF = 1 baseline differently -- pooled puts a
uniform library at EF = 0.5, per-strand at EF = 1.0 -- the threshold "EF > 1" does
not mean the same thing in Figure 1C as it does in Figure 2. Any statement about
"enrichment" must name which analysis it refers to.

## 4. Supplementary parameters -- OBTAINED

All four blocking parameters are now in `data/external/`, extracted from
`gkad1204_supplemental_files.zip`.

**Provenance.** PMC serves that zip only behind a JavaScript proof-of-work bot
check, which was not circumvented, and the legacy PMC FTP dataset files were
removed in August 2026. It was retrieved instead from NLM's sanctioned AWS Open
Data mirror, which requires no credentials:

    https://pmc-oa-opendata.s3.amazonaws.com/PMC10853789.1/gkad1204_supplemental_files.zip

The article is CC-BY. Files are legacy BIFF `.xls`, so parsing needs `xlrd`
(already in `envs/rnmp-analysis.yaml`).

| file | from | contents |
|---|---|---|
| `table_S1_adapters.tsv` | S1 | Adapter.L1-L8, their 3-nt tags and reverse complements |
| `table_S2_libraries.tsv` | S2 | per-library barcode, fragmentation, rNMP total, composition |
| `table_S3_common_rez.tsv` | S3 | the 13 common REZs with published EF -- validation target |
| `table_S4_common_hotspots_ef.tsv` | S4 | 48 hotspot positions x 32 libraries of published EF |
| `table_S7_max_counts.tsv` | S7 | max per-strand rNMP count per library (Fig 1B scaling) |

`config/samples_xu2024.tsv` is now fully populated for all 32 ribose-seq
libraries. The reported rNMP totals sum to 23,089,962.

### The UMI pattern

Adapter.L1-L8 all have the 5' structure `P-NN <3-nt tag> NNNNNN AGATCGGAAGAGC...`,
i.e. 2 + 3 + 6 = 11 nt before the Illumina sequence. Read 1 reads that region in
the opposite orientation, giving `pattern='NNNNNNXXXNN'` for `umi_tools extract`.
Two independent checks support this:

- Ribose-Map's own shipped example config (`lib/SRR11364933.config`) uses exactly
  `pattern='NNNNNNXXXNN'` with `barcode='TCA'`, and `TCA` is the reverse
  complement of Adapter.L2's `TGA`. That example is a ribose-seq library from
  this same lab and adapter set.
- `TCA` also appears as the Table S2 barcode for library FS264.

### Table S2's barcode column is not internally consistent

Every S2 barcode is an Adapter.L tag **or its reverse complement**, but the table
mixes the two forms: the liver libraries (FS257-FS264) report the reverse
complement while every other library reports the tag as synthesised. Six of the
32 are ambiguous either way, because `AGC` and `GCT` are each other's reverse
complements (Adapter.L5 / Adapter.L6).

This is recorded per library in the `barcode_form` column
(`tag` / `revcomp` / `ambiguous`). **`03_ribosemap.sh` must not trust the column.**
It determines the barcode empirically -- tallying the 3-mer at read-1 positions
7-9, which `pattern='NNNNNNXXXNN'` places at the barcode -- and fails with the
observed distribution if the top 3-mer disagrees with S2.

The failure mode if this is wrong is at least loud rather than silent:
Ribose-Map's `alignment.sh` demultiplexes with `grep -B 1 -A 2 ^$barcode`, so a
wrong barcode yields an empty FASTQ and zero aligned reads.

### Fragmentation sets

S2 gives the enzyme combination verbatim for every library; these are normalised
into a `fragmentation` code plus a semicolon-separated `enzymes` list (which is
what the RE-site background filter actually consumes):

| code | n | enzymes |
|---|---|---|
| `dsF` | 10 | NEBNext dsDNA Fragmentase only |
| `dsF+RE1` | 8 | Fragmentase + RE1 |
| `RE1` | 4 | HpyCH4V + Hpy166II, Eco53KI, RsaI, StuI |
| `RE3` | 4 | CviKI-1 + MlyI, MscI, MslI |
| `RE4` | 4 | AluI + Hpy166II, PmlI, StuI + Cac8I, HpyCH4V, PmeI |
| `RE2` | 2 | AleI, AluI, PvuII + DraI, HaeIII, SspI |

Two discrepancies against the main text, both harmless once known:

1. **`RE4` is not described in the main text.** The Methods list only RE1, RE2
   and RE3, but the four liver libraries FS261/FS262/FS263/FS264 used a fourth
   combination. Any RE-site masking for those libraries must come from S2, not
   from the Methods.
2. The Methods name the RE1 enzyme **Hpy116II**; S2 and the enzyme catalogues
   say **Hpy166II**. S2 is correct.

### Trinucleotide null expectation

The Methods give 0.0625 for the trinucleotide Mann-Whitney null; Supplementary
Table S5's footnote gives 0.125. The Methods are correct -- sheet `human_nnr`
carries all 64 trinucleotides, i.e. 16 neighbour contexts per rNMP, and the
per-rNMP normalisation makes the null 1/16.

### Dinucleotide orientation

Fixed from the S5 sheet layouts rather than inferred: `human_nr` is keyed
`AA, CA, GA, TA, AC, ...` (rNMP is the **second** character, neighbour is 5'),
`human_rn` is keyed `AA, AC, AG, AT, CA, ...` (rNMP **first**, neighbour 3').
Recorded in `config/params.yaml` under `patterns:`.

## 5. Scope decision for the extension (settled)

**Berglund et al. 2017 is retained as the extension dataset, with technique
controls treated as first-class results rather than supplements.** Concretely:

1. Process all 22 KOH libraries with `technique='Alk-HydEn-seq'`, then call rNMPs
   by per-position background subtraction against the matched KCl library named
   in the `background_run` column.
2. Every cross-dataset panel must be accompanied by:
   - the **KCl-only strand ratio** -- pure technique/background bias with no rNMP
     signal. If this already favours one strand, the KOH strand ratio cannot be
     read as biology.
   - **`random` vs `HincII`** for the same DNA, as an internal replicate. A REZ
     or hotspot that appears in only one fragmentation mode is a fragmentation
     artefact.
   - the **origin/HincII-proximal mask**, because 5'-end pile-up at the
     replication origins overlaps the paper's control-region REZs.
   - the Xu-side equivalent: the **13 WGS strand-bias control runs** in
     PRJNA941970 (Methods, "Compare light and heavy strand sequencing bias").
3. The HeLa strand-bias discrepancy is the analytical centrepiece, not a failed
   reproduction. Xu et al. report light-strand bias in most cell types; Berglund
   et al. report none in HeLa but light-strand enrichment in normal fibroblasts.
   The deliverable is a defensible attribution of that difference to capture
   chemistry, cell type, or both -- with the KCl track as the discriminator.

Rejected alternatives, for the record: HeLa-only (drops the fibroblast arm, which
is the half that agrees with Xu et al.); switching to human emRiboSeq from
GSE64521 (unclaimed ground, but loses the KOH/KCl internal control); running both
(doubles compute and risks thinning the reproduction half).

## 6. Library FS303 is mislabelled in SRA

SRA's `SampleName` for run SRR23726642 (library **FS303**) is `HEK293T-WT`.
Supplementary Table S2 assigns FS303 to **hESC-H9**, and Table S4's column header
reads `FS303_hESC-H9`. The paper's Figure 1A settles which is right: it reports
hESC-H9 with N = 5 and HEK293T with N = 3, and only the Table S2 assignment gives
those counts (hESC-H9 = FS197, FS198, FS199, FS201, FS303; HEK293T WT = FS203,
FS305, FS326).

`config/samples_xu2024.tsv` therefore derives `cell_category` from Table S2, not
from SRA metadata, and carries `cell_type_s2` / `genotype_s2` columns so the
provenance is visible. Left uncorrected this silently moves one library between
two cell types in every per-cell-type figure; with 5 and 3 libraries in the
affected groups it shifts both means materially.

## 7. The published significance levels in Figure 1A are not attainable

The Figure 1A caption specifies a two-tailed Mann-Whitney U test "for all cell
types with N >= 3", where N is the number of libraries (printed under each bar and
confirmed against Table S2). The exact test has a hard floor at those group sizes:

| group size | smallest attainable two-tailed exact p | best possible stars |
|---|---|---|
| 3 vs 3 | 0.1000 | ns |
| 4 vs 4 | 0.0286 | * |
| 5 vs 5 | 0.0079 | ** |

Against the asterisks printed in the panel:

| cell type | N | paper | floor at that N | attainable? |
|---|---|---|---|---|
| hESC-H9 | 5 | **** | ** | no |
| TLTB | 4 | **** | * | no |
| DLTB | 4 | ** | * | no |
| CD4+T | 5 | ** | ** | yes -- reproduced |
| RNH2A-KO T3-8 | 4 | * | * | yes (we get p = 0.057, just outside) |
| HEK293T | 3 | ns | ns | yes -- reproduced |
| RNH2A-KO T3-17 | 4 | ns | ns | yes -- reproduced |

### Is it just a different test variant? No.

The natural objection is that they used the normal approximation rather than the
exact test. Checked in R directly, on maximally separated groups:

| group size | `exact=TRUE` | normal + continuity | normal, `correct=FALSE` |
|---|---|---|---|
| 3 vs 3 | 0.1000 | 0.0809 | 0.0495 |
| 4 vs 4 | 0.0286 | 0.0304 | 0.0209 |
| 5 vs 5 | 0.0079 | 0.0122 | 0.0090 |

Even the most permissive variant bottoms out at p = 0.0090 at N = 5 and 0.0209 at
N = 4 -- still `**` and `*`. Because the smallest attainable p is bounded by
z_max = mu/sigma = n*sqrt(3/(2n+1)), reaching p < 1e-4 needs **n >= 11 per group**
and p < 1e-5 needs **n >= 14**. The largest cell type here has N = 5. R also
selects `exact=TRUE` by default at these sizes, so the approximation would not
even be reached without forcing it.

One detail supports the reading that they did use an exact test somewhere: at
N = 3 the uncorrected approximation gives p = 0.0495, which would print `*`, but
the paper shows `ns` for HEK293T. Their non-significant calls match exact-test
behaviour while their strongest calls are unattainable under any variant.

So the test cannot have been run on per-library strand percentages. The unit must
have had far more observations -- most plausibly per-bin or per-nucleotide values
pooled within each cell type -- but the Methods do not say. The asterisks in
`figures/fig1A_strand_bias_by_celltype.png` are **ours**, computed on the unit the
caption names; direction and ordering reproduce, the exact levels do not and
cannot. Resolving this needs the authors' analysis scripts
(Zenodo 10.5281/zenodo.8121711).

## 8. Restriction-enzyme background subtraction -- IMPLEMENTED

Previously bracketed rather than reproduced, because masking whole recognition
sites undershot the published counts (4,443 vs 4,765 for FS310) while masking
nothing overshot them (4,798). The bracket was resolved by reading the authors'
own code rather than tuning a mask to hit the number.

**Source.** The Data availability section names it: "The scripts used to remove
rNMPs at restriction enzyme locations are available under the GNU GPL V3.0
license on Zenodo (https://doi.org/10.5281/zenodo.8121711)", and the Methods add
that the protocol follows Balachander et al., Nat Commun 2020;11:2447 (ref. 21).
Zenodo was returning 504 on every endpoint when this was done, so the code was
obtained from the GitHub repository the deposit snapshots:
`xph9876/ArtificialRiboseDetection` (Penghao Xu, first author).

**The rule.** Far narrower than whole-site masking. For each match of a
recognition pattern, take the *midpoint* of the matched site and mask exactly two
positions, split by strand into two classes (0-based BED coordinates):

    mid = (match_start + match_end) // 2

    noda : (mid-1, '+') and (mid, '-')   removed unconditionally
    da   : (mid-1, '-') and (mid, '+')   removed only if the read shows a dA tail

`noda` positions are the 3'-terminal base of a blunt-cut fragment end on each
strand -- an undigested end that survives T5 exonuclease and is miscalled as an
rNMP. `da` positions are one base further in, where the called rNMP lands when
the dA-tailing step added a non-templated dA to that same end; the confirmation
pattern is `'T' + recognition_pattern[cut:]` matched at the read start.

**Port fidelity.** `src/rnmp/re_sites.py` is a port of their `re_build()`. Being
a derivative of GPL-v3 code, that module is GPL v3; `data/external/res_all.list`
is copied verbatim from the same repository. The port was verified by running
their function directly on GRCh38 chrM: identical position sets for all 18
enzymes used in the study. `tests/test_re_sites.py` pins this.

Two faithfulness decisions, both deliberate:

1. The `cut-after` column is used ONLY for the dA pattern; masked positions come
   from the recognition midpoint regardless. For the 17 blunt cutters here
   midpoint == cut site, but MlyI (GAGTCNNNNN, cuts after base 10) cuts outside
   its recognition sequence, so for MlyI it does not. Reproduced as written --
   the published numbers came from this code, and "correcting" it would change
   them. A test guards against a well-meaning fix.
2. The authors skip `d = 11` leading bases (3-nt barcode + 8-nt UMI) before
   matching the dA pattern. The deposited reads are already demultiplexed and
   UMI-extracted (S2), so `DA_OFFSET = 0` here.

**Effect.** Cohort total 23,252,983 -> 23,073,140 against 23,089,962 published:
**+0.71% -> -0.07%**, a residual of 16,822 in 23 million. Mean per-library
|deviation| 6.10% -> 4.90%; composition deviation vs Table S2 falls from a mean
of 5.11 pp to 2.97 pp. Removed 153,767 `noda` and 26,076 dA-confirmed rNMPs.

**Filter order does not matter.** The Methods order is circular merge -> RE
removal -> mismatch removal; this pipeline runs mismatch removal first. Verified
on FS310: the kept set is byte-identical either way (4,679 rNMPs, same index
set). Only the *attribution* differs -- running mismatch first absorbs the
dA-tailed reads, because a non-templated dA makes the called base mismatch the
reference. That is why `re_filter_audit.csv` understates the dA class for some
libraries (FS310 shows 0 dA removals, but 15 when RE runs first).

**What it does not fix.** Four libraries stay above +10%, all hESC-H9. The
filter is a no-op on two of them by construction. That residual is not part of
this finding and is not explained by it -- see §14, where it is scoped as an open
problem.

**A consistency check the filter passes.** None of the 48 published hotspot
positions falls inside any library's mask: 0 of 1,056 (position, library) pairs.
That is required rather than fortunate -- the authors ran this filter *before*
calling hotspots, so a mask that covered a published hotspot would be wrong.

## 9. The "70% of common hotspots are rCMPs" figure

The Results state that "70% of the common-rNMP hotspots are rCMPs". That number
is not reproducible from the paper's own supplementary data:

| set | rCMP fraction |
|---|---|
| Suppl. Table S4, all 48 common hotspots | **32/48 = 66.7%** |
| Suppl. Table S4, its own top 20 by median EF | 13/20 = 65% |
| our 51 common hotspots | 33/51 = 64.7% |
| the 46 hotspots common to both sets | 32/46 = **69.6%** |

The direction is solid and not in question -- rCMP is the dominant hotspot base at
roughly two-thirds by every slicing, against a chrM background where C is not the
most abundant base. But the specific "70%" appears to be the shared-subset figure
(69.6%) or a rounding of 66.7%, not a value the published table yields. Treated
here as confirmed in substance, flagged in the exact figure.

## 10. Hotspot EF: the counts reproduce, the denominators do not

Validating our single-nucleotide EF against all 1,536 published values in
Suppl. Table S4 (48 positions x 32 libraries):

| metric | before RE filter | after RE filter |
|---|---|---|
| EF relative deviation, our totals | median 1.50% | median **1.30%** |
| EF relative deviation, published totals | median 0.000% | median 0.000% |
| rNMP counts identical to published | 1295/1536 = 84.3% | 1295/1536 = 84.3% |
| counts within +/-1 | 96.1% | 96.1% |
| Pearson r | 0.99399 | 0.99477 |

The hotspot counts are unchanged by the filter, and that is structural: no
published hotspot position is masked in any library (S8), and none of our five
extra hotspots is masked either. The hotspot set is identical before and after --
51 positions, 46 shared with Table S4.

### CORRECTION: r = -0.995 is an identity, not evidence of cause

An earlier version of this document read the correlation between a library's
count deviation and its median hotspot EF deviation (r = -0.995) as evidence that
the missing RE-site filter was the *single cause* of the residual. That inference
was wrong, and implementing the filter is what exposed it: after the filter the
correlation is **r = -0.997** -- it did not fall at all.

The reason is algebraic. At a single nucleotide `EF = (R_pos / R_total,L) x 16569`.
If our count at a position equals theirs, then

    EF_ours / EF_paper = T_paper / T_ours = 1 / (1 + delta)

so the relative EF deviation is forced to `-delta / (1 + delta)`, a deterministic
monotone function of that library's count deviation. Checked against the data:
observed median EF deviation matches that expression to **0.0000 pp for all 32
libraries**, correlation 1.000000.

So r ~= -1 localises the error to the *denominator* and says nothing about why
the denominator is wrong. It cannot fall while any library total deviates,
whatever the cause. The quantities that actually moved are the magnitudes: median
|EF deviation| 1.50% -> 1.30%, mean |library-total deviation| 6.10% -> 4.90%,
cohort total +0.71% -> -0.07%.

## 11. Common-hotspot set: 51 vs the published 48

Applying the Methods definition (>= 1 rNMP at that position in *all* 32
libraries) gives 51 positions, overlapping the published 48 at 46. We call five
the paper does not -- (+, 386), (+, 1781), (+, 13134), (-, 6799), (-, 14909) --
and miss two it does: (+, 7018) and (-, 14923). Consistent with the count
differences in S10 rather than a definitional difference: the position/strand
conventions are confirmed exactly (S12).

## 12. Position, strand and reference conventions confirmed jointly

Suppl. Table S4 carries an `rnmp` column (the ribonucleotide identity) and a
`Pattern(NNNRNNN)` column (the 7-nt genomic context). Both were recomputed from
GRCh38 chrM using our coordinate and strand conventions:

- rNMP base matches at **48/48** positions
- 7-nt context matches at **48/48** positions

These are not redundant. The base identity jointly tests the coordinate, the
strand assignment and the reference; the context additionally pins the
orientation, since a reverse-complement error would break the 7-mer while leaving
the central base correct. Together they are the strongest available confirmation
that the ribose-seq coordinate convention in `src/rnmp/pipeline.py` is right.

## 13. The two missing hotspots are neither a filter nor a circularity artefact

Two of the 48 published common hotspots are still not called: (+, 7018) and
(-, 14923). Both candidate explanations were tested and both fail.

**Not circular recovery.** They sit 7,017 and 1,646 nt from the chrM 1/16569
junction. A junction-spanning read can only contribute rNMPs within one read
length (<= 151 nt) of the junction, so recovering reads across position 1 cannot
place an rNMP at either position. Wiring in circular recovery would not touch
them -- worth stating because it removes the risk of crediting that fix with a
gain it cannot produce.

**Not the RE-site filter either.** Each position fails on exactly one library:
FS310 at (+, 7018) and FS303 at (-, 14923). In those libraries our data contains
**no aligned read at all** at that position -- not a filtered one, not a
reference-mismatched one -- while the paper's own EF inverts to 1 rNMP for FS310
and 5 for FS303. The filter only removes rNMPs, so it cannot supply a read we
never had. Confirmed after implementing it: the hotspot set is unchanged.

The residual cause is therefore upstream of every filter -- reads present in
their alignment and absent from ours. Both are single-digit counts in one library
each, so it is a handful of reads, not a systematic offset: candidates are the
circular-merge step contributing reads elsewhere in the genome, a difference in
Bowtie 2 invocation, or the deposited FASTQs differing slightly from what was
aligned for publication. Not resolved.

## 14. UNRESOLVED: the FS197 / FS198 excess has no identified cause

This is the largest unexplained discrepancy in the reproduction. It is recorded
separately from §8 on purpose: implementing the restriction-site filter did not
explain it, and folding it into that finding would overstate what the filter
accounts for.

| library | cell type | fragmentation | enzymes | deviation before | after |
|---|---|---|---|---|---|
| FS197 | hESC-H9 | dsDNA Fragmentase only | none | +13.63% | **+13.63%** |
| FS198 | hESC-H9 | dsDNA Fragmentase only | none | +26.82% | **+26.82%** |
| FS199 | hESC-H9 | RE1 (5 enzymes) | 5 | +40.57% | +32.57% |
| FS201 | hESC-H9 | RE1 (5 enzymes) | 5 | +40.18% | +34.59% |

**Why FS197/FS198 are the informative pair.** Both were fragmented with dsDNA
Fragmentase alone and used no restriction enzymes at all, so there is no RE-site
background in them to remove -- the filter removes exactly zero rNMPs from each,
and their deviation is bit-identical before and after. Yet they recover 13.6% and
26.8% *more* rNMPs than published. Whatever inflates them is therefore a
different mechanism from the one §8 fixes, and it cannot be diagnosed by
improving that filter.

**What is known.**

- It is not a cohort-wide effect. 23 of 32 libraries sit within 3% of their own
  published count, and the cohort total is -0.07%. This is four libraries in one
  cell type out of ten.
- It is not the reference-mismatch filter. That filter is applied and its removal
  rate in these libraries is unremarkable.
- It is not a cell-type-assignment error. The fifth hESC-H9 library, FS303, is at
  -1.59%. If hESC-H9 were mislabelled, FS303 would deviate too.
- It is not confined to counts. hESC-H9 is also the only cell type whose Figure 1A
  light-strand percentage misses the published panel by more than 2.3 pp (4.49 pp).
- It does not propagate to the hotspot analysis. Hotspot counts at the 48
  published positions are identical integers in 84.3% of comparisons including
  these libraries; the excess is distributed across the genome, not concentrated
  at enriched positions.

**Candidate explanations, none tested.** An extra per-library filter applied by
the authors and not described in the Methods; a difference in how duplicate reads
were handled for these specific libraries (no UMI-based deduplication is possible
from the deposited reads, §1); or the deposited FASTQs differing from what was
aligned for publication. Distinguishing these needs either the authors' per-library
logs or their intermediate BAMs, neither of which is deposited.

**Consequence for the reproduction.** Every hESC-H9 aggregate in this repository
should be read as having a wider error bar than the other nine cell types. The
Figure 1A panel and `per_library_rnmp.csv` both carry the per-library deviations,
so a reader can see which points are affected.


## 15. hESC-H9 top-20 hotspots: 18/20 rCMP vs the paper's 19/20

Table S4 has a per-subset `hESC-H9` sheet (560 common hotspots) in addition to the
`All libraries` sheet, so this is directly checkable rather than inferred.

The one-position difference is **neither a base-assignment error nor a ranking tie
at the top-20 boundary**. It is a set-membership difference well inside the top 20:

- **(-, 13354), an rAMP**, ranks **3rd** in our hESC-H9 set by median EF
  (EF 30.29) and is absent from the paper's hESC-H9 hotspot set entirely.
- Our top 20 also omits their (+, 7037) and (+, 12488), both rCMP, and adds
  (+, 5929), rCMP. Net effect: -1 rCMP, +1 rAMP -- exactly the 19 -> 18 change.

Base assignment is not the issue: across all 536 hotspots both sets share, the
rNMP identity agrees 536/536 and the 7-nt context 536/536. The single apparent
mismatch, (-, 5373), is notation -- the paper writes `U` where this repository
writes `T` (`TGAUTGA` vs `TGATTGA`) -- and it sits far outside the top 20
(EF 1.75).

Worth recording even though it is not the cause here: the paper's own top-20
boundary is tie-dense. Six positions tie at exactly EF-median 19.2094, filling
their ranks 15-20, and rank 21 is 19.12. Any reordering there is arbitrary, so
the "19/20" figure is not robust to a change of tie-break rule even in their own
data. Our boundary is similarly tied (four positions at 15.1466).

Not chased further: why (-, 13354) is common in our hESC-H9 libraries and not
theirs. It is one position in one subset and does not affect the all-libraries
hotspot set (§11), where our five extras are a different five positions.

## 16. Circular recovery -- implemented, and it forces a two-count-set design

Ribose-Map treats chrM as a linear contig, so a read crossing position 1/16569
cannot align end-to-end and is lost with its rNMP. `scripts/04_circularize.py`
implements the Methods' fix (second alignment against a rotated chrM, merge).

**No FASTQ is re-fetched.** Bowtie 2 runs end-to-end, so junction-spanning reads
are present in every existing `aligned.sam` with flag 4 rather than discarded.
Cohort-wide that is 205,507 reads. They are extracted, aligned against
`chrM_rot8284`, and mapped back with
`canonical = ((rotated - 1 + 8284) mod 16569) + 1`.

**Only junction-crossing alignments are kept.** In practice this filter removes
nothing -- and that is a consistency check, not a wasted guard. A read that
aligns to the rotated reference but not the linear one *must* cross the junction,
because the two references contain identical sequence everywhere else. If the
filter had removed anything, the rotation or the coordinate map would be wrong.
Recovered rNMPs then pass the reference-mismatch and restriction-site filters.

**Effect:** +184,035 rNMPs (0.8% of the cohort), every one within 135 nt of the
junction, zero changes beyond 200 nt.

### The evidence points both ways, and that is the finding

| quantity | pre-circular | post-circular |
|---|---|---|
| cohort total vs Table S2 | **-0.07%** | +0.72% |
| libraries moved away from published count | -- | 22 of 32 |
| REZ5H library count (published **26/32**) | 20/32 (miss) | **26/32** (exact) |
| published REZs recovered (Table S3) | 12/13 | **13/13** |
| total REZ library-count deviation, 13 zones | 11 | **5** |
| Table S4 hotspot EF, median abs deviation | 1.30% | 1.24% |
| hESC-H9 top-20 rCMP (paper 19/20) | 18/20 | 17/20 |

REZ5H is the decisive case. Before recovery it was the only missed REZ, at 20 of
32 libraries against the paper's 26 -- below the 80% threshold. After recovery it
is at **exactly 26/32**: the published number, not merely a number past the
threshold, and all twelve other zones are unchanged. A fix that moves precisely
the quantity it was meant to move and nothing else is about as well-validated as
a reproduction step gets.

But the same set pushes 22 of 32 library totals *away* from Table S2 and the
cohort from -0.07% to +0.72%. Both cannot describe one count set, so the
conclusion is about the paper rather than about this pipeline: **Table S2's rNMP
counts do not include circularly-recovered reads, while the REZ analysis does.**
The paper's published outputs are not all drawn from the same stage of its own
pipeline.

This repository therefore keeps both files and uses each where the paper's tables
place it -- stated explicitly rather than silently:

    results/coordinates/per_position_counts_refiltered.csv.gz
        totals, composition, hotspots (B1, B2, B4)   -- matches Table S2 / S4
    results/coordinates/per_position_counts_circular.csv.gz
        REZ detection (B3, Figure 1C)                -- matches Table S3

Deliberately NOT done: tuning which reads to recover so that both the totals and
the REZs agree. Any such rule would be fitted to the answer, and there is no
statement in the Methods to support one.

**Consistency check requested and passed.** The two missing common hotspots
(+ 7018) and (- 14923) sit 7,017 and 1,646 nt from the junction. After recovery
both remain at zero in the library that fails them (FS310, FS303 respectively),
and the common-hotspot set is unchanged at 51 positions with 46 shared. §13
predicted this before the step was written, so the REZ gain cannot be
misattributed and this fix cannot be credited with closing the hotspot gap.

## 17. Figures 1B and 5 -- reproduced, with one scale ambiguity

**Figure 1B** (two circular panels, CD4+T FS185 and DLTB FS259) reproduces
visually. Panel identity is confirmed by the totals: 143,269 against the paper's
printed 143,927 for FS185, and 4,632,030 against 4,665,409 for FS259. Radial
scaling follows the caption -- a full-height spike is HALF the library's maximum
count on that strand, from Suppl. Table S7. No per-position values are published
for these panels, so this is a visual reproduction and is labelled as one.

**Figure 5** reproduces every qualitative claim. Dominant rNMP type by cell type:
rAMP for CD4+T, DLTB (0.756), TLTB (0.757) and HCT116; rCMP for hESC-H9 (0.522),
both whole-blood samples and HEK293T; rGMP for both RNASEH2A knockouts (0.40,
0.43 -- the cohort's highest). rUMP is minor everywhere (0.013-0.051).

Panel C reproduces the paper's striking structure: the contribution heatmap is
essentially one column per cell type, so a cell type's strand bias is carried by
a single nucleotide rather than by a general excess.

**The one ambiguity.** The paper's colour scale runs 0-30. Under the Methods as
written the largest per-library contribution here is 55.0, and cell-type means
reach 34.5 (DLTB) and 41.5 (TLTB), so those cells saturate. The arithmetic is
internally exact -- each library's four contributions sum to its strand-bias
level to floating-point precision -- so the discrepancy is in the scale, not the
formula. Either their scale saturates too, or a final normalisation step differs
from the Methods text. Not resolvable from the paper; left as written rather than
rescaled to fit.

**An asymmetry the panel exposes.** The type dominating a cell type's
*composition* is not always the type carrying its *strand bias*.
RNH2A-KO-T3-8 is rGMP-dominant by composition (0.40) but rCMP-dominant by
contribution; only RNH2A-KO-T3-17 is rGMP on both.

### A duplication caught by the existing test suite

`composition_normalized_frequency`, `dinucleotide_normalized_frequency` and
`strand_bias_contribution` were already implemented in `src/rnmp/metrics.py`
from the Methods, with tests. A second set was written for this round and
appended to the same module, where the new `strand_bias_contribution` shadowed
the original and broke its test. The duplicates were removed and
`src/rnmp/composition.py` now calls the original functions; the only thing kept
from the second attempt is `background_counts`, which is genuinely new (the
existing functions take background dicts as arguments and nothing computed them
from the reference). All 17 original metric tests pass again, and the figures are
byte-identical either way -- the two implementations agreed, which is why only
the test caught it.

**A reference detail this surfaced.** GRCh38 chrM carries a single `N` at
position 3107 -- the rCRS placeholder kept so GRCh38 chrM numbering matches rCRS.
It is excluded from every background, so the four base counts sum to 16,568 and
not 16,569. No rNMP is called there in any library (no read can align to an `N`,
and the reference-mismatch filter would drop it anyway), so it is harmless -- but
it is pinned by a test so the off-by-one is never "fixed" into a bug.

## 18. Figure 6A -- and a SECOND significance floor

The NR dinucleotide is the rNMP plus its 5' dNMP neighbour **on the rNMP's own
strand**: `pos-1` on the light strand, complement of `pos+1` on the heavy strand.
This is silent if wrong -- it produces a complete, plausible heatmap -- so the
orientation is taken from Suppl. Table S5's own column keys (`AA, CA, GA, TA, AC,
...`, i.e. the rNMP is the second character) and pinned by a test.

**Result.** C-C is the most enriched pattern for rCMP in nearly every library
(mean 0.352 against the 0.25 expectation), with T-A for rAMP (0.354), A-G for
rGMP (0.310) and C-U for rUMP (0.319). Only two of the four patterns per rNMP
exceed expectation in each case, and the liver biopsies separate from every other
cell type.

**The statistical finding.** Table S5 publishes p-values rather than
frequencies, so it validates the paper's conclusions, not its numbers. Calls
agree for **87 of 112** cell-type x pattern tests, and **45 of 48** wherever the
group has more than four libraries. All 25 disagreements run one way -- the paper
calls significance, this reproduction does not -- and 22 of them are in groups of
**N = 4**:

| test | smallest attainable one-tailed p at N = 4 |
|---|---|
| Wilcoxon signed-rank (correct one-sample test) | 1/16 = **0.0625** |
| Mann-Whitney U vs a constant (as the Methods state) | 1/5 = **0.2** |

The paper reports **p = 0.0105** for those patterns. Neither test can produce
that value from four observations, so the calls are unreachable rather than
merely different. Our frequencies at the disagreeing cells sit just above 0.25
(DLTB T-A at 0.364, C-C at 0.272) -- precisely where a low-powered test should
fail to resolve them, which is evidence that the frequencies are not the problem.

**This is the second instance of the same issue.** Section 10 (Finding 2) showed
the Figure 1A significance stars are likewise below the attainable floor at
N = 4 and N = 5. Two independent panels, one cause. Recorded, and per the
project's standing instruction not chased further: identifying which test
produces 0.0105 at N = 4 would require the authors' analysis scripts.

## 19. Fresh-clone test -- what it found

The repository was validated the way a stranger would encounter it: extract into
an empty directory, delete everything `.gitignore` excludes (leaving exactly what
`git clone` delivers), build both conda environments from the committed YAMLs,
run the tests, execute the notebook. No manual intervention.

**Result: passes.** 94 files, 9.0 MB. `pytest` -> **71 passed**. The notebook
executes **32 code cells, 10 inline figures, 6 rendered tables**, regenerating
all ten PNGs with no network access.

Seven defects were found and fixed. All seven were invisible from inside the
development directory, which is the point of the exercise.

**1. The notebook only ran from one directory.** Its setup cell computed
`ROOT = pathlib.Path.cwd().parent`, so the command documented in the README
(`python scripts/execute_notebook.py notebooks/...` from the repo root) resolved
ROOT to the repo's *parent* and failed. Fixed two ways: the notebook now finds its
root by walking up for `pyproject.toml`, and `execute_notebook.py` chdirs to the
notebook's own directory. Both invocations now behave identically.

**2. A clone contained none of the notebook's inputs.** `.gitignore` had
`results/**` and `data/reference/**`, so `chrM.fa`, both count sets and all 17
result tables were excluded -- and the README's claim that steps 1-3 need no
network was false for anyone cloning from GitHub. The count files are the boundary
between the expensive half of the pipeline (747 MB of downloads, hours of
alignment) and the cheap half (figures in ~10 s), so they are now committed: 3.4
MB gzipped, against the 11 GB of SAM/BED intermediates which stay ignored.

**3. Phase A could not run on a clone.** The pilot-validation section reads
FS310's `aligned.sam` from `results/ribosemap/`, correctly excluded as a pipeline
intermediate -- so the notebook died at the first analysis cell. FS310 is the
study's smallest library, and its alignment is 1.6 MB raw but **78 kB gzipped**,
with the coordinate BEDs another 114 kB. Those are now committed and
`pipeline.py` reads `.gz` transparently, falling back to a `.gz` sibling when the
plain file is absent. Two tests guard it, asserting the shipped file still
reproduces 4,866 aligned reads and 4,798 rNMPs.

**4. `fastqc` was a declared dependency that does not work.** The bioconda
package is a Java wrapper and the environment pulls in no JRE, so the executable
fails with "Unable to locate a Java Runtime". Nothing in the pipeline calls it --
Trim Galore only invokes it under `--fastqc`, which is never passed -- so it was
removed from `envs/rnmp-pipeline.yaml` with a comment saying what to add if QC
reports are wanted. Declaring a dependency that fails on first use is worse than
not declaring it.

**5. The reference directory shipped duplicates and two empty files.** Ribose-Map's
modules write byproducts into `data/reference/` during a run, and eight were being
committed. `chrM-chrM.fa` is a **byte-identical copy** of `chrM.fa` (16,852
bytes); `chrM-chromosomes.fa` is **0 bytes**, which is worse than absent because
an empty FASTA looks like data; the `.txt` unit lists are single-line files. None
is read anywhere in this repository. All are now ignored, and the tracked
reference is just the two FASTAs, their `.fai` indexes, the two `.chrom.sizes`,
the pinned Ribose-Map commit and the rotation offset.

**6. An empty gene annotation was masking its own failure.** `chrM_genes.bed` was
committed at 0 bytes -- the Ensembl GTF fetch in `00_setup_reference.sh` had
produced no MT records. The step guarded on `[[ ! -f ... ]]`, so on a fresh clone
the committed empty file would satisfy that test **forever** and the fetch would
never be retried. The guard is now `[[ ! -s ... ]]` (regenerate when absent *or*
empty), the step fails loudly and deletes the empty output instead of leaving it,
and the file is untracked. Nothing in the current analysis consumes it -- it was
added for transcription-overlap work that was never done -- so its absence
changes no result.

**7. Thirty-two regenerable config files were being committed.** `config/` held a
per-library Ribose-Map `.config` for each of the 32 libraries -- 845 bytes each,
identical apart from the library and run name. They are written by
`scripts/run_library.sh` at run time, never read back as inputs, and derived
deterministically from `config/samples_xu2024.tsv` and `config/params.yaml`
(both committed), so they carried no information the repository did not already
hold. Untracked; `config/` is now the three hand-authored files.

Also confirmed: `git init && git add -A && git commit` succeeds with 94 tracked
files and no empty directories (git does not track those, and one -- `analysis/`
-- was an empty leftover, now removed). Every pinned binary in
`envs/rnmp-pipeline.yaml` runs on osx-arm64: bowtie2 2.5.4, samtools 1.21,
bedtools 2.31.1, Trim Galore 0.6.10 with cutadapt 5.2, umi_tools 1.1.5, seqtk
1.4, mawk 1.3.4, sra-tools 3.2.1, R 4.3.3 with ggplot2 3.5.2 / ggseqlogo 0.2 /
data.table 1.17.8.

**Not covered by this test:** steps 4-6 (reference setup, alignment of the 32
libraries) were not re-run from scratch, because they need ~747 MB from ENA and a
GitHub clone of Ribose-Map. The environment that runs them is verified to solve
and every binary in it executes, but the end-to-end alignment path was last
exercised during the original cohort sweep, not in this clean directory.
