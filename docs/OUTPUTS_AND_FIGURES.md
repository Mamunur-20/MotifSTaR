# Outputs and figures

MotifSTaR writes results only after the input audit and calling steps have completed. Use a new output directory for each scientific run so that the table, manifest, and figures always describe the same analysis.

## Output tree

```text
output_directory/
├── all_loci.calls.csv
├── all_loci.calls.xlsx
├── all_loci.qc.csv
├── input_validation.csv
├── run_manifest.json
├── run.log
├── calls_by_locus/
│   └── <locus>.calls.csv
├── read_level/
│   └── read_evidence.tsv.gz
├── figure_data/
│   ├── allele_level.csv
│   ├── locus_summary.csv
│   ├── interruption_motif_spectrum.csv
│   └── sequence_composition_frequency.csv
└── figures/
    ├── cohort_interruption_prevalence.<format>
    ├── cohort_qc_overview.<format>
    ├── cohort_interruption_motif_spectrum.<format>
    └── locus/
        ├── <locus>.summary.<format>
        └── <locus>.sequence_composition_frequency.<format>
```

Some files are omitted when the relevant feature is disabled or when there is no data for a particular plot.

## Call tables

### `all_loci.calls.csv`

This is the main machine-readable result. It contains batch, sample, locus/component identity, allele structures, reconstructed sizes, read-support fields, source SVG, and QC status/reasons.

The most commonly used fields are:

| Field | Meaning |
|---|---|
| `batch` | Input collection/batch provenance derived from the layout or file map. |
| `sampleID` | Sample identifier used by MotifSTaR. |
| `locus` | Reported locus or component identifier. |
| `parent_locus` | Catalogue locus containing the selected component. |
| `locus_structure` | Catalogue structure supplied for the locus. |
| `allele1`, `allele2` | Run-length encoded motif paths, or `NA`/an explicit sentinel when no structure is reported. |
| `allele1_size`, `allele2_size` | Total reconstructed repeat units in each reported motif path. |
| `read_depth_allele1`, `read_depth_allele2` | Eligible reads matching the final reported position-wise consensus structure. |
| `consensus_eligible_read_depth_alleleN` | Clean, full-spanning reads of the expected allele length that entered the relevant consensus/grouping step. |
| `full_spanning_read_depth_alleleN` | Rows that geometrically span the left flank, repeat, and right flank. |
| `qc_status` | `PASS`, `REVIEW`, or `ERROR`. |
| `primary_qc_flag` | One principal machine-readable reason selected for summary. |
| `review_reason` | Plain-language explanation for review. |
| `source_svg` | Input file used for the result. Check for private absolute paths before sharing a table publicly. |

The compact files under `calls_by_locus/` contain the seven main report fields:

```text
sampleID,allele1,allele2,allele1_size,allele2_size,read_depth_allele1,read_depth_allele2
```

### `all_loci.calls.xlsx`

This workbook is designed for assisted manual review. It presents the calls, flag details, and a flag legend with stable colours. Excel formatting is not preserved in CSV, so use the workbook when reviewers need filtering and visual flag cues.

### `all_loci.qc.csv`

This is the detailed evidence/QC table. It preserves component mapping, genotype source, orientation, evidence scopes, clustering counts, interruption thresholds, every QC flag, and technical notes. Keep it with the main calls table even if it is not part of a figure.

## Read-depth hierarchy

MotifSTaR distinguishes several denominators:

1. **Displayed rows:** every read row rendered in the relevant SVG evidence scope.
2. **Full-spanning rows:** rows that traverse the repeat and required left/right flanks.
3. **Callable full-spanning rows:** full-spanning rows whose target sequence is unambiguous and divisible into complete motif units.
4. **Consensus-eligible rows:** callable rows that also have the expected allele length and, for same-sized calls, can enter the selected sequence group.
5. **Supporting rows:** eligible rows matching the final reported structure.

These are not interchangeable. Report the named denominator with any support percentage.

## QC status

- **PASS:** no substantive ambiguity was detected. Informational provenance flags may still be present.
- **REVIEW:** the call requires inspection because evidence was insufficient or ambiguous, an expected allele was missing, an indel dominated the target, orientation/component mapping was unresolved, or another configured review condition occurred.
- **ERROR:** the SVG could not be parsed or processed.

A withheld allele is not automatically an incorrect motif prediction. It is an explicit no-call under the configured evidence requirements.

## Read-level evidence

`read_level/read_evidence.tsv.gz` records each displayed row, reconstructed sequence/motif path when possible, panel/component context, eligibility status, and exclusion reason. It supports audit and debugging but can be large. Disable it with `--no-write-read-level` only when storage is constrained and row-level auditability is not required.

## Reproducibility files

### `run_manifest.json`

The manifest records:

- MotifSTaR and Python versions;
- complete threshold/configuration values;
- resolved input paths;
- SHA-256 hashes of catalogue and annotation files;
- selected loci and number of workers;
- input, output, QC, read-evidence, figure, and elapsed-time counts;
- orientation, compound-pooling, and allele-order policy notes.

### `input_validation.csv`

This file shows how each input was configured and highlights unconfigured loci or identity problems.

### `run.log`

The log records progress, warnings, errors, and final call counts. A first-time Matplotlib font-cache message is expected on some systems.

## How to read the figures

### Cohort interruption prevalence

Each locus is summarized by the proportion of resolved alleles containing at least one motif other than the catalogue repeat unit. Error bars are Wilson 95% confidence intervals. The denominator is resolved alleles, not all enrolled participants; check callability alongside prevalence.

### Cohort QC overview

The callability panel shows the fraction of expected alleles for which MotifSTaR reported a structure. The support panel shows median supporting full-spanning read depth. A low callability value indicates that many allele structures were withheld under the QC rules; it does not by itself imply a high false-call rate.

### Cohort interruption motif spectrum

This plot counts non-catalogue motif occurrences across the reported allele motif paths. For example, a label such as `TBP — ACA: 1308` means the motif token `ACA` appeared 1,308 times across the resolved allele paths included in the plot. It does **not** mean 1,308 participants and does not mean 1,308 independent interruption-positive alleles unless that exact denominator is calculated separately.

### Per-locus summary

Panel A displays reconstructed repeat length in repeat units. A **pure** path contains only the catalogue repeat unit; an **interrupted** path contains at least one different motif token. Catalogue normal-maximum or pathogenic-minimum lines, when available, provide descriptive context only and are not standalone clinical classifications.

Panel B lists the most frequent exact resolved motif paths. Bar length is the number of resolved alleles with that complete path, not the number of participants.

### Sequence-composition frequency

Each row is one exact resolved motif path. The text gives the ordered run-length encoded path; the motif charts summarize catalogue and interruption units; repeat units give total reconstructed length; and the horizontal bar gives the number of resolved alleles with that exact path.

## Adding public example output

Generate a fresh run from only the public example SVGs, then review every file before committing it. In particular:

- replace or remove identifying sample IDs;
- remove private absolute paths from tables, logs, and the manifest;
- verify that no source SVG contains sensitive metadata;
- include `run_manifest.json` and `figure_data/` with the figures;
- state the MotifSTaR version and command in `examples/output/README.md`.

Because `.gitignore` excludes generated output by default, add only the reviewed public files explicitly with `git add -f`.
