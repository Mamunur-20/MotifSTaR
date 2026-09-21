# Command-line reference

Run `python motifstar.py --help` to see the options available in the installed version. This page explains what each public option means and records the version 1.4.7 defaults.

## Recommended reproducible command

```bash
python -u motifstar.py \
  --svg-root examples/input/svgs \
  --catalog config/variant_catalog_without_offtargets.GRCh38.json \
  --annotation config/gene_strands_companion_locus.xlsx \
  --output examples/output/motifstar_run \
  --workers 4 \
  --component-mode primary \
  --min-spanning-reads 3 \
  --min-cluster-reads 3 \
  --heterozygous-interruption-min-fraction 0.50 \
  --same-size-min-cluster-reads 2 \
  --same-size-min-cluster-fraction 0.20 \
  --same-size-min-top-two-coverage 0.90 \
  --same-size-within-cluster-consensus-min-fraction 0.50 \
  --same-size-max-reported-clusters 2 \
  --same-size-max-assignment-distance 2 \
  --figure-formats svg,pdf,png
```

`python -u` requests unbuffered output so progress messages appear promptly, which is useful in HPC logs. It does not change the scientific analysis.

## Required paths

| Option | Meaning |
|---|---|
| `--svg-root PATH` | Root directory searched recursively for `.svg` and `.svg.gz` REViewer files. |
| `--catalog PATH` | ExpansionHunter variant catalogue JSON. |
| `--annotation PATH` | XLSX/XLSM/CSV/TSV companion table containing locus/gene and reporting strand/orientation. |
| `--output PATH` | Destination directory. It should be new or empty unless `--overwrite` is deliberately used. |

## Input interpretation

| Option | Default | Meaning |
|---|---:|---|
| `--annotation-sheet NAME` | first sheet | Selects a worksheet in an XLSX/XLSM annotation file. |
| `--genotypes PATH` | none | Supplies curated genotype/allele-size overrides. This changes the expected allele-size branch but does not provide motif composition. |
| `--file-map PATH` | none | Supplies explicit SVG path, sample, locus, and optional batch metadata when filenames/folders do not follow the normal convention. |
| `--locus LOCUS` | all | Limits processing to one or more loci/components. Repeat the option or use comma-separated values, for example `--locus ATXN7 --locus FXN` or `--locus ATXN7,FXN`. |
| `--component-mode {combined,primary,all}` | `combined` | Controls compound loci. `combined` concatenates catalogue repeat components into one allele path; `primary` reports only the primary catalogue target; `all` emits separate component records for diagnostic work. `primary` is recommended for the biological report used in this project. |

### Component-mode warning

The program retains `combined` as its legacy default. A production command should state the intended mode explicitly. In the supplied analysis, `--component-mode primary` means that the reported allele follows the catalogue `RepeatUnit` at the primary reference region rather than including adjacent repeat blocks.

## Processing and general support

| Option | Default | Meaning |
|---|---:|---|
| `--workers N` | scheduler allocation or at most 4 | Number of SVG files processed in parallel. If omitted, MotifSTaR uses `PBS_NCPUS`, then `SLURM_CPUS_PER_TASK`, or otherwise `min(4, available CPUs)`. Do not exceed the CPUs allocated by a scheduler. |
| `--min-spanning-reads N` | `3` | Minimum clean full-spanning support needed to avoid an insufficient-depth review condition and to report a supported allele structure. |
| `--min-cluster-reads N` | `3` | Minimum reads supporting a motif path for standard strong-cluster/allele-support assessment. |
| `--max-near-motif-distance N` | `0` | Maximum motif-unit Hamming distance used by optional near-path clustering. `0` keeps exact complete paths and is recommended unless an alternative has been validated. |
| `--near-cluster-parent-ratio X` | `3.0` | Minimum parent-to-minor support ratio for near-path assignment when near clustering is enabled. It has no effect on the recommended exact-clustering setting of distance `0`. |
| `--strong-cluster-fraction X` | `0.0` | Optional minimum fraction for a standard strong cluster in addition to its read count. `0` disables this extra requirement. |

## Different-sized allele consensus

| Option | Default | Meaning |
|---|---:|---|
| `--heterozygous-interruption-min-fraction X` | `0.50` | Within a different-sized panel, a particular non-reference nucleotide is incorporated at a repeat position only when its support fraction is **strictly greater than** this threshold. At `0.50`, exactly 50% is not enough. |

Different-sized panels are normally evaluated separately. MotifSTaR also compares the complete ordered vector of displayed repeat-block sizes at compound loci. Panels are pooled as same-sized evidence only when those complete vectors are equal, even if the selected primary blocks happen to have the same size.

## Same-sized allele inference

| Option | Default | Meaning |
|---|---:|---|
| `--same-size-min-cluster-reads N` | `2` | Minimum exact-path reads required for a candidate same-sized allele group. |
| `--same-size-min-cluster-fraction X` | `0.20` | A candidate group must also represent at least this fraction of eligible same-sized evidence. |
| `--same-size-min-top-two-coverage X` | `0.90` | Selected same-sized group(s) must collectively explain at least this fraction of eligible evidence. Below the threshold, the result is retained for review or withheld according to the available support. |
| `--same-size-within-cluster-consensus-min-fraction X` | `0.50` | Within each selected group, a non-reference nucleotide is incorporated only when its support is strictly greater than this fraction. |
| `--same-size-max-reported-clusters {1,2}` | `2` | Maximum number of supported same-sized sequence structures that can be reported. |
| `--same-size-max-assignment-distance N` | `2` | Maximum motif-unit differences for assigning a minor path to one uniquely nearest selected group. Ambiguous nearest assignments are not forced. |

REViewer SVGs do not supply read identifiers for molecular phasing. Same-sized panels are therefore pooled for sequence-group inference. When distinct groups have strong, unique panel-specific support, MotifSTaR preserves top/bottom panel order. Otherwise it uses a stable reporting order without claiming parental or physical haplotype phase.

## Deduplication and strand policy

| Option | Default | Meaning |
|---|---:|---|
| `--homozygous-deduplication {none,cross_panel_exact}` | `none` | Controls an optional exact-row heuristic across panels. Because SVGs have no read IDs, the safe default does not claim that identical displayed rows are duplicate molecules. |
| `--unresolved-strand-policy {skip,technical,error}` | `skip` | `skip` withholds a biological motif call when reporting orientation is unknown; `technical` reports in technical SVG/reference orientation; `error` treats the condition as an error. `skip` is the recommended conservative policy. |

## Output controls

| Option | Default | Meaning |
|---|---:|---|
| `--write-read-level` / `--no-write-read-level` | enabled | Enables or disables `read_level/read_evidence.tsv.gz`. Disable only when the storage saving is more important than row-level auditability. |
| `--figures` / `--no-figures` | enabled | Enables or disables all figure generation and figure-source tables. |
| `--per-locus-figures` / `--no-per-locus-figures` | enabled | Enables or disables the per-locus figures while retaining cohort figures. Disabling them can reduce runtime and file count for large panels. |
| `--figure-formats LIST` | `png,pdf,svg` | Comma-separated formats chosen from `png`, `pdf`, and `svg`, for example `--figure-formats svg` or `--figure-formats svg,pdf,png`. |
| `--dpi N` | `300` | Resolution for raster PNG output. It does not change vector SVG/PDF geometry. |

## Validation, exit behaviour, and overwrite safety

| Option | Default | Meaning |
|---|---:|---|
| `--validate-only` | off | Audits files and configuration, writes validation/manifest information, and stops before SVG calling. Use this before every new dataset. |
| `--strict` | off | Returns a non-zero exit status when errors or unconfigured loci are present. Useful in automated pipelines. |
| `--fail-on-review` | off | Returns exit code 3 when any result has `qc_status=REVIEW`. The result files are still written. |
| `--overwrite` | off | Allows reuse of a non-empty output directory and removes only recognized MotifSTaR-owned outputs before writing. A new output directory is safer for reproducibility. |
| `--allow-duplicate-sample-locus` | off | Explicitly permits the same sample–locus key in multiple batches. Use only when those repeated observations are intentional and batch provenance will be retained. |
| `--log-level {DEBUG,INFO,WARNING,ERROR}` | `INFO` | Controls log detail. `DEBUG` can assist development but may create a large log. |
| `--version` | — | Prints the software version and exits. |
| `-h`, `--help` | — | Prints command-line help and exits. |

## Practical command variants

### One locus

```bash
python -u motifstar.py \
  --svg-root examples/input/svgs \
  --catalog config/variant_catalog_without_offtargets.GRCh38.json \
  --annotation config/gene_strands_companion_locus.xlsx \
  --output examples/output/zfhx3_only \
  --locus ZFHX3 \
  --component-mode primary \
  --workers 2 \
  --figure-formats svg
```

### Cohort tables without figures or read-level output

```bash
python -u motifstar.py \
  --svg-root /path/to/svgs \
  --catalog config/variant_catalog_without_offtargets.GRCh38.json \
  --annotation config/gene_strands_companion_locus.xlsx \
  --output /path/to/results/motifstar_tables \
  --component-mode primary \
  --workers 4 \
  --no-figures \
  --no-write-read-level
```

### Cohort figures only, without per-locus plots

```bash
python -u motifstar.py \
  --svg-root /path/to/svgs \
  --catalog config/variant_catalog_without_offtargets.GRCh38.json \
  --annotation config/gene_strands_companion_locus.xlsx \
  --output /path/to/results/motifstar_cohort \
  --component-mode primary \
  --workers 4 \
  --no-per-locus-figures \
  --figure-formats svg,pdf
```

## Exit codes

- `0`: run completed under the selected policy.
- `1`: fatal input, configuration, processing, or unexpected error.
- `2`: strict-mode failure caused by errors or unconfigured loci.
- `3`: `--fail-on-review` was used and at least one review call was present.

Always read `run.log` and `run_manifest.json`; an exit code alone does not explain biological callability.
