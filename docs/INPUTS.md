# Inputs

MotifSTaR requires three inputs and accepts two optional mapping tables.

## Required inputs

### 1. REViewer SVG root (`--svg-root`)

The directory is searched recursively for `.svg` and `.svg.gz` files. The usual layout is:

```text
svg_root/
  batch_or_collection/
    sample_001/
      sample_001.eh_realigned.reviewer.ZFHX3.svg
    sample_002/
      sample_002.eh_realigned.reviewer.ZFHX3.svg.gz
```

The normal filename suffix is `.reviewer.<LOCUS>.svg` or `.reviewer.<LOCUS>.svg.gz`. MotifSTaR derives the locus from that suffix and normally uses the parent directory as `sampleID`. The directory above the sample can be retained as `batch` provenance.

Use `--file-map` when these assumptions do not describe the data. Run `--validate-only` before analysis to inspect how every file was identified.

### 2. ExpansionHunter variant catalogue (`--catalog`)

Provide a JSON list of locus objects. MotifSTaR uses fields including:

- `LocusId` and, when present, `Gene`;
- `LocusStructure`;
- `ReferenceRegion`;
- `RepeatUnit`;
- `VariantId` and `VariantType` for compound loci;
- `MainReferenceRegion` to identify the primary biological target;
- optional motif and disease annotations used for QC or descriptive figures.

Catalogue fields may be scalar values or arrays, depending on whether a locus has one or several components.

For the recommended biological report, use `--component-mode primary`. This reports the catalogue `RepeatUnit` for the primary region and avoids attaching flanking repeat components to the main allele structure.

### 3. Reporting-orientation annotation (`--annotation`)

Provide an `.xlsx`, `.xlsm`, `.csv`, or `.tsv` table containing a locus/gene identifier and strand or reporting orientation. The first workbook sheet is used unless `--annotation-sheet` is supplied.

Minimum example:

| Gene | strand |
|---|---|
| ZFHX3 | + |
| C9orf72 | - |

Common accepted identifier headings include `LocusId`, `locus`, `Gene`, `gene_symbol`, and `symbol`. Strand values can include `+`, `-`, `forward`, `reverse`, `plus`, `minus`, `sense`, and `antisense`.

An explicit `reporting_orientation` column may contain `reference` or `reverse_complement`. A `reported_repeat_unit` column may be used only when a curated reporting motif must override the oriented catalogue motif. Unresolved orientation is withheld by default rather than guessed.

## Optional inputs

### File map (`--file-map`)

Use a CSV, TSV, or XLSX file to assign identities when the directory or filename convention is unsuitable.

```csv
svg_path,sampleID,locus,batch
svgs/public_001.svg,public_001,ZFHX3,public_demo
svgs/public_002.svg.gz,public_002,ZFHX3,public_demo
```

Accepted path headings include `svg_path`, `path`, `file`, and `filename`. Accepted sample headings include `sampleID`, `sample_id`, and `sample`. Relative paths are resolved under `--svg-root`.

### Genotype override (`--genotypes`)

Use a CSV, TSV, or XLSX file containing `sampleID`, `locus`, and either a combined genotype or two allele-size columns.

Combined form:

```csv
sampleID,locus,genotype
public_001,ZFHX3,18|21
```

Two-column form:

```csv
sampleID,locus,allele1_size,allele2_size
public_001,ZFHX3,18,21
```

Accepted genotype headings include `genotype` and `eh_genotype`.

## Public example data

The repository deliberately does not contain private cohort SVGs. Before adding a public example:

1. Confirm that redistribution is permitted by the source licence or data-access terms.
2. Remove or replace participant identifiers and private paths.
3. Include the source, download date, dataset version, and licence in `examples/input/README.md`.
4. Keep the example small enough for a normal Git clone.
5. Run MotifSTaR on exactly those files and include the corresponding example output, configuration, and manifest.

Git ignores `examples/input/svgs/` and generated `examples/output/` by default to reduce the risk of accidental data publication. Add reviewed public files explicitly with `git add -f`.

## Input audit

Use a fresh output directory:

```bash
python -u motifstar.py \
  --svg-root examples/input/svgs \
  --catalog config/variant_catalog_without_offtargets.GRCh38.json \
  --annotation config/gene_strands_companion_locus.xlsx \
  --output examples/output/input_check \
  --component-mode primary \
  --validate-only
```

Review `input_validation.csv`, especially:

- files whose sample or locus could not be parsed;
- loci missing from the catalogue;
- loci missing a reporting orientation;
- unexpected duplicate sample–locus combinations;
- catalogue/SVG component mismatches.

Add `--strict` if an automated workflow should return a non-zero exit code when errors or unconfigured loci are found.
