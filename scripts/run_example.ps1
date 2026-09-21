$ErrorActionPreference = "Stop"
$repoDir = Split-Path -Parent $PSScriptRoot
Set-Location $repoDir

python -u motifstar.py `
  --svg-root examples/input/svgs `
  --catalog config/variant_catalog_without_offtargets.GRCh38.json `
  --annotation config/gene_strands_companion_locus.xlsx `
  --output examples/output/motifstar_run `
  --workers 4 `
  --component-mode primary `
  --min-spanning-reads 3 `
  --min-cluster-reads 3 `
  --heterozygous-interruption-min-fraction 0.50 `
  --same-size-min-cluster-reads 2 `
  --same-size-min-cluster-fraction 0.20 `
  --same-size-min-top-two-coverage 0.90 `
  --same-size-within-cluster-consensus-min-fraction 0.50 `
  --same-size-max-reported-clusters 2 `
  --same-size-max-assignment-distance 2 `
  --figure-formats svg,pdf,png
