# Validation and benchmark provenance

This document separates the two validation strategies and makes their denominators explicit.

## Broad multi-locus review

Three reviewers examined 1,600 sample–locus outputs from 20 samples across 80 configured STR loci representing 77 gene annotations. The difference between 80 loci and 77 gene annotations is expected because some genes have more than one configured locus. Reviewer consensus indicated that all reported structures were supported by the corresponding REViewer SVG evidence.

This review assesses whether the MotifSTaR-reported structure is supported by the displayed evidence. It is not a measurement of clinical sensitivity or specificity.

## Quantitative `ZFHX3` comparison

The quantitative comparison contained 1,000 manually annotated cases:

| Dataset | Cases | Complete two-allele MotifSTaR calls | Initial exact agreement among complete calls | Adjudicated agreement among complete calls |
|---|---:|---:|---:|---:|
| MinE controls | 800 | 741 (92.6%) | 706/741 (95.3%) | 741/741 (100%) |
| Patient subset | 200 | 174 (87.0%) | 165/174 (94.8%) | 174/174 (100%) |
| Combined | 1,000 | 915 (91.5%) | 871/915 (95.2%) | 915/915 (100%) |

An initial exact match required both allele motif structures and reconstructed sizes to agree after allowing allele-order exchange in the historical comparison. Forty-four complete discordant cases were manually re-reviewed; all 44 supported the tool interpretation. Eighty-five cases did not provide complete two-allele calls under the configured evidence requirements and remained in the overall callability denominator.

The appropriate headline measures are therefore:

- two-allele callability: 915/1,000 = 91.5%;
- initial exact agreement among complete calls: 871/915 = 95.2%;
- adjudicated agreement among complete calls: 915/915 = 100%;
- complete-and-adjudicated-supported two-allele yield across all cases: 915/1,000 = 91.5%.

The 100% value must always be described as adjudicated agreement among complete calls, not as overall biological accuracy.

## Performance

A reported four-worker Linux HPC run processed 4,000 REViewer SVGs and generated tables, QC outputs, and figures in 142 seconds. The supplied 1.4.7 result manifest records a separate four-worker run of 4,000 SVGs with 84.414 seconds of core processing time; its full wall-clock duration was approximately 140 seconds when input audit and figure/output generation are included.

For a paper or release note, state exactly which timer is used:

- **core processing time** from `run_manifest.json`; or
- **end-to-end wall-clock time** from the first and final log timestamps.

Do not mix the two values. Hardware, filesystem, enabled read-level output, per-locus figures, and figure formats can materially change runtime.

## Version attribution

The public release is version 1.4.7. The 1.4.7 update changed compound-panel pooling and top/bottom allele handling relative to earlier development versions. Before the manuscript attributes validation numbers to the public release, archive:

1. the exact 1.4.7 code commit;
2. the catalogue and annotation hashes;
3. the complete validation command and run manifest;
4. the adjudication table and reviewer procedure;
5. the exact version-specific validation output.

## Interpretation boundary

Both validation strategies use expert interpretation of REViewer SVG evidence. They demonstrate reproducibility and concordance with that interpretation. An orthogonal sequencing assay would be required to estimate analytical sensitivity and specificity independently of the SVG representation.
