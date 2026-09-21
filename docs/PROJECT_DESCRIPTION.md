# Project title and abstract

## Title

**MotifSTaR: allele-level reconstruction and analysis of short tandem repeat motif composition**

## Abstract

Repeat expansions are established causes of neurological disease, but routine short-read workflows commonly report allele size while leaving nucleotide-level repeat composition to visual interpretation. We present MotifSTaR, an original tool that converts ExpansionHunter REViewer SVGs into strand-normalized, allele-level motif structures using catalogue-aware repeat selection, fully spanning reads, genotype-aware consensus and explicit quality-control rules. MotifSTaR was evaluated using two complementary approaches. First, three reviewers manually examined 1,600 sample–locus outputs from 20 samples across 80 configured STR loci representing 77 gene annotations. Reviewer consensus indicated that all reported structures were supported by the corresponding SVG evidence. Second, a quantitative comparison was performed using 1,000 manually annotated *ZFHX3* cases comprising 800 MinE controls and 200 patient samples. MotifSTaR generated complete two-allele structures for 915 cases. Initial manual annotations agreed exactly with 871 (95.2%) of these complete calls. Manual re-review of all 44 discordant calls supported the MotifSTaR result in every case. The tool therefore achieved 100% concordance with the adjudicated SVG interpretation among complete calls, with an overall two-allele callability of 91.5% (915/1,000). In a four-worker high-performance-computing run, MotifSTaR processed 4,000 SVGs and generated analysis-ready tables, quality-control outputs and figures in 142 seconds. These validation results demonstrate concordance with expert interpretation of SVG evidence.

**Availability and implementation:** MotifSTaR is implemented in Python 3.10 or later as a command-line application and can run on a standalone computer or in a Linux-based high-performance-computing environment with the required Python packages. The required inputs are REViewer SVG/SVG.GZ files, an ExpansionHunter JSON variant catalogue and an XLSX/CSV strand-annotation table. Source code, documentation, tests, a Conda environment specification and archived releases will be available from the project repository and DOI-linked archive.

## Short GitHub description

Reconstruct strand-normalized, allele-level STR motif composition from ExpansionHunter REViewer SVGs, with explicit read support, QC, assisted-review tables, and publication-ready figures.

## Suggested GitHub topics

`short-tandem-repeats`, `repeat-expansions`, `expansionhunter`, `reviewer`, `motif-interruptions`, `bioinformatics`, `genomics`, `python`
