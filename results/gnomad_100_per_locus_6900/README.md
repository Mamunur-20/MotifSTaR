\# MotifSTaR results for the gnomAD STR population-reference dataset



This directory contains MotifSTaR outputs derived from 6,900 publicly

available ExpansionHunter REViewer SVGs: 100 SVGs per locus across 69 STR

loci. These are 6,900 sample–locus observations and should not be interpreted

as 6,900 confirmed unique participants.



\## Run summary



\- MotifSTaR version: 1.4.7

\- Input SVGs: 6,900

\- Observed loci: 69

\- Unconfigured loci: 0

\- Component calls: 6,900

\- PASS: 5,617

\- REVIEW: 1,283

\- ERROR: 0

\- Workers: 4

\- Processing time: 136.1 seconds



REVIEW is a quality-control outcome indicating that one or more configured

evidence requirements were not satisfied. It is not a processing error.



\## Input data



The source SVGs were obtained from the public gnomAD STR resource. The complete

input archive is not stored in this Git repository.



Input archive:

https://drive.google.com/file/d/15PRshfYcfhfw19yqNvKGCMIB7pOLEiDr/view?usp=sharing



Public source:

https://storage.googleapis.com/gnomad-str-public/



The gnomAD terms of use and citation requirements apply to the source data.

This dataset is provided as a population-reference software demonstration and

should not be described as a cohort of neurologically healthy controls.



\## Contents



\- `all\_loci.calls.csv`: complete sample–locus call table

\- `all\_loci.calls.xlsx`: formatted workbook for review

\- `all\_loci.qc.csv`: detailed evidence and QC table

\- `calls\_by\_locus/`: locus-specific call tables

\- `figure\_data/`: data underlying the generated figures

\- `figures/`: cohort-level and locus-level figures

\- `read\_level/read\_evidence.tsv.gz`: compressed read-level evidence

\- `run\_manifest.json`: software version, parameters and input hashes

\- `input\_validation.csv`: catalogue and annotation audit

