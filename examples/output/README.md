# example output

Generate this directory from the redistributable SVGs in `examples/input/svgs/` using the documented recommended command.

Before committing output:

1. confirm that all sample IDs are public-safe;
2. remove or replace private absolute paths in CSV, XLSX, JSON, and log files;
3. retain the exact command, MotifSTaR version, catalogue/annotation hashes, and `run_manifest.json`;
4. retain `figure_data/` with every published visualization;
5. state the SVG source and redistribution terms;
6. open the workbook and figures to confirm that they render correctly.

Generated files are ignored by default. Add reviewed files explicitly with `git add -f` rather than removing the safety rule globally.
