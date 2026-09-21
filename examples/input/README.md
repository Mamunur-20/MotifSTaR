# Example input

Add a small, redistributable set of ExpansionHunter REViewer `.svg` or `.svg.gz` files under `svgs/`.

For every example dataset, record:

- source repository or dataset URL;
- dataset/release version and download date;
- licence or redistribution permission;
- reference genome and ExpansionHunter/REViewer version, when known;
- any renaming, filtering, or de-identification performed;
- the expected sample and locus identifiers.

Do not add private participant data. Git ignores the SVG directory by default; after checking every file, use `git add -f examples/input/svgs/<reviewed-file>` to add only the public examples you intend to release.
