# Public-release checklist

Use this checklist before making the GitHub repository public or citing it in a manuscript.

## Scientific identity

- [ ] Repository name is `MotifSTaR`.
- [ ] Title is **MotifSTaR: allele-level reconstruction and analysis of short tandem repeat motif composition**.
- [ ] The manuscript, README, figures, DOI record, and software all use the same name and version.
- [ ] `python motifstar.py --version` reports the intended release version.
- [ ] The public commit corresponding to the manuscript is tagged, for example `v1.4.7`.

## Required legal decisions

- [ ] The code copyright owner has approved a software licence.
- [ ] A `LICENSE` file containing the exact approved licence is added.
- [ ] Institutional intellectual-property and publication requirements have been checked.
- [ ] Third-party catalogue/annotation licences and attributions are documented.
- [ ] Redistribution of every example SVG is explicitly permitted.

The repository currently does not contain a software licence because that legal choice must be made by the copyright owner. Until a licence is added, public readers can view the code but do not automatically receive permission to reuse it. This must be resolved before the manuscript says the software is freely reusable.

## Privacy and provenance

- [ ] No private participant SVGs are included.
- [ ] No private sample IDs are present in tables, workbooks, figures, manifests, or logs.
- [ ] No absolute `/g/data`, `/scratch`, home-directory, or Windows user paths are present.
- [ ] No credentials, access tokens, email secrets, or internal URLs are present.
- [ ] Public example data record source URL, version, date, genome build, and licence.
- [ ] The catalogue and annotation record source, modifications, reviewer, and hashes.

Useful pre-release searches from the repository root:

```bash
git grep -n -E '/g/data|/scratch|C:\\Users|token|password|secret'
git status --short
```

Inspect binary XLSX files manually as text searches will not cover their contents reliably.

## Reproducibility

- [ ] A clean environment can be created from `environment.yml`.
- [ ] `python -m py_compile motifstar.py` succeeds.
- [ ] All regression tests pass.
- [ ] `--validate-only` succeeds for the public example.
- [ ] The documented example command completes from a fresh output directory.
- [ ] Example figures open correctly in a browser/PDF viewer.
- [ ] `run_manifest.json` matches the released commit, configuration, and public example inputs.
- [ ] Figure-source CSVs are included with selected example figures.

## Validation and claims

- [ ] Validation numbers are linked to the exact public commit and configuration hashes.
- [ ] Callability is reported with conditional/adjudicated concordance.
- [ ] The 100% result is described as adjudicated agreement among complete SVG-derived calls, not overall clinical accuracy.
- [ ] Runtime is labelled as core processing time or end-to-end wall time, not a mixture of both.
- [ ] Any version change after validation has been regression-tested and disclosed.

## GitHub publication

Create an empty public repository named `MotifSTaR`, without adding a second README or `.gitignore` in the GitHub form. Then, from this local folder:

```bash
git init
git add .
git commit -m "Prepare MotifSTaR 1.4.7 public release"
git branch -M main
git remote add origin https://github.com/Mamunur-20/MotifSTaR.git
git push -u origin main
```

Before `git add .`, recheck the privacy and licence items above.

After pushing:

- [ ] Add the short description and topics from `docs/PROJECT_DESCRIPTION.md`.
- [ ] Enable GitHub Issues if public bug reports are wanted.
- [ ] Create a signed or annotated version tag and GitHub Release.
- [ ] Connect the repository to Zenodo and archive the tagged release.
- [ ] Add the final GitHub and Zenodo DOI links to the manuscript and `CITATION.cff`.
- [ ] Protect or retain the exact release for the journal's required availability period.

## Suggested licence discussion

Ask the copyright owner or institutional research office whether the intended release should use a standard open-source licence such as MIT, BSD-3-Clause, Apache-2.0, or GPL-3.0, or whether a separately approved non-commercial licence is required. Do not describe the repository as open source if the selected licence prohibits uses that do not meet an open-source definition.
