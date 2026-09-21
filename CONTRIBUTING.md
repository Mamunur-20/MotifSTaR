# Contributing to MotifSTaR

Thank you for helping improve MotifSTaR.

## Before opening an issue

1. Confirm the problem occurs with the current release.
2. Run `python motifstar.py --version` and include the reported version.
3. Check `run.log`, `input_validation.csv`, `all_loci.qc.csv`, and `run_manifest.json`.
4. Remove participant identifiers, private filesystem paths, and other sensitive information.

## A useful bug report includes

- a short description of the expected and observed behaviour;
- the complete command, with private paths replaced by neutral placeholders;
- the relevant QC flags and review reason;
- a minimal redistributable SVG, catalogue entry, and annotation row when possible;
- operating system, Python version, and package versions.

## Code changes

Create a branch, keep the change focused, add or update a regression test, and run:

```bash
python -m unittest discover -s tests -p "test_*.py" -v
python -m py_compile motifstar.py
```

Do not change a scientific default without documenting the reason, expected effect, validation evidence, and backward-compatibility impact. Generated calls should be compared at allele, QC, and read-support levels.

## Sensitive data

Never commit private participant SVGs, manually curated clinical tables, absolute institutional paths, credentials, or access tokens. Public test data must have documented redistribution permission.
