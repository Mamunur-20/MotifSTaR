# Publishing MotifSTaR on GitHub

Complete `RELEASE_CHECKLIST.md` before making the repository public, especially the licence, third-party data, and privacy checks.

## Recommended repository details

- **Repository name:** `MotifSTaR`
- **Description:** Reconstruct strand-normalized, allele-level STR motif composition from ExpansionHunter REViewer SVGs, with explicit read support, QC, assisted-review tables, and publication-ready figures.
- **Visibility:** Public, once the code/data owners have approved release.
- **Topics:** `short-tandem-repeats`, `repeat-expansions`, `expansionhunter`, `reviewer`, `motif-interruptions`, `bioinformatics`, `genomics`, `python`

## Option A: upload with Git (recommended)

### 1. Create the GitHub repository

1. Sign in to GitHub.
2. Select **New repository**.
3. Choose the owner `Mamunur-20` and repository name `MotifSTaR`.
4. Do not ask GitHub to add a README, `.gitignore`, or licence at this stage, because the local folder already contains the first two and the licence still needs owner approval.
5. Create the repository.

### 2. Review the local folder

Open a terminal in the prepared `MotifSTaR` folder and check that no sensitive files are present:

```bash
git grep -n -E '/g/data|/scratch|C:\\Users|token|password|secret' -- . || true
```

This search cannot reliably inspect binary XLSX files, so open and review the companion workbook manually as well.

### 3. Create the first commit

```bash
git init
git add .
git status
git commit -m "Prepare MotifSTaR 1.4.7 public release"
git branch -M main
```

Read the `git status` list before committing. The private full-cohort output is not inside this prepared repository and should remain outside it.

### 4. Connect and push

```bash
git remote add origin https://github.com/Mamunur-20/MotifSTaR.git
git push -u origin main
```

If GitHub requests authentication, use the approved browser/device flow, GitHub CLI, SSH key, or a personal access token. Do not put a token in a command saved in shell history or in any repository file.

### 5. Tag the manuscript version

After tests and release checks pass:

```bash
git tag -a v1.4.7 -m "MotifSTaR 1.4.7"
git push origin v1.4.7
```

Create a GitHub Release from that tag and describe the configuration/validation version clearly.

## Option B: GitHub browser upload

The browser is suitable for a first small upload, but Git is safer for later updates.

1. Create the empty `MotifSTaR` repository.
2. Select **uploading an existing file**.
3. Drag the **contents** of the local `MotifSTaR` folder into the upload area. Do not upload the outer ZIP as the only repository file.
4. Confirm that hidden files such as `.gitignore` and the `.github/workflows/` directory are included.
5. Review the complete file list and commit to `main`.

GitHub's browser may not preserve empty directories, which is why the example input/output directories contain README files.

## Updating an existing `STR-Interruption` repository

If the earlier repository already contains useful history, it can be renamed to `MotifSTaR` in GitHub repository settings instead of creating an unrelated repository. After the rename, update the local remote:

```bash
git remote set-url origin https://github.com/Mamunur-20/MotifSTaR.git
git remote -v
```

Then replace or merge the prepared files carefully and commit the branding change. GitHub usually redirects the old repository URL, but the manuscript, Zenodo record, `CITATION.cff`, and badges should use the new canonical URL.

## Add public example SVGs and output later

The safety `.gitignore` excludes the SVG and generated-output areas. Add only reviewed public files explicitly:

```bash
git add -f examples/input/svgs/sample_001/sample_001.eh_realigned.reviewer.ZFHX3.svg
git add -f examples/output/public_demo/run_manifest.json
git add -f examples/output/public_demo/figure_data
git add -f examples/output/public_demo/figures
git status
git commit -m "Add redistributable example input and output"
git push
```

Do not add a private cohort output directory and then try to clean it later; keep it outside the public repository from the start.

## Archive with Zenodo

After the public GitHub release is final:

1. Link the GitHub repository to Zenodo under the approved account/organization.
2. Enable archiving for `MotifSTaR`.
3. Create or re-publish the GitHub release so Zenodo archives that tag.
4. Add the Zenodo DOI to the manuscript, README, and `CITATION.cff`.
5. Keep the tagged GitHub commit unchanged; later improvements should use a new version tag.
