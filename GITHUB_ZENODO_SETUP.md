# GitHub and Zenodo Setup for bulknull v0.1.0

## Step 1: Initialize and Push to GitHub

From the project root (`~/dev/bulknull`):

```bash
git init
git add .
git commit -m "Initial commit: bulknull v0.1.0 - post-hoc diagnostic for dilution in cell-type-specific expression"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/bulknull.git
git push -u origin main
```

Replace `YOUR_USERNAME` with your GitHub username.

## Step 2: Configure GitHub Repository (via Web)

Visit https://github.com/YOUR_USERNAME/bulknull/settings and:

- Add description: "Post-hoc diagnostic for whether a bulk-level null result reflects true dilution of a cell-type-specific signal"
- Add topics: `bulk-rnaseq`, `single-cell`, `statistical-power`, `bioinformatics`, `dilution-analysis`
- Enable "Include in the public source code index"

## Step 3: Enable Zenodo-GitHub Integration

1. Go to https://zenodo.org and sign in (or create account)
2. Go to https://zenodo.org/account/settings/github
3. Click the ON/OFF toggle for the `bulknull` repository to enable integration
4. Zenodo will now automatically mint a DOI for each release

## Step 4: Create and Push Release Tag

From the project root:

```bash
git tag -a v0.1.0 -m "Release v0.1.0: bulknull - post-hoc diagnostic for dilution in cell-type-specific expression

This release includes:
- Bayesian dilution_score() function (Dilution Discordance Score)
- frequentist diagnosability_index() function (statistical power)
- applicability gate for bulk-null filtering
- 52 automated test assertions
- vignette with synthetic validation + DPP9/M8 case study"
git push origin v0.1.0
```

## Step 5: Create GitHub Release (via Web)

1. Go to https://github.com/YOUR_USERNAME/bulknull/releases
2. Click "Create a new release"
3. Select tag `v0.1.0`
4. Title: "bulknull v0.1.0: Post-hoc Diagnostic for Dilution"
5. Description: Copy from tag message above
6. Click "Publish release"

Zenodo will automatically create a DOI within minutes.

## Step 6: Record DOIs

After Zenodo processes the release (usually within 5 minutes):

1. Visit https://zenodo.org/search?q=bulknull (or search by username)
2. Find the v0.1.0 release
3. Record two DOIs:
   - **Version DOI** (e.g., 10.5281/zenodo.XXXXXXX): For citing this specific release
   - **Concept DOI** (e.g., 10.5281/zenodo.YYYYYYY): For citing all versions; use this in manuscripts

Update `README.md`, `APPLICATION_NOTE.md`, and `DESCRIPTION` with the Concept DOI:

```
Zenodo DOI: 10.5281/zenodo.YYYYYYY
GitHub: https://github.com/dritikaarora/bulknull
```

## Troubleshooting

- **Zenodo integration not triggering:** Verify the repository is public and the Zenodo integration is enabled in settings
- **DOI not appearing:** Zenodo processes releases asynchronously; check again after 10 minutes
- **Tag already exists:** Delete locally and remotely first:
  ```bash
  git tag -d v0.1.0
  git push origin :refs/tags/v0.1.0
  ```
