# bulknull v0.1.0 Release Checklist

## Pre-Release: Local Setup
```bash
cd ~/dev/bulknull

# Initialize git (if not already done)
git init
git add .
git commit -m "Initial commit: bulknull v0.1.0

- Bayesian dilution_score() function
- Frequentist diagnosability_index() power metric
- Applicability gate for bulk-null filtering
- 52 automated test assertions
- DPP9/M8 case study (k=5 tissue-only: 5 lung IPF series)
- Synthetic validation across cluster fractions 0.01-0.50
"

git branch -M main
```

## Release Steps (PERFORM IN ORDER)

### 1. Push to GitHub (requires credentials)
```bash
git remote add origin https://github.com/YOUR_USERNAME/bulknull.git
git push -u origin main
```

### 2. Enable Zenodo-GitHub Integration
1. Go to https://zenodo.org/account/settings/github
2. Toggle ON for bulknull repository
3. Verify connection status

### 3. Create and Push Release Tag
```bash
git tag -a v0.1.0 -m "Release v0.1.0: bulknull - post-hoc diagnostic for dilution

Framework:
- DDS (Dilution Discordance Score): posterior probability of dilution
- DI (Diagnosability Index): statistical power to detect dilution
- Applicability gate: bulk FDR > 0.05 requirement

Case study: DPP9 in M8 (myeloid) cells from IPF lung tissue (k=5 tissue-only)
- n_cluster = 1332 cells (3.3% of 19,175 total)
- sc_zscore = 2.077409 (strong single-cell signal)
- bulk_beta = -0.084210, SE = 0.062774 (null in bulk)
- DDS = 0.454221 (ambiguous)
- DI = 0.064973 one-sided (very low power)

Validation: Synthetic data (10,000 simulations) confirms DDS ranking and DI power prediction.

Package: 52 test assertions, roxygen2 documentation, vignette, CRAN-compatible."

git push origin v0.1.0
```

### 4. Create GitHub Release
1. Go to https://github.com/YOUR_USERNAME/bulknull/releases
2. Click "Create a new release"
3. Select tag v0.1.0
4. Paste the tag message above as release notes
5. Click "Publish release"

### 5. Zenodo DOI Assignment (automatic, 5-10 minutes)
- Zenodo processes the GitHub release
- Two DOIs generated:
  - **Version DOI**: Specific to v0.1.0 (use for this release)
  - **Concept DOI**: Use this in manuscripts (links all versions)

### 6. Update Package Metadata
After receiving DOIs, update:
- README.md Availability section
- DESCRIPTION (if DOI field exists)
- APPLICATION_NOTE.md Availability section

## Expected DOI Format
- Version DOI: 10.5281/zenodo.XXXXXXX
- Concept DOI: 10.5281/zenodo.YYYYYYY

Record both before manuscript submission.
