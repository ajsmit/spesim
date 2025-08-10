# spesim GitHub Workflow & Best Practices

*A friendly, end‑to‑end guide to maintaining `main`, developing on `dev`, using feature branches, and shipping clean releases for the **spesim** R package.*

---

## 0) TL;DR Quick Reference

- **Default branch**: `main` (stable, installable).
- **Staging branch**: `dev` (integration of finished features).
- **Work branches**: `feature/<short-desc>` or `fix/<short-desc>` from `dev`.
- **Flow**: `feature` → PR into `dev` → test on `dev` → PR into `main` → tag and release.
- **Versioning**: Semantic Versioning (MAJOR.MINOR.PATCH), e.g. `0.1.0`.
- **Releases**: Tag (`vX.Y.Z`), GitHub Release notes, NEWS.md entry.
- **CI**: GitHub Actions run R CMD check on PRs to `dev` and `main`.

---

## 1) Branching Model

### Main Branch (`main`)

- Always **stable** and **installable** via `devtools::install_github()`.
- Protected: require PR review + all checks pass.
- Only receives merges from `dev` (release merges) and emergency hotfix PRs.

### Development Branch (`dev`)

- Integration branch for features that are **done** (reviewed & passing checks).
- Can break occasionally, but aim to keep it green.
- Frequently merged forward/back with `main` to reduce drift.

### Feature Branches

- `feature/<short-description>` or `fix/<short-description>`.
- Created from `dev`. Merge back into `dev` via PR when done.

```bash
# Ensure you're up to date locally
git checkout dev
git pull origin dev

# Create and switch to a new feature branch
git checkout -b feature/add-thomas-process
```

---

## 2) Daily Development Cycle

1) **Work locally** on your feature branch.  
2) **Commit often** with clear, atomic messages.  
3) **Run checks** early and often.  
4) **Push** and open a PR into `dev` when ready.

```bash
# Stage and commit
git add -A
git commit -m "Implement Thomas/Neyman–Scott process and tests"

# Run R checks locally
R -q -e "devtools::test(); devtools::check(document = FALSE)"
# or
R -q -e "rcmdcheck::rcmdcheck('--as-cran')"

# Push to remote
git push -u origin feature/add-thomas-process
```

---

## 3) Pull Request (PR) Workflow

- Open PR **from** your feature branch **into** `dev`.
- Fill out the PR template (what/why/how, tests, docs, screenshots if relevant).
- CI must be green (R CMD check, linting, etc.).
- Request review; address comments with follow-up commits.
- **Squash & merge** (default) to keep history clean.

**PR title style**

- `feat: add Neyman–Scott (Thomas) process`
- `fix: repair Voronoi clipping when domain has holes`
- `docs: expand README with usage examples`

**PR checklist**

- [ ] Unit tests added/updated
- [ ] Roxygen docs updated + `devtools::document()` run
- [ ] `NEWS.md` entry added
- [ ] CI green on Linux/macOS/Windows
- [ ] No unchecked `TODO`/`FIXME` in diff

---

## 4) Merging `dev` → `main` (Release)

1) Prepare a release PR from `dev` into `main`:

```bash
git checkout dev
git pull
# optional: sync with latest main
git checkout main && git pull && git checkout dev && git merge --no-ff main
```

2) Bump version in `DESCRIPTION` and update `NEWS.md`:

```bash
# In DESCRIPTION
Version: 0.2.0
```

3) Rebuild docs and run full checks:

```r
devtools::document()
devtools::check()
```

4) Open PR `dev` → `main`. After review & green CI, **merge**.

5) Tag and create a GitHub Release:

```bash
git checkout main
git pull
git tag -a v0.2.0 -m "spesim 0.2.0"
git push origin v0.2.0
```

- Draft release notes summarizing changes (pull from `NEWS.md`).

---

## 5) Hotfixes

For urgent fixes to `main`:

```bash
git checkout main
git pull
git checkout -b hotfix/fix-cran-error
# commit fix
git push -u origin hotfix/fix-cran-error
# PR into main, merge after CI passes
```

Then **back-merge** to `dev` to keep it in sync:

```bash
git checkout dev
git pull
git merge --no-ff main
git push
```

---

## 6) Versioning & Change Log

- **Semantic Versioning**: `MAJOR.MINOR.PATCH`
  - **MAJOR**: breaking API changes.
  - **MINOR**: backward-compatible features.
  - **PATCH**: bug fixes, no API changes.
- Maintain `NEWS.md` with headings per version:

```markdown
# spesim 0.2.0
- feat: Thomas process for clustered species
- fix: graceful failure when no quadrats fit domain
- perf: faster rank-abundance calculation
```

---

## 7) Issue Tracking & Labels

Use GitHub issues with consistent labels:

- `bug`, `enhancement`, `documentation`, `question`, `good first issue`, `help wanted`, `performance`.
- Assign milestones (`v0.2.0`) to group issues by release.
- Use Project boards (Kanban style) for planning if helpful.

---

## 8) CI/CD (GitHub Actions)

Recommended checks on PRs and pushes to `dev`/`main`:

- R CMD check on Linux/macOS/Windows.
- `lintr` (optional).
- Build vignettes and roxygen docs (as needed).

Minimal example workflow (save as `.github/workflows/R-CMD-check.yaml`):

```yaml
name: R-CMD-check
on:
  push:
    branches: [ main, dev ]
  pull_request:
    branches: [ main, dev ]
jobs:
  R-CMD-check:
    runs-on: ${{ matrix.config.os }}
    name: ${{ matrix.config.os }} (${{ matrix.config.r }})
    strategy:
      fail-fast: false
      matrix:
        config:
          - {os: ubuntu-latest, r: 'release'}
          - {os: macos-latest,  r: 'release'}
          - {os: windows-latest, r: 'release'}
    steps:
      - uses: actions/checkout@v4
      - uses: r-lib/actions/setup-r@v2
        with:
          r-version: ${{ matrix.config.r }}
      - uses: r-lib/actions/setup-r-dependencies@v2
        with:
          extra-packages: any::rcmdcheck
          needs: check
      - uses: r-lib/actions/check-r-package@v2
        with:
          args: 'c("--no-manual")'
```

---

## 9) Local Developer Hygiene

- **Keep branches small** and focused: 1 PR ≈ 1 feature/fix.
- **Rebase** local feature branches on `dev` if diverged:

```bash
git checkout feature/add-thomas-process
git fetch origin
git rebase origin/dev
# Resolve conflicts, then:
git push --force-with-lease
```

- **Avoid large binary files**; use `inst/extdata/` sparingly and consider Git LFS if needed.
- **Run `devtools::document()`** whenever roxygen changes.

---

## 10) Working with Tags & Releases

- Tag releases as `vX.Y.Z` from `main` only.
- Annotated tags preferred (`-a` with message).
- Always push tags: `git push --tags` (or the single tag).

---

## 11) Forking Workflow (if you work across clones)

```bash
# Clone your fork
git clone git@github.com:<you>/spesim.git
cd spesim

# Add the upstream (original repo)
git remote add upstream git@github.com:<owner>/spesim.git

# Keep your dev in sync with upstream dev
git checkout dev
git pull upstream dev
git push origin dev
```

---

## 12) Commit Message Conventions (suggested)

Use a short **type** prefix:

- `feat:`, `fix:`, `docs:`, `perf:`, `refactor:`, `test:`, `build:`, `ci:`

Examples:

- `feat: implement Strauss process for non-dominant species`
- `fix: avoid zero-length sfc in quadrat placement`
- `docs: add end-to-end simulation example to README`

---

## 13) Release Process Checklist

1. Merge all features into `dev` and ensure CI is green.  
2. Bump version in `DESCRIPTION`; update `NEWS.md`.  
3. Open PR `dev` → `main`, get review, ensure CI green.  
4. Merge to `main`.  
5. Tag: `git tag -a vX.Y.Z -m "spesim X.Y.Z"` → `git push origin vX.Y.Z`.  
6. Draft a GitHub Release with highlights, link to closed issues/PRs.  
7. (Optional) Submit to CRAN; create `cran-comments.md` and follow checks.

---

## 14) Troubleshooting & Common Commands

```bash
# See current branches and the one you are on
git branch

# Switch branches
git checkout dev

# Stash local changes (e.g., before switching branches)
git stash
git checkout main
git stash pop   # to re-apply

# Undo last commit but keep changes staged
git reset --soft HEAD~1

# Remove a local branch
git branch -d feature/old-branch
# Remove a remote branch
git push origin --delete feature/old-branch

# See remote URLs
git remote -v
```

---

## 15) Housekeeping Files

- **`README.md`**: clear install/use instructions, examples, badges.
- **`NEWS.md`**: human-readable changes per version.
- **`CONTRIBUTING.md`**: PR/issue guidelines, coding style.
- **`CODE_OF_CONDUCT.md`**: community expectations.
- **`LICENSE`**: e.g., MIT; ensure `DESCRIPTION: License:` matches.
- **`.Rbuildignore`**: keep non-package files out of builds.

---

## 16) Glossary

- **PR (Pull Request)**: proposed changes from a branch to another (e.g., `feature` → `dev`).
- **CI (Continuous Integration)**: automated checks verifying the code builds and passes tests.
- **Tag**: a named pointer to a specific commit, used for releases.
- **Semantic Versioning**: a versioning scheme signaling compatibility with `MAJOR.MINOR.PATCH`.

---

Happy shipping! 🚀
