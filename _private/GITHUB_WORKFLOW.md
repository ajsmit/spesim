# GitHub Workflow for `spesim` Development

This document describes the branching strategy, development cycle, and best practices for maintaining and extending the `spesim` package.

---

## 1. Branch Structure

### Main Branch (`main`)

- `main` is the root ancestor for everything.
- It’s the “official” stable release branch
  - always stable, installable, and passing checks (`R CMD check` with 0 errors/warnings/notes if possible).
- Protected (require PR reviews before merging)
  - only tested, review-approved changes go here.
- Tags/releases (like v0.1.0) are based off main.

### Development Branch (`dev`)

- `dev` is created from `main` and acts as the integration branch.
- It’s where multiple feature branches get merged before main
  - used for integrating new features before release.
- Can occasionally break, but aim to keep it functional
  - you might occasionally sync `dev` with `main` to keep it updated with the latest release fixes.
- Frequently merged into `main` for new versions.

### Feature Branches

- Feature branches (e.g., `feature/add-new-model`) are created from `dev`
  - named `feature/<short-description>` or `fix/<short-description>`.
- This is where all active development happens.
- Once a feature is done and tested, it’s merged back into `dev` (not `main` directly).
- Merged back into `dev` when complete.


```
main ──●───────●─────────●───────────────●─────────→ (stable releases)
        ^       ↑                         ↑
        |       |                         |
        |       |                         merge dev into main for release
        |       |
        |    dev ──●──────●───────────●───●────────→ (integration branch)
        |            ↑    ^               |
        |            |    |               merge feature branch into dev
        |            |    |
        |   feature/add-new-model ─●──●───┘ (active development)
        |            (work, commits, tests here)
```

---

## 2. Creating a New Feature Branch

```bash
# Ensure you are up to date
git checkout dev
git pull origin dev

# Create and switch to a new branch
git checkout -b feature/add-new-simulation-model
```

---

## 3. Regular Development Cycle

1. **Work locally** in your feature branch.

2. **Commit often**, with clear commit messages:

```bash
git add .
git commit -m "Implement Neyman–Scott simulation model"
```

3. **Push to GitHub** regularly:

```bash
git push -u origin feature/add-new-simulation-model
```

---

## 4. Updating Your Branch with Latest Changes from `dev`

To keep your branch in sync with `dev`:

```bash
git checkout dev
git pull origin dev
git checkout feature/add-new-simulation-model
git merge dev
```

Resolve any conflicts, then:

```bash
git add .
git commit -m "Merge latest dev into feature branch"
git push
```

---

## 5. Creating a Pull Request (PR)

1. Go to your repository on GitHub.
2. Click **Compare & pull request**.
3. Set:
   - **Base branch** = `dev`
   - **Compare branch** = your feature branch
4. Fill in a clear PR title and description.
5. Request a review from a collaborator.
6. Address feedback, push changes, and re-request review.

---

## 6. Merging PRs

- **Feature branch → dev**: Once approved, squash-merge or rebase-merge into `dev`.
- **dev → main**: Only when preparing a release. Ensure all checks pass.

```bash
git checkout main
git pull origin main
git merge dev
git push origin main
```

Then tag the release (see Section 7).

---

## 7. Tagging a Release

```bash
# From main
git tag -a v0.1.0 -m "First CRAN-ready release"
git push origin v0.1.0
```

Then create a GitHub Release from the tag.

---

## 8. Best Practices

- Keep `main` always stable.
- Never commit directly to `main` or `dev` — always via PRs.
- Keep PRs small and focused.
- Update documentation (`README.md`, vignettes) alongside code changes.
- Run `R CMD check` locally before opening a PR.
- Use `.Rbuildignore` to keep private or dev files out of the build.
