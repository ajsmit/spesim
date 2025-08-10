# SPESIM GitHub Workflow & Best Practices

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

This guide documents how to work on the **spesim** repository using a simple and reliable branch model:

- **`main`** is always stable and releasable.
- **`dev`** is the integration branch where features are merged and tested together.
- Short‑lived **feature branches** branch off `dev` and merge back into `dev` via pull requests.
- Releases are cut by merging `dev` → `main` and tagging a version.

---

## Branch Relationships (Diagram)

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

### Mermaid (rendered on GitHub)

```{mermaid}
gitGraph
    commit id: "v0.1.0 Release"
    branch dev
    commit id: "Dev setup"
    branch feature/add-new-model
    commit id: "Implement Neyman–Scott"
    commit id: "Add Strauss process"
    checkout dev
    merge feature/add-new-model
    commit id: "Integrate feature into dev"
    checkout main
    merge dev id: "Release v0.2.0"
```

### ASCII (for environments that can’t render Mermaid)

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

**Interpretation**

- `dev` starts from `main` after each release.
- New work branches off `dev` as `feature/...` (or `fix/...`).
- When a feature is ready, open a PR into `dev`, get reviews, squash-merge.
- When `dev` is ready to ship, open a PR into `main`, tag a release, and draft notes.

---

## Branch Model Summary

- **`main`**: Always green (passing checks). Only release merges land here. Protected branch.
- **`dev`**: Integrates completed features. May be less stable than `main`, but must pass CI before merge.
- **`feature/*`**: Short-lived branches for focused changes, e.g., `feature/thomas-process`.

---

## Day‑to‑Day Workflow

1. **Sync your local repo**
   ```bash
   git switch dev
   git pull --ff-only origin dev
   ```

2. **Create a feature branch**
   ```bash
   git switch -c feature/short-description
   ```

3. **Do the work**
   - Commit small, logical units.
   - Run local checks often:
     ```bash
     R CMD check --as-cran
     # or using devtools
     R -q -e "devtools::check(document = TRUE, cran = FALSE)"
     ```

4. **Push & open a PR into `dev`**
   ```bash
   git push -u origin feature/short-description
   # then open a PR from feature/... -> dev on GitHub
   ```

5. **Address review & merge**
   - Squash merge with a clear message.
   - Delete the feature branch once merged.

---

## Cutting a Release (merge `dev` → `main`)

1. **Stabilize `dev`**
   - Ensure CI passes, NEWS is updated, DESCRIPTION version is set.

2. **Bump version**
   - Example:
     ```r
     # DESCRIPTION
     Version: 0.2.0
     ```

3. **Merge into `main`**
   - Open PR: `dev` → `main` (Title: `Release v0.2.0`).
   - On merge, create a Git tag:
     ```bash
     git checkout main
     git pull --ff-only
     git tag -a v0.2.0 -m "spesim v0.2.0"
     git push origin v0.2.0
     ```

4. **Draft GitHub release**
   - Use the tag `v0.2.0`.
   - Paste highlights from `NEWS.md` / changelog.
   - Attach built artifacts if relevant.

5. **Restart the cycle**
   - Make sure `dev` branches from the new `main`:
     ```bash
     git switch dev
     git merge --ff-only main
     git push origin dev
     ```

---

## Hotfixes on `main`

Use **hotfix branches** only for urgent fixes:
```bash
git switch main
git pull --ff-only
git switch -c hotfix/urgent-cran-issue
# commit fix
git push -u origin hotfix/urgent-cran-issue
# open PR -> main, merge once green
git tag -a v0.2.1 -m "Hotfix: CRAN issue"
git push origin v0.2.1

# forward-port to dev
git switch dev
git pull --ff-only
git merge --no-ff main
git push origin dev
```

---

## Versioning & Changelog

- Follow **SemVer**: `MAJOR.MINOR.PATCH`.
- Maintain a `NEWS.md` with sections:
  - **Added**, **Changed**, **Fixed**, **Deprecated**, **Removed**, **Security**.

Example entry:
```markdown
## v0.2.0 — 2025-08-10
### Added
- Thomas process for clustered dominant species.
- Strauss & Geyer processes for non-dominants with diagnostics.

### Fixed
- Robust sf-binding for heterogeneous point-process outputs.
```

---

## CI & Required Checks (recommended)

- Run `R CMD check` on pushes to `dev` and PRs into `dev`/`main`.
- Optionally enforce:
  - `lintr` / `styler` for code style.
  - `roxygen2` docs checks.
  - Unit tests via `testthat`.

---

## Pull Request Conventions

- One logical change per PR.
- Title: imperative mood (e.g., “Add Thomas process sampler”).
- Body:
  - Context/motivation
  - Summary of changes
  - Testing notes (how you verified)
  - Breaking changes / migration if any
- Ensure “Allow edits by maintainers” is checked.

---

## Common Workflows

**Start new work from fresh `dev`:**
```bash
git switch dev && git pull --ff-only origin dev
git switch -c feature/my-feature
```

**Update your feature branch with latest `dev`:**
```bash
git switch feature/my-feature
git fetch origin
git rebase origin/dev   # or: git merge origin/dev
# resolve conflicts if any
git push --force-with-lease
```

**Abandon a stale feature branch (locally + remote):**
```bash
git branch -D feature/old-idea
git push origin :feature/old-idea
```

---

## Maintenance Tips

- Keep PRs small; review early and often.
- Write tests alongside features.
- Prefer squash merges to keep a tidy history.
- Tag every release; keep `main` deployable at all times.

---

## FAQ

**Q: Where do I branch from?**  
A: Always branch from `dev` for new work.

**Q: When do I merge to `main`?**  
A: Only for releases (PR: `dev` → `main`).

**Q: What about quick fixes after a release?**  
A: Use a hotfix branch off `main`, then forward-port to `dev`.

---

## Glossary

- **main**: stable, releasable branch.
- **dev**: integration branch for upcoming release.
- **feature/**: short-lived topic branches off `dev`.
- **hotfix/**: urgent fix branches off `main`.
- **tag**: immutable pointer to a specific release commit.

---

## Templates

### PR Title
```
Add <feature>: <short reason or impact>
```

### PR Body
```
## Motivation
<why this change is needed>

## Changes
- <bullet of changes>

## Tests
- <how you tested>

## Notes
- Breaking changes? Deprecations?
```

---

Happy collaborating! Keep `main` clean, `dev` green, and ship often 🚀
