# Release checklist

This repository is private: releases have no Zenodo archival (no `.zenodo.json`, no DOI) and no documentation website (no GitHub Pages). Releases are made directly from `main`.

1. **Bump version** in `CITATION.cff` (`version`, `date-released`) and `nextflow.config` (`manifest.version`).
2. **Update the changelog** — add a section to `CHANGELOG.md`.
3. **Commit & tag** on `main`:
   ```bash
   git add <files>
   git commit -m "Release vX.Y.Z"
   git tag vX.Y.Z
   git push origin main
   git push origin vX.Y.Z
   ```
4. **Create the GitHub release** (stays within the private repository):
   ```bash
   gh release create vX.Y.Z --title "vX.Y.Z" --notes "..."
   ```
   Alternatively, on GitHub: *Releases → Draft a new release*, choose the tag and paste the changelog section as notes.
