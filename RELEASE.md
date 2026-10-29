# Release checklist

1. **Bump version** in `CITATION.cff` (`version`, `date-released`) and `nextflow.config` (`manifest.version`).
2. **Update docs** — turn `## Unreleased` in `docs/changelog.md` into `## vX.Y.Z — YYYY-MM-DD`, and add a short section to `docs/index.md` (Release notes; also bump the version in *How to cite*).
3. **Commit & tag** on `dev-version`:
   ```bash
   git add <files>
   git commit -m "Release vX.Y.Z"
   git tag vX.Y.Z
   git push origin dev-version
   git push origin vX.Y.Z
   ```
4. **Merge `dev-version` into `main`** — required to trigger the GitHub Pages rebuild:
   ```bash
   git checkout main
   git merge --ff-only dev-version
   git push origin main
   git checkout dev-version
   ```
5. **Create the GitHub release** (triggers Zenodo archival and the `update-latest-tag` workflow):
   ```bash
   gh release create vX.Y.Z --title "vX.Y.Z" --notes "..."
   ```
6. **Verify** with `gh run list --limit 5` — expect two green runs: *Deploy docs* and *Update latest tag*.

## One-time setup

Needed once before the first release made with this checklist. Zenodo archival and a public documentation site require a public repository.

1. **Make the repository public** (*Settings → General → Danger Zone → Change visibility*).
2. **GitHub Pages** — *Settings → Pages → Build and deployment → Source*: **GitHub Actions**. The *Deploy docs* workflow (`.github/workflows/docs.yml`) then publishes `docs/` to <https://crolllab.github.io/genomedenovo_nf/> on every push to `main` that touches `docs/` or `mkdocs.yml`; it can also be started by hand (*Actions → Deploy docs to GitHub Pages → Run workflow*).
3. **Workflow permissions** — *Settings → Actions → General → Workflow permissions*: **Read and write**, so that *Update latest tag* can move the `latest` tag.
4. **Zenodo** — on <https://zenodo.org/account/settings/github/>, sync and switch on `crolllab/genomedenovo_nf`. Metadata comes from `.zenodo.json`. Release archives leave out `docs/`, `.github/` and `mkdocs.yml` (`.gitattributes`).
5. **After the first archived release**, add the DOI: the concept DOI (all versions) to `CITATION.cff` (`doi:` and `identifiers:`, as in genomepanel_nf), and the DOI badge and citation to `docs/index.md` (placeholders are in HTML comments there).

## Preview the docs locally

```bash
pip install -r docs/requirements.txt
mkdocs serve            # http://127.0.0.1:8000
mkdocs build --strict   # the check the Deploy docs workflow runs
```
