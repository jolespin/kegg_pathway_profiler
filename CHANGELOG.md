# Daily Log: 
* [2025.12.18] - Added step-level binary coverage tracking and `step_coverage.tsv.gz` output with MultiIndex columns (pathway, step)
* [2025.12.5] - Added `parse_attribute_from_gff`
* [2025.12.4] - Added parallelization to `profile-pathway-coverage.py`
* [2025.12.4] - Changed the way `__init__.py` loads in submodules
* [2025.12.4] - Added optional `progressbar` to `pathway_coverage_wrapper`
* [2025.10.21] - Add `get-kos-from-pykofamsearch.py`
* [2024.10.16] - Add `pyexeggutor` dependency and replace utils
* [2024.9.21] - Fixed `--download` in `build-pathway-databse.py`
* [2024.8.23] - Initial pre-release

# Pending:
- Use entrypoints instead of scripts in `bin/`