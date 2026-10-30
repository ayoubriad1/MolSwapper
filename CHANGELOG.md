# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2026-10-30

### Added
- Flask web console: single-page UI, 3Dmol.js rendering, JSON endpoint `/api/swap`
- Core test suite: 29 tests in nine classes, run in CI across Python 3.10 to 3.13
- `tests/conftest.py` so `pytest tests/` works in a fresh checkout without installing

### Changed
- Dependencies pinned exactly (`rdkit==2025.03.3`, `Pillow==11.2.1`, `Flask==3.1.3`, `pytest==9.1.1`) instead of open ranges

### Fixed
- CI installed requirements but never the package, so every job failed on `import molswapper`
- CI coverage target named `MolSwapper`; the package is `molswapper`
- `tests/test_basic.py.py` had a double extension and aborted pytest collection

---

## [1.0.0] - 2026-05-14

### Added
- Initial release of MolSwapper
- Bond-cut fragmentation method for rotatable bond cutting
- BRICS retrosynthetic fragmentation method
- Fragment recombination engine for generating chimeric molecules
- Molecular property calculations (MW, LogP, HBD, HBA, TPSA, rotatable bonds, ring count)
- Lipinski's Rule of Five filtering
- 3D conformer generation with energy minimization (MMFF and UFF force fields)
- CSV export with comprehensive molecular descriptors
- SDF export with 3D coordinates and metadata
- Individual .mol file export
- Grid image visualization of first 24 products
- Combination space analysis and theoretical maximum calculation
- CLI with extensive options for fragmentation, filtering, and output
- Library import support for custom workflows
- Support for SMILES input, .txt, .csv, and .smi file formats
- Self-recombination mode for single molecule inputs
- Deduplication by canonical SMILES
- Parallel processing for combination space analysis

### Documentation
- Comprehensive README with usage examples
- Project structure and module connection diagrams
- Contributing guidelines
- Installation instructions
- CLI flag reference and examples
- Library usage examples
- Output format documentation

---

## Planning for Future Releases

### [1.2.0] - Planned
- Performance optimizations for large libraries
- Additional filtering options
- Custom fragmentation callbacks
- Batch processing improvements

### [2.0.0] - Planned
- GUI interface
- Database integration
- Advanced visualization tools
- Expanded force field options
- Constraint-based filtering
