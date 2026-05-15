# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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

### [1.1.0] - Planned
- Unit tests and integration tests
- Test coverage reporting
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
