# Development Guide

Detailed instructions for setting up a development environment and working on MolSwapper.

## Prerequisites

- Python 3.9 or higher
- Git
- pip or conda (for package management)

## Initial Setup

### 1. Clone the Repository

```bash
git clone https://github.com/yourusername/MolSwapper.git
cd MolSwapper
```

### 2. Create a Virtual Environment

Using venv (recommended):
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

Or using conda:
```bash
conda create -n molswapper python=3.11
conda activate molswapper
```

### 3. Install Dependencies

```bash
pip install -r requirements.txt
```

For development (with testing tools):
```bash
pip install -r requirements.txt
pip install pytest pytest-cov black flake8
```

## Running MolSwapper Locally

### From the Repository Root

```bash
python -m MolSwapper --help
python -m MolSwapper --smiles "c1ccccc1CCO" "c1ccccc1CCNC"
```

### As a Module in Your Own Code

```python
from MolSwapper.parser import parse_molecules
from MolSwapper.fragmenter import fragment_by_bonds
from MolSwapper.recombiner import bond_cut_recombine

mols = parse_molecules(smiles_list=["c1ccccc1CCO", "c1ccccc1CCNC"])
results = bond_cut_recombine(mols)
```

## Project Structure

```
MolSwapper/
├── __init__.py              # Package metadata
├── __main__.py              # CLI entry point
├── main.py                  # CLI orchestration
├── parser.py                # Input parsing
├── fragmenter.py            # Fragmentation engines
├── recombiner.py            # Fragment recombination
├── properties.py            # Molecular descriptors
├── conformer.py             # 3D conformer generation
├── combinator.py            # Combination analysis
├── output.py                # Export functions
├── requirements.txt         # Dependencies
├── setup.py                 # Package configuration
├── README.md                # User documentation
├── CONTRIBUTING.md          # Contribution guidelines
├── DEVELOPMENT.md           # This file
├── CHANGELOG.md             # Release history
├── LICENSE                  # MIT License
├── .gitignore               # Git exclusions
├── .gitattributes           # Git attributes
└── .github/
    └── workflows/
        └── tests.yml        # CI/CD pipeline
```

## Module Overview

### parser.py
- `parse_molecules()` — Load molecules from SMILES, files
- Handles .txt, .csv, .smi formats
- Returns dict of {name: RDKit Mol object}

### fragmenter.py
- `fragment_by_bonds()` — Cut at rotatable bonds
- `fragment_by_brics()` — Use BRICS retrosynthetic rules
- Returns fragment dictionary with attachment point info

### recombiner.py
- `bond_cut_recombine()` — Swap bond-cut fragments
- `brics_recombine()` — Swap BRICS fragments
- `self_recombine()` — Single-molecule self-combination
- Returns list of result dicts with SMILES, properties, etc.

### properties.py
- `enrich_results()` — Calculate molecular descriptors
- `deduplicate()` — Remove duplicate SMILES
- `filter_by_lipinski()` — Apply Lipinski's Rule of Five

### conformer.py
- `generate_conformers_batch()` — 3D structure generation
- `generate_input_conformers()` — For input molecules
- Supports MMFF and UFF force fields

### output.py
- `save_csv()` — Export results table
- `save_sdf()` — Export 3D structures
- `save_individual_mol_files()` — Export single molecules
- `save_grid_image()` — Visualize products

### combinator.py
- `calculate_bond_combinations()` — Theoretical analysis
- `print_combination_report()` — Display report

## Testing

### Running Tests

```bash
pytest tests/ -v
```

### With Coverage

```bash
pytest tests/ --cov=MolSwapper --cov-report=html
```

### Adding New Tests

Create test files in `tests/` directory:
```python
# tests/test_parser.py
import pytest
from MolSwapper.parser import parse_molecules

def test_parse_smiles():
    mols = parse_molecules(smiles_list=["c1ccccc1CCO"])
    assert len(mols) == 1
    assert "SMILES_0" in mols
```

## Code Style

### Format Code with Black

```bash
black MolSwapper/
```

### Check Style with Flake8

```bash
flake8 MolSwapper/ --max-line-length=100
```

### Before Committing

```bash
black MolSwapper/
flake8 MolSwapper/
pytest tests/
```

## Common Development Tasks

### Adding a New Feature

1. Create a feature branch: `git checkout -b feature/my-feature`
2. Make changes to relevant modules
3. Add or update tests in `tests/`
4. Update docstrings and comments
5. Update CHANGELOG.md
6. Submit a pull request

### Fixing a Bug

1. Create a bug fix branch: `git checkout -b bugfix/issue-name`
2. Write a test that reproduces the bug
3. Fix the bug
4. Ensure test passes
5. Update CHANGELOG.md with fix note
6. Submit a pull request

### Adding a New CLI Option

1. Modify `main.py` to add the argument
2. Update relevant module functions
3. Update README.md with flag documentation
4. Add example in README examples section
5. Write tests for the new option

### Optimizing Performance

1. Profile the code: `python -m cProfile -s cumtime -m MolSwapper --smiles ... 2>&1 | head -20`
2. Identify bottlenecks
3. Make changes
4. Benchmark before/after
5. Document improvements in commit message

## Debugging

### Verbose Output

```bash
python -m MolSwapper --verbose --smiles "c1ccccc1CCO"
```

### Python Debugger

```python
import pdb; pdb.set_trace()  # Add breakpoint in code
```

### Logging

Add to your test script:
```python
import logging
logging.basicConfig(level=logging.DEBUG)
```

## Documentation

### Docstring Format

```python
def my_function(param1: str, param2: int) -> list:
    """
    Brief description of what the function does.
    
    Longer description explaining the details, if needed.
    
    Args:
        param1: Description of param1
        param2: Description of param2
        
    Returns:
        Description of return value
        
    Example:
        >>> result = my_function("example", 42)
        >>> print(result)
        [...]
    """
    pass
```

## Release Process

1. Update version in `setup.py` and `__init__.py`
2. Update `CHANGELOG.md` with release notes
3. Commit: `git commit -m "Release v1.x.x"`
4. Tag: `git tag -a v1.x.x -m "Release v1.x.x"`
5. Push: `git push origin main --tags`
6. Create GitHub release from the tag

## Troubleshooting

### RDKit Installation Issues

```bash
# If RDKit fails to install, try conda:
conda install -c conda-forge rdkit
```

### ImportError for Pillow

```bash
pip install --upgrade Pillow
```

### Memory Issues with Large Libraries

Use `--count-only` first to estimate output size:
```bash
python -m MolSwapper -f large_library.csv --count-only
```

---

For questions or issues, open a GitHub issue or discussion.
