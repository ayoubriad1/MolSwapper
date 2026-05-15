# MolSwapper - Quick Reference

Fast reference guide for commands, options, and common workflows.

## Installation

```bash
pip install -r requirements.txt
```

## Basic Commands

```bash
# Two-molecule example
python -m MolSwapper --smiles "c1ccccc1CCO" "c1ccccc1CCNC"

# From CSV file with filters
python -m MolSwapper -f molecules.csv --method both --deduplicate --lipinski -v

# Preview before generating
python -m MolSwapper -f molecules.csv --count-only

# Complete pipeline with all outputs
python -m MolSwapper -f molecules.csv --method both --deduplicate \
  --sdf output.sdf --mol-dir mol_files/ --image grid.png -o results.csv -v
```

## Input Options

| Command | Description |
|---------|-------------|
| `--smiles "SMI1" "SMI2"` | Direct SMILES strings on command line |
| `-f FILE` or `--file FILE` | Input file (.txt, .csv, .smi) |
| Both together | Combine SMILES strings and files |

### Input File Formats

SMILES file (.txt or .smi):
```
c1ccccc1CCO phenylethanol
c1ccccc1CCNC phenylethylamine
CC(=O)Oc1ccccc1C(=O)O aspirin
```

CSV file:
```csv
smiles,name
c1ccccc1CCO,phenylethanol
c1ccccc1CCNC,phenylethylamine
```

## Fragmentation Options

| Flag | Effect |
|------|--------|
| `--method bonds` | Rotatable bond cutting (default) |
| `--method brics` | BRICS retrosynthetic rules |
| `--method both` | Both methods combined |
| `--max-products N` | Max BRICS products per pair (default: 50) |

## Filtering Options

| Flag | Effect |
|------|--------|
| `--deduplicate` | Remove duplicate SMILES |
| `--lipinski` | Apply Lipinski's Rule of Five |
| `--min-weight N` | Minimum molecular weight |
| `--max-weight N` | Maximum molecular weight |
| `--min-atoms N` | Minimum heavy atom count |
| `--max-atoms N` | Maximum heavy atom count |
| `--max-logp N` | Maximum LogP value |

Lipinski's Rule of Five (when --lipinski is used):
- MW <= 500
- LogP <= 5
- HBD <= 5
- HBA <= 10
- Allows 1 violation

## 3D Conformer Options

| Flag | Effect |
|------|--------|
| `--sdf FILE.sdf` | Generate 3D structures in SDF format |
| `--mol-dir DIR` | Save individual .mol files to directory |
| `--num-conformers N` | Number of conformers per molecule (default: 1) |
| `--force-field MMFF\|UFF` | Energy minimization method (default: MMFF) |
| `--no-optimize` | Skip energy minimization |
| `--include-inputs` | Include original molecules in SDF output |

## Analysis and Output Options

| Flag | Effect |
|------|--------|
| `--count-only` | Calculate theoretical max only, do not generate |
| `-o FILE` or `--output FILE` | CSV output path (auto-generated if omitted) |
| `--image FILE.png` | Save grid image of first 24 products |
| `-v` or `--verbose` | Print detailed results table to console |

## Common Workflows

### Workflow 1: Quick Test
```bash
python -m MolSwapper --smiles "c1ccccc1CCO" "c1ccccc1CCNC" -v
```
Output: CSV results in console and file

### Workflow 2: Drug-Like Filtering
```bash
python -m MolSwapper -f molecules.csv --lipinski --deduplicate -v
```
Output: Only Lipinski-compliant products

### Workflow 3: Full Export
```bash
python -m MolSwapper -f molecules.csv --method both --deduplicate \
  --sdf products.sdf --mol-dir mol_files/ --image grid.png \
  --include-inputs -o results.csv -v
```
Output: CSV, SDF with 3D, individual MOL files, grid image

### Workflow 4: Large Library Preview
```bash
python -m MolSwapper -f large_library.csv --count-only
```
Output: Estimated total products without generation

### Workflow 5: Weight and LogP Constraints
```bash
python -m MolSwapper -f molecules.csv \
  --min-weight 150 --max-weight 400 --max-logp 3.0 -v
```
Output: Products within specified ranges

### Workflow 6: Single Molecule Self-Recombination
```bash
python -m MolSwapper --smiles "CC(=O)Oc1ccccc1C(=O)O" \
  --sdf output.sdf --image products.png -v
```
Output: Self-generated chimeras with 3D visualization

### Workflow 7: Multiple Conformers
```bash
python -m MolSwapper --smiles "c1ccccc1CCO" "c1ccccc1CCNC" \
  --sdf output.sdf --num-conformers 5 --force-field UFF
```
Output: 5 conformers per product, UFF minimization

## CSV Output Column Reference

| Column | Description |
|--------|-------------|
| product_id | Unique identifier tracing parent fragments |
| product_smiles | Generated molecule SMILES |
| source_a | First parent molecule name |
| source_b | Second parent molecule name |
| frag_a | Fragment from molecule A |
| frag_b | Fragment from molecule B |
| swap_type | Recombination arrangement (A1-B2, etc.) |
| num_atoms | Heavy atom count |
| mol_weight | Molecular weight in Daltons |
| logp | LogP partition coefficient |
| hbd | Hydrogen bond donor count |
| hba | Hydrogen bond acceptor count |
| tpsa | Topological polar surface area |
| rotatable_bonds | Rotatable bond count |
| num_rings | Ring system count |
| lipinski_violations | Number of rule violations |

## Performance Optimization

| Tip | Benefit |
|-----|---------|
| Use `--count-only` first | Estimate output before generation |
| Start with `--method bonds` | Faster than BRICS method |
| Apply filters early | Reduces 3D generation overhead |
| Use `--no-optimize` | Skip energy minimization for speed |
| Reduce `--num-conformers` | Default 1 is usually sufficient |
| Split large input | Process multiple smaller batches |

## Python API Quick Start

```python
from MolSwapper.parser import parse_molecules
from MolSwapper.fragmenter import fragment_by_bonds
from MolSwapper.recombiner import bond_cut_recombine
from MolSwapper.properties import enrich_results, filter_by_lipinski
from MolSwapper.conformer import generate_conformers_batch
from MolSwapper.output import save_csv, save_sdf

# Load molecules
mols = parse_molecules(smiles_list=["c1ccccc1CCO", "c1ccccc1CCNC"])

# Generate products
results = bond_cut_recombine(mols)

# Enrich and filter
results = enrich_results(results)
results = filter_by_lipinski(results)

# Generate 3D structures
conformer_data = generate_conformers_batch(results, num_conformers=1)

# Save results
save_csv(results, "output.csv")
save_sdf(conformer_data, "output.sdf")
```

## Troubleshooting

| Issue | Solution |
|-------|----------|
| RDKit not found | Run: `pip install rdkit` |
| No products generated | Verify SMILES validity; check input file format |
| Memory error with large files | Use `--count-only` first; split input; skip 3D |
| Slow performance | Use `--method bonds` instead of `both` |
| Module import error | Run from parent directory of MolSwapper/ |

## File Format Reference

| Extension | Description |
|-----------|-------------|
| .smi | SMILES file format |
| .txt | Text file format |
| .csv | Comma-separated values |
| .sdf | 3D structure file format |
| .mol | Individual molecule file format |
| .png | Image file format |

## Default Configuration

| Parameter | Default Value |
|-----------|---------------|
| Fragmentation method | bonds |
| Force field | MMFF |
| Number of conformers | 1 |
| Max BRICS products | 50 |
| Output filename | swapped_molecules_YYYYMMDD_HHMMSS.csv |

## Help and Documentation

```bash
# Show all available options
python -m MolSwapper --help

# View detailed documentation
cat README.md

# Visual guide
cat VISUAL_GUIDE.md
```

---

Reference: [README.md](README.md) | [VISUAL_GUIDE.md](VISUAL_GUIDE.md) | [DEVELOPMENT.md](DEVELOPMENT.md)
