# MolSwapper

A modular toolkit for molecular fragment swapping and chimeric structure generation. MolSwapper takes a set of molecules, fragments each one, and recombines pieces across molecules to generate novel chimeric structures with full 3D conformers and molecular property calculations.

## Overview

MolSwapper automates the computational process of:

1. Parsing molecules from SMILES strings or files (`.txt`, `.csv`, `.smi`)
2. Analyzing the combinatorial space to estimate output before generation
3. Fragmenting molecules by rotatable bonds or BRICS retrosynthetic rules
4. Recombining fragments across all molecule pairs
5. Computing molecular properties (MW, LogP, HBD, HBA, TPSA, etc.)
6. Filtering products by Lipinski's Rule of Five, molecular weight, and other criteria
7. Generating energy-minimized 3D conformers using ETKDG and force fields
8. Exporting results to CSV, SDF (3D structures), individual .mol files, and grid images

For single molecule inputs, MolSwapper operates in self-recombination mode, fragmenting and recombining the molecule's own pieces in different arrangements.

## Core Pipeline Architecture

```mermaid
graph TD
    A["Input<br/>SMILES or Files"]
    B["Parse Molecules"]
    C["Analyze Combination Space"]
    D["Fragment by Bonds or BRICS"]
    E["Recombine Fragments"]
    F["Calculate Properties"]
    G["Apply Filters"]
    H["Generate 3D Conformers"]
    I["Export Results"]
    
    A --> B
    B --> C
    C --> D
    D --> E
    E --> F
    F --> G
    G --> H
    H --> I
    
    style A fill:#E8F4F8,stroke:#1a4d5c,stroke-width:2px,color:#1a4d5c
    style B fill:#E8F4F8,stroke:#1a4d5c,stroke-width:2px,color:#1a4d5c
    style C fill:#FFF4E6,stroke:#4d3319,stroke-width:2px,color:#4d3319
    style D fill:#E8F4F8,stroke:#1a4d5c,stroke-width:2px,color:#1a4d5c
    style E fill:#F0E8F4,stroke:#3d1a4d,stroke-width:2px,color:#3d1a4d
    style F fill:#E8F4F8,stroke:#1a4d5c,stroke-width:2px,color:#1a4d5c
    style G fill:#FFF4E6,stroke:#4d3319,stroke-width:2px,color:#4d3319
    style H fill:#E8F4F8,stroke:#1a4d5c,stroke-width:2px,color:#1a4d5c
    style I fill:#F0E8F4,stroke:#3d1a4d,stroke-width:2px,color:#3d1a4d
```

## Module Architecture

```mermaid
graph TD
    Main["main.py<br/>CLI Orchestrator"]
    
    Parser["parser.py<br/>Load Molecules"]
    Combo["combinator.py<br/>Estimate Space"]
    Frag["fragmenter.py<br/>Fragment Molecules"]
    Recomb["recombiner.py<br/>Swap Fragments"]
    Props["properties.py<br/>Calculate & Filter"]
    Conf["conformer.py<br/>Generate 3D"]
    Output["output.py<br/>Export Results"]
    
    Main --> Parser
    Parser --> Combo
    Combo --> Frag
    Frag --> Recomb
    Recomb --> Props
    Props --> Conf
    Conf --> Output
    
    style Parser fill:#E8F4F8,stroke:#1a4d5c,stroke-width:2px,color:#1a4d5c
    style Combo fill:#FFF4E6,stroke:#4d3319,stroke-width:2px,color:#4d3319
    style Frag fill:#F0E8F4,stroke:#3d1a4d,stroke-width:2px,color:#3d1a4d
    style Recomb fill:#E8F4F8,stroke:#1a4d5c,stroke-width:2px,color:#1a4d5c
    style Props fill:#FFF4E6,stroke:#4d3319,stroke-width:2px,color:#4d3319
    style Conf fill:#F0E8F4,stroke:#3d1a4d,stroke-width:2px,color:#3d1a4d
    style Output fill:#E8F4F8,stroke:#1a4d5c,stroke-width:2px,color:#1a4d5c
    style Main fill:#F0E8F4,stroke:#3d1a4d,stroke-width:2px,color:#3d1a4d
```

## Fragment Swapping Mechanism

```mermaid
graph TD
    MolA["Molecule A<br/>c1ccccc1CCO"]
    MolB["Molecule B<br/>c1ccccc1CCNC"]
    
    MolA --> CutA["Cut at<br/>Rotatable Bonds"]
    MolB --> CutB["Cut at<br/>Rotatable Bonds"]
    
    CutA --> FragA1["Fragment A1<br/>Left Half"]
    CutA --> FragA2["Fragment A2<br/>Right Half"]
    
    CutB --> FragB1["Fragment B1<br/>Left Half"]
    CutB --> FragB2["Fragment B2<br/>Right Half"]
    
    FragA1 --> P1["Product 1<br/>A1 + B1"]
    FragA1 --> P2["Product 2<br/>A1 + B2"]
    FragA2 --> P3["Product 3<br/>A2 + B1"]
    FragA2 --> P4["Product 4<br/>A2 + B2"]
    
    FragB1 --> P1
    FragB2 --> P2
    FragB1 --> P3
    FragB2 --> P4
    
    P1 --> Enrich["Calculate Properties<br/>MW, LogP, HBD, HBA, TPSA"]
    P2 --> Enrich
    P3 --> Enrich
    P4 --> Enrich
    
    Enrich --> Output["Export Results<br/>CSV + SDF + MOL"]
    
    style MolA fill:#E8F4F8,stroke:#1a4d5c,stroke-width:2px,color:#1a4d5c
    style MolB fill:#E8F4F8,stroke:#1a4d5c,stroke-width:2px,color:#1a4d5c
    style CutA fill:#FFF4E6,stroke:#4d3319,stroke-width:2px,color:#4d3319
    style CutB fill:#FFF4E6,stroke:#4d3319,stroke-width:2px,color:#4d3319
    style FragA1 fill:#F0E8F4,stroke:#3d1a4d,stroke-width:2px,color:#3d1a4d
    style FragA2 fill:#F0E8F4,stroke:#3d1a4d,stroke-width:2px,color:#3d1a4d
    style FragB1 fill:#F0E8F4,stroke:#3d1a4d,stroke-width:2px,color:#3d1a4d
    style FragB2 fill:#F0E8F4,stroke:#3d1a4d,stroke-width:2px,color:#3d1a4d
    style P1 fill:#E8F4F8,stroke:#1a4d5c,stroke-width:2px,color:#1a4d5c
    style P2 fill:#E8F4F8,stroke:#1a4d5c,stroke-width:2px,color:#1a4d5c
    style P3 fill:#E8F4F8,stroke:#1a4d5c,stroke-width:2px,color:#1a4d5c
    style P4 fill:#E8F4F8,stroke:#1a4d5c,stroke-width:2px,color:#1a4d5c
    style Enrich fill:#FFF4E6,stroke:#4d3319,stroke-width:2px,color:#4d3319
    style Output fill:#F0E8F4,stroke:#3d1a4d,stroke-width:2px,color:#3d1a4d
```

## Fragmentation Methods

```mermaid
graph TD
    Input["Input Molecule"]
    
    Input --> BC["Bond-Cut Method<br/>Rotatable Bonds"]
    Input --> BR["BRICS Method<br/>Retrosynthetic Rules"]
    
    BC --> BCDesc["Cuts all rotatable bonds<br/>not in rings<br/>Produces fragment pairs<br/>More combinations"]
    
    BR --> BRDesc["Uses synthetic rules<br/>Decomposition<br/>Synthetically valid<br/>Fewer combinations"]
    
    BCDesc --> BCOut["Output:<br/>More fragments<br/>Higher diversity"]
    
    BRDesc --> BROut["Output:<br/>Fewer fragments<br/>Better chemistry"]
    
    style Input fill:#E8F4F8,stroke:#1a4d5c,stroke-width:2px,color:#1a4d5c
    style BC fill:#FFF4E6,stroke:#4d3319,stroke-width:2px,color:#4d3319
    style BR fill:#F0E8F4,stroke:#3d1a4d,stroke-width:2px,color:#3d1a4d
    style BCDesc fill:#FFF4E6,stroke:#4d3319,stroke-width:2px,color:#4d3319
    style BRDesc fill:#F0E8F4,stroke:#3d1a4d,stroke-width:2px,color:#3d1a4d
    style BCOut fill:#FFF4E6,stroke:#4d3319,stroke-width:2px,color:#4d3319
    style BROut fill:#F0E8F4,stroke:#3d1a4d,stroke-width:2px,color:#3d1a4d
```

## Filtering Pipeline

```mermaid
graph TD
    Raw["Raw Products<br/>All Generated"]
    
    Dedup["Deduplication<br/>Remove Duplicates"]
    
    Valid["Validation<br/>Remove Invalid"]
    
    Lipinski["Lipinski Filter<br/>MW < 500, LogP < 5<br/>HBD < 5, HBA < 10"]
    
    MW["Molecular Weight<br/>Min/Max Constraints"]
    
    LogP["LogP Filter<br/>Max Value"]
    
    Atoms["Atom Count<br/>Min/Max Constraints"]
    
    Final["Final Products<br/>Filtered Set"]
    
    Raw --> Dedup
    Dedup --> Valid
    Valid --> Lipinski
    Lipinski --> MW
    MW --> LogP
    LogP --> Atoms
    Atoms --> Final
    
    style Raw fill:#E8F4F8,stroke:#1a4d5c,stroke-width:2px,color:#1a4d5c
    style Dedup fill:#FFF4E6,stroke:#4d3319,stroke-width:2px,color:#4d3319
    style Valid fill:#F0E8F4,stroke:#3d1a4d,stroke-width:2px,color:#3d1a4d
    style Lipinski fill:#E8F4F8,stroke:#1a4d5c,stroke-width:2px,color:#1a4d5c
    style MW fill:#FFF4E6,stroke:#4d3319,stroke-width:2px,color:#4d3319
    style LogP fill:#F0E8F4,stroke:#3d1a4d,stroke-width:2px,color:#3d1a4d
    style Atoms fill:#E8F4F8,stroke:#1a4d5c,stroke-width:2px,color:#1a4d5c
    style Final fill:#FFF4E6,stroke:#4d3319,stroke-width:2px,color:#4d3319
```

## Combination Space Analysis

```mermaid
graph TD
    Input["Input: N Molecules<br/>with Rotatable Bonds"]
    
    Count["Count Cuts<br/>Per Molecule"]
    
    Pairs["Identify All<br/>Molecule Pairs"]
    
    Calc["Calculate Cut-Pairs<br/>per Pair<br/>Cut_A × Cut_B"]
    
    Estimate["Estimate Max Products<br/>Cut_Pairs × 4"]
    
    Total["Total Theoretical Maximum<br/>Sum across all pairs"]
    
    Report["Combination Space<br/>Report"]
    
    Input --> Count
    Count --> Pairs
    Pairs --> Calc
    Calc --> Estimate
    Estimate --> Total
    Total --> Report
    
    style Input fill:#E8F4F8,stroke:#1a4d5c,stroke-width:2px,color:#1a4d5c
    style Count fill:#FFF4E6,stroke:#4d3319,stroke-width:2px,color:#4d3319
    style Pairs fill:#F0E8F4,stroke:#3d1a4d,stroke-width:2px,color:#3d1a4d
    style Calc fill:#E8F4F8,stroke:#1a4d5c,stroke-width:2px,color:#1a4d5c
    style Estimate fill:#FFF4E6,stroke:#4d3319,stroke-width:2px,color:#4d3319
    style Total fill:#F0E8F4,stroke:#3d1a4d,stroke-width:2px,color:#3d1a4d
    style Report fill:#E8F4F8,stroke:#1a4d5c,stroke-width:2px,color:#1a4d5c
```

## Requirements

- Python 3.9 or higher
- RDKit (for molecular operations)
- Pillow (for image generation)

## Installation

```bash
git clone https://github.com/yourusername/MolSwapper.git
cd MolSwapper
pip install -r requirements.txt
```

## Quick Start

```bash
# Simple example with two molecules
python -m MolSwapper --smiles "c1ccccc1CCO" "c1ccccc1CCNC" -v

# From a CSV file
python -m MolSwapper -f molecules.csv --method bonds --deduplicate -v

# Full pipeline with all outputs
python -m MolSwapper -f molecules.csv --method both --deduplicate --lipinski \
  --sdf products.sdf --mol-dir mol_files/ --image grid.png -o results.csv -v

# Preview without generating
python -m MolSwapper -f molecules.csv --count-only
```

## Usage

### Input Options

| Flag | Description |
|------|-------------|
| `--smiles "SMI1" "SMI2" ...` | One or more SMILES strings directly on the command line |
| `--file FILE` or `-f FILE` | Path to an input file (`.txt`, `.csv`, or `.smi`) |

**File Format Examples:**

`.txt` or `.smi` — one SMILES per line, optional name after whitespace:
```
c1ccccc1CCO phenylethanol
c1ccccc1CCNC phenylethylamine
CC(=O)Oc1ccccc1C(=O)O aspirin
```

`.csv` — auto-detects columns named `smiles`/`smi`/`molecule` and `name`/`id`/`compound`:
```csv
smiles,name
c1ccccc1CCO,phenylethanol
c1ccccc1CCNC,phenylethylamine
```

### Fragmentation Methods

| Flag | Method | Description |
|------|--------|-------------|
| `--method bonds` | Bond-cut (default) | Cuts every rotatable single bond not in a ring. Produces fragment pairs with dummy atom attachment points. |
| `--method brics` | BRICS | Uses retrosynthetic BRICS rules to decompose into synthetically meaningful building blocks. |
| `--method both` | Both | Runs bond-cut and BRICS, combines all results. |
| `--max-products N` | | Maximum BRICS products per molecule pair (default: 50). |

### Filtering Options

| Flag | Description |
|------|-------------|
| `--deduplicate` | Remove duplicate products by canonical SMILES |
| `--lipinski` | Keep only products passing Lipinski's Rule of Five |
| `--min-weight N` | Minimum molecular weight |
| `--max-weight N` | Maximum molecular weight |
| `--min-atoms N` | Minimum heavy atom count |
| `--max-atoms N` | Maximum heavy atom count |
| `--max-logp N` | Maximum LogP value |

### 3D Conformer Options

| Flag | Description |
|------|-------------|
| `--sdf PATH.sdf` | Generate 3D conformers and save all products to an SDF file with full metadata |
| `--mol-dir DIRECTORY` | Save each product as an individual `.mol` file |
| `--num-conformers N` | Number of 3D conformers to generate per molecule (default: 1) |
| `--force-field MMFF\|UFF` | Force field for energy minimization (default: MMFF) |
| `--no-optimize` | Skip energy minimization |
| `--include-inputs` | Include original input molecules in the SDF file |

### Combination Analysis

| Flag | Description |
|------|-------------|
| `--count-only` | Only calculate and print the theoretical combination space, then stop |

### Output Options

| Flag | Description |
|------|-------------|
| `--output FILE` or `-o FILE` | Output CSV file path. If omitted, auto-generates filename with timestamp |
| `--image FILE.png` | Save a grid image of the first 24 product molecules |
| `--verbose` or `-v` | Print detailed results table to the console |

## Examples

**Two molecules with basic settings:**
```bash
python -m MolSwapper --smiles "c1ccccc1CCO" "c1ccccc1CCNC"
```

**Full pipeline with 3D conformers in SDF:**
```bash
python -m MolSwapper --smiles "c1ccccc1CCO" "c1ccccc1CCNC" \
  --deduplicate --sdf products.sdf --include-inputs -v
```

**Large library with all output formats:**
```bash
python -m MolSwapper -f library.csv --method both --deduplicate --lipinski \
  --sdf all_products.sdf --mol-dir mol_files/ --image grid.png -o results.csv
```

**Preview combination space before generating:**
```bash
python -m MolSwapper -f big_library.smi --count-only
```

**Multiple conformers with UFF force field:**
```bash
python -m MolSwapper --smiles "c1ccccc1CCO" "c1ccccc1CCNC" \
  --sdf output.sdf --num-conformers 5 --force-field UFF
```

**Single molecule self-recombination:**
```bash
python -m MolSwapper --smiles "CC(=O)Oc1ccccc1C(=O)O" --image products.png -v
```

**Filter by molecular weight and LogP:**
```bash
python -m MolSwapper -f molecules.txt --min-weight 150 --max-weight 400 --max-logp 3.0
```

## Output Formats

### CSV Columns

| Column | Description |
|--------|-------------|
| `product_id` | Unique identifier tracing which fragments were combined |
| `product_smiles` | Canonical SMILES of the generated molecule |
| `source_a` | Name of the first parent molecule |
| `source_b` | Name of the second parent molecule |
| `frag_a` | Fragment from molecule A |
| `frag_b` | Fragment from molecule B |
| `swap_type` | Which fragment halves were joined |
| `num_atoms` | Heavy atom count |
| `mol_weight` | Molecular weight |
| `logp` | Calculated LogP |
| `hbd` | Hydrogen bond donors |
| `hba` | Hydrogen bond acceptors |
| `tpsa` | Topological polar surface area |
| `rotatable_bonds` | Number of rotatable bonds |
| `num_rings` | Ring count |
| `lipinski_violations` | Number of Lipinski violations |

### SDF File Contents

Each molecule entry in the SDF file contains:
- 3D atomic coordinates (energy-minimized)
- All CSV properties embedded as SD tags
- Fragment origin tracking (parent molecule and cut information)
- Molecule name set to product_id for easy lookup

When `--include-inputs` is used, original input molecules appear at the top with `Role: input_molecule` for comparison.

### Individual .mol Files

With `--mol-dir`, each product gets its own `.mol` file named after its product_id. These are useful for loading into visualization tools such as PyMOL, UCSF Chimera, or Avogadro.

## Project Structure

```
MolSwapper/
├── __init__.py           # Package metadata
├── __main__.py           # Enables running as: python -m MolSwapper
├── main.py               # CLI entry point and orchestration
├── parser.py             # Input parsing from various formats
├── fragmenter.py         # Fragmentation engines (bond-cut and BRICS)
├── recombiner.py         # Fragment swapping and joining
├── properties.py         # Molecular descriptors and filtering
├── conformer.py          # 3D conformer generation
├── combinator.py         # Combination space analysis
├── output.py             # Export functions (CSV, SDF, MOL, images)
├── requirements.txt      # Python dependencies
├── setup.py              # Package configuration
├── README.md             # This file
├── LICENSE               # MIT License
└── CONTRIBUTING.md       # Contribution guidelines
```

## Using as a Library

Each module can be imported for custom workflows:

```python
from MolSwapper.parser import parse_molecules
from MolSwapper.fragmenter import fragment_by_bonds, fragment_by_brics
from MolSwapper.recombiner import bond_cut_recombine, brics_recombine
from MolSwapper.properties import enrich_results, deduplicate, filter_by_lipinski
from MolSwapper.conformer import generate_conformers_batch
from MolSwapper.combinator import calculate_bond_combinations
from MolSwapper.output import save_csv, save_sdf, save_individual_mol_files

# Parse molecules
mols = parse_molecules(smiles_list=["c1ccccc1CCO", "c1ccccc1CCNC"])

# Check combination space
report = calculate_bond_combinations(mols)

# Fragment and recombine
results = bond_cut_recombine(mols)

# Enrich with properties and filter
results = enrich_results(results)
results = deduplicate(results)
results = filter_by_lipinski(results)

# Generate 3D conformers
conformer_data = generate_conformers_batch(results, num_conformers=1, optimize=True)

# Save outputs
save_csv(results, "output.csv")
save_sdf(conformer_data, "output.sdf")
```

## Combination Space Report

The combination analysis runs automatically and shows:

```
------------------------------------------------------------
COMBINATION SPACE ANALYSIS
------------------------------------------------------------
Input molecules: 3
Molecule pairs:  3

Cuts per molecule:
  phenylethanol: 2 rotatable-bond cuts
  phenylethylamine: 3 rotatable-bond cuts
  aspirin: 3 rotatable-bond cuts

Pair breakdown:
  phenylethanol x phenylethylamine: 2 x 3 cuts = 6 cut-pairs -> up to 24 products
  phenylethanol x aspirin: 2 x 3 cuts = 6 cut-pairs -> up to 24 products
  phenylethylamine x aspirin: 3 x 3 cuts = 9 cut-pairs -> up to 36 products

THEORETICAL MAXIMUM: 84 possible molecules
(actual count will be lower after validity checks and deduplication)
------------------------------------------------------------
```

Use `--count-only` to see this report without running full generation. This is useful for large datasets to estimate output size.

## Fragment Swapping Formula

For N molecules with c_i rotatable-bond cuts each:

```
Total theoretical max = sum over all pairs (i,j): c_i * c_j * 4
Number of pairs = N * (N-1) / 2
```

Each cut-pair produces 4 possible products:
- A1 + B1 (left + left)
- A1 + B2 (left + right)
- A2 + B1 (right + left)
- A2 + B2 (right + right)

## License

MIT

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on contributing to this project.

## Development

For development setup and guidelines, see [DEVELOPMENT.md](DEVELOPMENT.md).

For a quick reference of all commands, see [QUICK_REFERENCE.md](QUICK_REFERENCE.md).
