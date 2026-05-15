# MolSwapper - Visual Guide

Comprehensive reference guide to understanding MolSwapper's architecture and workflow.

## System Overview

```
┌────────────────────────────────────────────────────────────────┐
│                    MOLSWAPPER WORKFLOW                         │
├────────────────────────────────────────────────────────────────┤
│                                                                  │
│  INPUT              PROCESSING            OUTPUT                │
│  ─────              ──────────            ──────                │
│                                                                  │
│  SMILES ─┐                                                      │
│  Files   ├─ Parse ─ Analyze ─ Fragment ─ Recombine             │
│  CSV ────┘              │        │            │                │
│              Combination     Properties       │                │
│              Space           Calculation     │                 │
│                              │               │                 │
│                         Filtering ───────────┘                 │
│                              │                                 │
│                         3D Generation                          │
│                              │                                 │
│                         Export Results                         │
│                         (CSV, SDF, MOL, PNG)                  │
│                                                                  │
└────────────────────────────────────────────────────────────────┘
```

## Data Flow

```
Input: SMILES or File
          ↓
Load Molecules
          ↓
Analyze Combination Space
  - Count rotatable bonds
  - Calculate cut-pairs
  - Estimate theoretical max
          ↓
Fragment Molecules
  - Bond-cut method OR
  - BRICS method OR
  - Both methods
          ↓
Recombine Fragments
  - Attempt 4 combinations per cut-pair
  - Validate chemistry
  - Generate product SMILES
          ↓
Calculate Properties
  - MW, LogP, HBD, HBA, TPSA
  - Atom count, rings, rotatable bonds
  - Lipinski violations
          ↓
Apply Filters
  - Deduplicate
  - Lipinski rules
  - Weight range
  - LogP limits
  - Atom count
          ↓
Generate 3D Conformers
  - ETKDG coordinate generation
  - Energy minimization (MMFF/UFF)
  - Keep best conformer
          ↓
Export Results
  - CSV (tabular data)
  - SDF (3D structures)
  - MOL files (individual molecules)
  - PNG image (grid preview)
          ↓
Output: Final Product Set
```

## Fragment Swapping Step-by-Step

### Step 1: Input Molecules

```
Molecule A                  Molecule B
(Phenylethanol)            (Phenylethylamine)
   C - C - OH                C - C - NH
   |   |                      |   |
 [Ring] | (rotatable)       [Ring] | (rotatable)
        |                         |
```

### Step 2: Fragmentation

```
Molecule A:                 Molecule B:
   
Fragment A1                Fragment B1
   C                          C
   |                          |
 [Ring]                     [Ring]
 (Left)                     (Left)
   
Fragment A2                Fragment B2
C - OH                      C - NH
(Right)                     (Right)
```

### Step 3: Recombination (4 Products)

```
A1 + B1 ──→ C connected to C with [Ring] groups
A1 + B2 ──→ C connected to C-NH with [Ring] group
A2 + B1 ──→ C-OH connected to C with [Ring] group
A2 + B2 ──→ C-OH connected to C-NH
```

### Step 4: Processing Each Product

```
Product
   ↓
Validate Chemistry
   ↓
Calculate Properties
(MW: 234.5, LogP: 2.1, HBD: 2, HBA: 4, TPSA: 89.4)
   ↓
Apply Filters
(Check Lipinski, MW range, LogP)
   ↓
Generate 3D Structure
(ETKDG + energy minimization)
   ↓
Export to CSV, SDF, MOL
```

## Combination Space Calculation

For 3 molecules with rotatable bonds:

```
Molecule A: 2 rotatable bonds
Molecule B: 3 rotatable bonds
Molecule C: 3 rotatable bonds

Pair A × B:
  Cut combinations: 2 × 3 = 6 cut-pairs
  Products: 6 × 4 = 24 maximum products

Pair A × C:
  Cut combinations: 2 × 3 = 6 cut-pairs
  Products: 6 × 4 = 24 maximum products

Pair B × C:
  Cut combinations: 3 × 3 = 9 cut-pairs
  Products: 9 × 4 = 36 maximum products

Total Theoretical Maximum: 84 products
(Actual count will be lower after validation and deduplication)
```

## Module Dependencies

```
┌──────────────────────────┐
│      main.py             │
│   (CLI Orchestrator)     │
└──────────────────────────┘
         ↓
┌──────────────────────────┐
│     parser.py            │
│   (Load Molecules)       │
└──────────────────────────┘
         ↓
┌──────────────────────────┐
│   combinator.py          │
│  (Estimate Space)        │
└──────────────────────────┘
         ↓
┌──────────────────────────┐
│   fragmenter.py          │
│  (Fragment Molecules)    │
└──────────────────────────┘
         ↓
┌──────────────────────────┐
│   recombiner.py          │
│  (Swap Fragments)        │
└──────────────────────────┘
         ↓
┌──────────────────────────┐
│   properties.py          │
│  (Calculate & Filter)    │
└──────────────────────────┘
         ↓
┌──────────────────────────┐
│   conformer.py           │
│  (Generate 3D)           │
└──────────────────────────┘
         ↓
┌──────────────────────────┐
│   output.py              │
│  (Export Results)        │
└──────────────────────────┘
```

## Filtering Pipeline

```
Raw Products (84 estimated)
         ↓
Deduplicate by canonical SMILES
         ↓
Remove chemically invalid structures
         ↓
Apply Lipinski's Rule of Five
  - MW ≤ 500
  - LogP ≤ 5
  - HBD ≤ 5
  - HBA ≤ 10
         ↓
Filter by molecular weight range
  - Minimum MW
  - Maximum MW
         ↓
Filter by LogP threshold
         ↓
Filter by atom count
  - Minimum atoms
  - Maximum atoms
         ↓
Final Products (typically 20-40)
```

## Output File Structure

```
Working Directory/
├── swapped_molecules_YYYYMMDD_HHMMSS.csv
│   └─ Columns: product_id, SMILES, MW, LogP, HBD, HBA, 
│              TPSA, rotatable_bonds, rings, violations
│
├── products.sdf (if --sdf specified)
│   └─ 3D coordinates for all products
│      - XYZ atomic positions
│      - All CSV properties as SD tags
│      - Molecule naming: product_id
│
├── mol_files/ (if --mol-dir specified)
│   ├── mol_1_cut1_mol_2_cut1_A1-B1.mol
│   ├── mol_1_cut1_mol_2_cut2_A1-B2.mol
│   └── ... (one file per product)
│
└── grid.png (if --image specified)
    └─ Grid visualization of first 24 products
```

## Fragmentation Methods

### Bond-Cut Method

```
Algorithm:
1. Identify all rotatable single bonds
2. Exclude bonds in rings
3. Cut at each bond
4. Create fragment pairs with dummy atoms
5. Recombine in 4 arrangements per cut-pair

Characteristics:
- More fragments generated
- Higher combinatorial diversity
- Less chemically constrained
- Faster execution
- Best for: Exploring maximum chemical space
```

### BRICS Method

```
Algorithm:
1. Apply retrosynthetic BRICS rules
2. Decompose into building blocks
3. Use RDKit BRICSBuild for recombination
4. Follow synthetic guidelines

Characteristics:
- Fewer but more meaningful fragments
- Synthetically feasible products
- Chemically validated decomposition
- Slower execution
- Best for: Drug-like molecule generation
```

## Molecular Properties Reference

```
Property                Abbreviation  Range/Unit
────────────────────────────────────────────────
Molecular Weight        MW            150-600 Da (typical)
Partition Coefficient   LogP          -2 to 6 (log units)
Hydrogen Bond Donors    HBD           0-5 (count)
Hydrogen Bond Acceptors HBA           0-10 (count)
Topological Polar       TPSA          0-200 (Angstroms²)
Surface Area

Rotatable Bonds         RotBonds      0-30 (count)
Ring Count             Rings          0-10 (count)
Heavy Atoms            HeavyAtoms     0-100 (count)
Lipinski Violations    ViolationCt    0-4 (count)
```

## Performance Characteristics

Small Dataset (2-3 molecules):
- Parse time: ~0.1s
- Fragment time: ~0.1s
- Recombine time: ~0.5s
- Property calculation: ~0.2s
- 3D generation: ~2-5s
- Total: ~3-7 seconds

Medium Dataset (10-20 molecules):
- Parse time: ~0.1s
- Fragment time: ~0.5s
- Recombine time: ~10-30s
- Property calculation: ~1-2s
- 3D generation: ~30-60s
- Total: ~45-100 seconds

Large Dataset (50+ molecules):
- Use --count-only first to estimate
- Typical output: 1000-10000 products
- Runtime: minutes to hours
- Recommendation: Apply filters early

## Command Workflow Examples

### Minimal (Preview Only)
```
python -m MolSwapper -f molecules.csv --count-only
→ Output: Theoretical maximum estimate
```

### Basic (CSV Results)
```
python -m MolSwapper --smiles "SMILES1" "SMILES2"
→ Output: swapped_molecules_YYYYMMDD_HHMMSS.csv
```

### Standard (With 3D)
```
python -m MolSwapper -f molecules.csv --method bonds \
  --deduplicate --sdf output.sdf -v
→ Output: CSV + SDF with 3D coordinates
```

### Comprehensive (All Outputs)
```
python -m MolSwapper -f molecules.csv --method both \
  --deduplicate --lipinski --sdf products.sdf \
  --mol-dir mol_files/ --image grid.png -o results.csv -v
→ Output: CSV, SDF, individual MOL files, grid image
```

---

For complete usage documentation, see [README.md](README.md).
For quick command reference, see [QUICK_REFERENCE.md](QUICK_REFERENCE.md).
