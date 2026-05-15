"""
Molecule Parsing
================
Reads molecules from various input sources and returns a standardized
dictionary of {name: RDKit Mol} that all downstream modules consume.

Supported formats:
    - Direct SMILES strings via command line
    - .txt / .smi files (one SMILES per line, optional name after whitespace)
    - .csv files (auto-detects smiles/name columns from headers)
"""

import csv
import os
from typing import Dict, List, Optional

from rdkit import Chem


def parse_smiles_list(smiles_list: List[str], start_index: int = 1) -> Dict[str, Chem.Mol]:
    """
    Parse a list of SMILES strings into named Mol objects.
    Names are auto-generated as mol_1, mol_2, etc.
    """
    molecules = {}
    for idx, smi in enumerate(smiles_list):
        smi = smi.strip()
        if not smi:
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol:
            name = f"mol_{start_index + idx}"
            molecules[name] = mol
        else:
            print(f"  [parse] Warning: Could not parse SMILES '{smi}'")
    return molecules


def parse_text_file(filepath: str) -> Dict[str, Chem.Mol]:
    """
    Parse a .txt or .smi file.
    Each line: SMILES [optional_name]
    Lines starting with # are skipped.
    """
    molecules = {}
    counter = 1
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            smi = parts[0]
            name = parts[1] if len(parts) > 1 else f"mol_{counter}"
            mol = Chem.MolFromSmiles(smi)
            if mol:
                molecules[name] = mol
                counter += 1
            else:
                print(f"  [parse] Warning: Could not parse SMILES '{smi}' (line name: {name})")
    return molecules


def _detect_csv_columns(header: List[str]) -> tuple:
    """Find the SMILES and name column indices from a CSV header row."""
    smi_col = 0
    name_col = None

    smiles_keywords = {"smiles", "smi", "molecule", "structure", "canonical_smiles"}
    name_keywords = {"name", "id", "mol_name", "compound", "compound_id", "title"}

    for idx, h in enumerate(header):
        hl = h.lower().strip()
        if hl in smiles_keywords:
            smi_col = idx
        elif hl in name_keywords:
            name_col = idx

    if name_col is None and len(header) > 1 and smi_col == 0:
        name_col = 1

    return smi_col, name_col


def parse_csv_file(filepath: str) -> Dict[str, Chem.Mol]:
    """
    Parse a CSV file with auto-detected SMILES and name columns.
    Falls back to column 0 for SMILES, column 1 for name.
    """
    molecules = {}
    with open(filepath, "r") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        if not header:
            return molecules

        smi_col, name_col = _detect_csv_columns(header)

        for row_idx, row in enumerate(reader):
            if len(row) <= smi_col:
                continue
            smi = row[smi_col].strip()
            if name_col is not None and len(row) > name_col:
                name = row[name_col].strip()
            else:
                name = f"mol_{row_idx + 1}"

            mol = Chem.MolFromSmiles(smi)
            if mol:
                molecules[name] = mol
            else:
                print(f"  [parse] Warning: Could not parse SMILES '{smi}' for {name}")

    return molecules


def parse_molecules(
    smiles_list: Optional[List[str]] = None,
    filepath: Optional[str] = None,
) -> Dict[str, Chem.Mol]:
    """
    Unified entry point: parse molecules from any supported source.
    Can combine both a file and a SMILES list in one call.

    Returns dict of {name: Mol} consumed by fragmenter, recombiner, etc.
    """
    molecules = {}

    if filepath:
        ext = os.path.splitext(filepath)[1].lower()
        if ext == ".csv":
            molecules.update(parse_csv_file(filepath))
        else:
            molecules.update(parse_text_file(filepath))

    if smiles_list:
        start = len(molecules) + 1
        molecules.update(parse_smiles_list(smiles_list, start_index=start))

    return molecules
