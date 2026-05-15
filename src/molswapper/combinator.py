"""
Combination Calculator
======================
Pre-computes the theoretical number of possible chimeric molecules
before actually generating them.

Given N input molecules, this module calculates:
    - How many rotatable-bond cuts each molecule has
    - How many fragment pairs exist per molecule pair
    - How many total swap combinations are possible (4 per cut-pair)
    - The grand total across all molecule pairs

This lets you know upfront what to expect, especially with large datasets,
and decide whether to apply limits or filters before running generation.
"""

import itertools
import math
from typing import Dict, List

from rdkit import Chem

from . import fragmenter


def count_cuts_per_molecule(mol_data: Dict[str, Chem.Mol]) -> Dict[str, int]:
    """
    Count how many rotatable-bond cuts each molecule produces.

    Returns:
        {name: number_of_cuts}
    """
    counts = {}
    for name, mol in mol_data.items():
        frags = fragmenter.fragment_by_bonds(mol)
        counts[name] = len(frags)
    return counts


def count_brics_fragments_per_molecule(mol_data: Dict[str, Chem.Mol]) -> Dict[str, int]:
    """
    Count how many BRICS fragments each molecule decomposes into.

    Returns:
        {name: number_of_brics_fragments}
    """
    counts = {}
    for name, mol in mol_data.items():
        frags = fragmenter.fragment_by_brics(mol)
        counts[name] = len(frags)
    return counts


def calculate_bond_combinations(mol_data: Dict[str, Chem.Mol]) -> dict:
    """
    Calculate the full combinatorial space for bond-cut swapping.

    For each molecule pair (A, B):
        combinations = cuts_A * cuts_B * 4
        (because each cut-pair produces 4 swap arrangements: A1-B1, A1-B2, A2-B1, A2-B2)

    Returns a report dict with per-pair and total counts.
    """
    cuts = count_cuts_per_molecule(mol_data)
    mol_names = list(cuts.keys())
    n_mols = len(mol_names)
    n_pairs = math.comb(n_mols, 2) if n_mols >= 2 else 0

    pair_details = []
    total_theoretical = 0

    if n_mols == 1:
        name = mol_names[0]
        n_cuts = cuts[name]
        # Self-recombination: C(n_cuts, 2) pairs * 4 arrangements
        self_pairs = math.comb(n_cuts, 2)
        self_total = self_pairs * 4
        total_theoretical = self_total
        pair_details.append({
            "mol_a": name,
            "mol_b": name,
            "cuts_a": n_cuts,
            "cuts_b": n_cuts,
            "cut_pairs": self_pairs,
            "max_products": self_total,
            "note": "self-recombination",
        })
    else:
        for i, j in itertools.combinations(range(n_mols), 2):
            name_a = mol_names[i]
            name_b = mol_names[j]
            ca = cuts[name_a]
            cb = cuts[name_b]
            cut_pairs = ca * cb
            max_products = cut_pairs * 4
            total_theoretical += max_products
            pair_details.append({
                "mol_a": name_a,
                "mol_b": name_b,
                "cuts_a": ca,
                "cuts_b": cb,
                "cut_pairs": cut_pairs,
                "max_products": max_products,
            })

    return {
        "num_molecules": n_mols,
        "num_pairs": n_pairs if n_mols >= 2 else 1,
        "cuts_per_molecule": cuts,
        "pair_details": pair_details,
        "total_theoretical_max": total_theoretical,
        "note": (
            "Theoretical max before deduplication and validity checks. "
            "Actual products will be fewer because some joins produce "
            "invalid chemistry or duplicate SMILES."
        ),
    }


def print_combination_report(report: dict):
    """Print a human-readable combination report to the console."""
    print("\n" + "-" * 60)
    print("COMBINATION SPACE ANALYSIS")
    print("-" * 60)
    print(f"Input molecules: {report['num_molecules']}")
    print(f"Molecule pairs:  {report['num_pairs']}")
    print()

    print("Cuts per molecule:")
    for name, n in report["cuts_per_molecule"].items():
        print(f"  {name}: {n} rotatable-bond cuts")
    print()

    print("Pair breakdown:")
    for p in report["pair_details"]:
        note = f"  ({p['note']})" if "note" in p else ""
        print(f"  {p['mol_a']} x {p['mol_b']}: "
              f"{p['cuts_a']} x {p['cuts_b']} cuts = {p['cut_pairs']} cut-pairs "
              f"-> up to {p['max_products']} products{note}")
    print()

    print(f"THEORETICAL MAXIMUM: {report['total_theoretical_max']} possible molecules")
    print(f"  (actual count will be lower after validity checks and dedup)")
    print("-" * 60)
