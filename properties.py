"""
Molecular Properties & Filtering
=================================
Computes descriptors for generated products and provides filtering
to keep only drug-like or user-specified ranges.

Works on the result dicts produced by the recombiner module,
enriching them with additional properties or removing unwanted entries.
"""

from typing import Dict, List, Optional

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


def compute_properties(smiles: str) -> Optional[Dict]:
    """
    Compute a standard set of molecular descriptors for one SMILES.
    Returns None if the SMILES is invalid.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    return {
        "mol_weight": round(Descriptors.MolWt(mol), 2),
        "logp": round(Descriptors.MolLogP(mol), 2),
        "hbd": rdMolDescriptors.CalcNumHBD(mol),
        "hba": rdMolDescriptors.CalcNumHBA(mol),
        "tpsa": round(Descriptors.TPSA(mol), 2),
        "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "num_rings": rdMolDescriptors.CalcNumRings(mol),
        "num_atoms": mol.GetNumHeavyAtoms(),
    }


def enrich_results(results: List[dict]) -> List[dict]:
    """
    Add full property columns to every result dict.
    Operates in-place and also returns the list for chaining.
    """
    for r in results:
        props = compute_properties(r["product_smiles"])
        if props:
            r.update(props)
    return results


def deduplicate(results: List[dict]) -> List[dict]:
    """Remove duplicate products by canonical SMILES."""
    seen = set()
    unique = []
    for r in results:
        smi = r["product_smiles"]
        if smi not in seen:
            seen.add(smi)
            unique.append(r)
    return unique


def filter_by_lipinski(results: List[dict]) -> List[dict]:
    """Keep only products passing Lipinski's Rule of Five."""
    filtered = []
    for r in results:
        mw = r.get("mol_weight", 0)
        logp = r.get("logp", 0)
        hbd = r.get("hbd", 0)
        hba = r.get("hba", 0)
        violations = sum([
            mw > 500,
            logp > 5,
            hbd > 5,
            hba > 10,
        ])
        if violations <= 1:
            r["lipinski_violations"] = violations
            filtered.append(r)
    return filtered


def filter_by_range(
    results: List[dict],
    min_weight: Optional[float] = None,
    max_weight: Optional[float] = None,
    min_atoms: Optional[int] = None,
    max_atoms: Optional[int] = None,
    max_logp: Optional[float] = None,
) -> List[dict]:
    """Filter results by user-specified property ranges."""
    filtered = []
    for r in results:
        if min_weight and r.get("mol_weight", 0) < min_weight:
            continue
        if max_weight and r.get("mol_weight", 0) > max_weight:
            continue
        if min_atoms and r.get("num_atoms", 0) < min_atoms:
            continue
        if max_atoms and r.get("num_atoms", 0) > max_atoms:
            continue
        if max_logp and r.get("logp", 0) > max_logp:
            continue
        filtered.append(r)
    return filtered
