"""
Molecule Fragmentation
======================
Breaks molecules into fragments using two strategies:

1. Bond-cut: Cuts rotatable single bonds not in rings, producing pairs of
   fragments with dummy atom attachment points ([*]).

2. BRICS: Uses retrosynthetic BRICS rules to decompose molecules into
   synthetically meaningful building blocks.

Both methods return standardized fragment data that the recombiner module
can directly consume.
"""

from typing import Dict, List, Tuple

from rdkit import Chem
from rdkit.Chem import BRICS, rdmolops


# --- Bond-cut fragmentation ---

def _get_rotatable_bond_matches(mol: Chem.Mol) -> List[Tuple[int, int]]:
    """Find rotatable single bonds not in rings."""
    pattern = Chem.MolFromSmarts(
        "[!$([NH]!@C(=O))&!D1]-&!@[!$([NH]!@C(=O))&!D1]"
    )
    return list(mol.GetSubstructMatches(pattern))


def fragment_by_bonds(mol: Chem.Mol) -> List[Tuple[str, str, int]]:
    """
    Cut every rotatable bond and return the two SMILES fragments produced
    by each cut, along with the bond index.

    Returns:
        List of (frag1_smiles, frag2_smiles, bond_index) tuples.
        Each fragment contains a dummy atom [*] at the cut site.
    """
    bond_matches = _get_rotatable_bond_matches(mol)
    fragments = []
    seen = set()

    for a1, a2 in bond_matches:
        bond = mol.GetBondBetweenAtoms(a1, a2)
        if bond is None:
            continue
        bond_idx = bond.GetIdx()

        try:
            fragmented = rdmolops.FragmentOnBonds(mol, [bond_idx], addDummies=True)
            frag_smiles = Chem.MolToSmiles(fragmented)
            parts = frag_smiles.split(".")

            if len(parts) == 2:
                key = tuple(sorted(parts))
                if key not in seen:
                    seen.add(key)
                    f1 = Chem.MolFromSmiles(parts[0])
                    f2 = Chem.MolFromSmiles(parts[1])
                    if f1 and f2:
                        fragments.append((parts[0], parts[1], bond_idx))
        except Exception:
            continue

    return fragments


# --- BRICS fragmentation ---

def fragment_by_brics(mol: Chem.Mol) -> List[str]:
    """
    Decompose a molecule using BRICS retrosynthetic rules.

    Returns:
        Sorted list of BRICS fragment SMILES (with numbered dummy atoms
        like [16*], [8*], etc. marking attachment types).
    """
    return sorted(BRICS.BRICSDecompose(mol))


# --- Batch fragmentation ---

def fragment_all_by_bonds(mol_data: Dict[str, Chem.Mol]) -> Dict[str, List[Tuple[str, str, int]]]:
    """
    Fragment every molecule in the dict by bond-cutting.

    Args:
        mol_data: {name: Mol} from the parser module.

    Returns:
        {name: [(frag1_smi, frag2_smi, bond_idx), ...]}
    """
    result = {}
    for name, mol in mol_data.items():
        result[name] = fragment_by_bonds(mol)
    return result


def fragment_all_by_brics(mol_data: Dict[str, Chem.Mol]) -> Dict[str, List[str]]:
    """
    Fragment every molecule in the dict using BRICS.

    Args:
        mol_data: {name: Mol} from the parser module.

    Returns:
        {name: [brics_frag_smi, ...]}
    """
    result = {}
    for name, mol in mol_data.items():
        result[name] = fragment_by_brics(mol)
    return result
