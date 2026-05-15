"""
Fragment Recombination
======================
Takes fragment data produced by the fragmenter module and joins pieces
from different molecules to create novel chimeric structures.

Three modes:
    - bond_cut_recombine:  Swap bond-cut fragments between molecule pairs.
    - brics_recombine:     Use RDKit BRICS build to reassemble BRICS fragments.
    - self_recombine:      Recombine fragments within a single molecule.

All functions return a list of result dicts with standardized keys
that the output module can directly save.
"""

import itertools
from typing import Dict, List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import BRICS, Descriptors

from . import fragmenter


def _join_fragments(frag1_smi: str, frag2_smi: str) -> Optional[str]:
    """
    Connect two fragments at their dummy atoms ([*]).
    Returns canonical SMILES of the product, or None if the join fails.
    """
    frag1 = Chem.MolFromSmiles(frag1_smi)
    frag2 = Chem.MolFromSmiles(frag2_smi)
    if frag1 is None or frag2 is None:
        return None

    dummy = Chem.MolFromSmarts("[#0]")
    f1_matches = frag1.GetSubstructMatches(dummy)
    f2_matches = frag2.GetSubstructMatches(dummy)
    if not f1_matches or not f2_matches:
        return None

    try:
        combo = Chem.CombineMols(frag1, frag2)
        ed = Chem.RWMol(combo)

        d1_idx = f1_matches[0][0]
        d2_idx = f2_matches[0][0] + frag1.GetNumAtoms()

        n1 = ed.GetAtomWithIdx(d1_idx).GetNeighbors()
        n2 = ed.GetAtomWithIdx(d2_idx).GetNeighbors()
        if not n1 or not n2:
            return None

        real_1 = n1[0].GetIdx()
        real_2 = n2[0].GetIdx()

        ed.AddBond(real_1, real_2, Chem.BondType.SINGLE)
        ed.RemoveAtom(max(d1_idx, d2_idx))
        ed.RemoveAtom(min(d1_idx, d2_idx))

        result = ed.GetMol()
        Chem.SanitizeMol(result)
        return Chem.MolToSmiles(result)
    except Exception:
        return None


def _make_product_record(
    smiles: str,
    product_id: str,
    source_a: str,
    source_b: str = "",
    frag_a: str = "",
    frag_b: str = "",
    swap_type: str = "",
) -> Optional[dict]:
    """Build a standardized result dict for one product."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        "product_id": product_id,
        "product_smiles": smiles,
        "source_a": source_a,
        "source_b": source_b,
        "frag_a": frag_a,
        "frag_b": frag_b,
        "swap_type": swap_type,
        "num_atoms": mol.GetNumHeavyAtoms(),
        "mol_weight": round(Descriptors.MolWt(mol), 2),
    }


# --- Bond-cut recombination across molecule pairs ---

def bond_cut_recombine(
    mol_data: Dict[str, Chem.Mol],
) -> List[dict]:
    """
    Fragment each molecule by cutting rotatable bonds, then swap every
    fragment pair between every molecule pair.

    Uses fragmenter.fragment_all_by_bonds() to get fragments, then
    tries all 4 combinations (A1-B1, A1-B2, A2-B1, A2-B2) per cut pair.
    """
    all_frags = fragmenter.fragment_all_by_bonds(mol_data)
    mol_names = list(all_frags.keys())
    results = []

    for i, j in itertools.combinations(range(len(mol_names)), 2):
        name_a = mol_names[i]
        name_b = mol_names[j]
        frags_a = all_frags[name_a]
        frags_b = all_frags[name_b]

        for fa_idx, (fa1, fa2, _) in enumerate(frags_a):
            for fb_idx, (fb1, fb2, _) in enumerate(frags_b):
                combos = [
                    (fa1, fb1, "A1-B1"),
                    (fa1, fb2, "A1-B2"),
                    (fa2, fb1, "A2-B1"),
                    (fa2, fb2, "A2-B2"),
                ]
                for frag_x, frag_y, label in combos:
                    product = _join_fragments(frag_x, frag_y)
                    if product:
                        pid = f"{name_a}_cut{fa_idx+1}_{name_b}_cut{fb_idx+1}_{label}"
                        rec = _make_product_record(
                            smiles=product,
                            product_id=pid,
                            source_a=name_a,
                            source_b=name_b,
                            frag_a=f"{fa1}.{fa2} (cut {fa_idx+1})",
                            frag_b=f"{fb1}.{fb2} (cut {fb_idx+1})",
                            swap_type=label,
                        )
                        if rec:
                            results.append(rec)

    return results


# --- BRICS recombination across molecule pairs ---

def brics_recombine(
    mol_data: Dict[str, Chem.Mol],
    max_per_pair: int = 50,
) -> List[dict]:
    """
    Decompose molecules with BRICS, then use RDKit's BRICSBuild to
    reassemble fragments from different molecules.

    Uses fragmenter.fragment_all_by_brics() for decomposition.
    """
    all_brics = fragmenter.fragment_all_by_brics(mol_data)
    mol_names = list(all_brics.keys())
    results = []

    for i, j in itertools.combinations(range(len(mol_names)), 2):
        name_a = mol_names[i]
        name_b = mol_names[j]

        frag_mols = []
        for s in all_brics[name_a] + all_brics[name_b]:
            m = Chem.MolFromSmiles(s)
            if m:
                frag_mols.append(m)

        if len(frag_mols) < 2:
            continue

        try:
            products = list(itertools.islice(BRICS.BRICSBuild(frag_mols), max_per_pair))
        except Exception:
            products = []

        for idx, prod in enumerate(products):
            try:
                smi = Chem.MolToSmiles(prod)
                rec = _make_product_record(
                    smiles=smi,
                    product_id=f"{name_a}_x_{name_b}_brics_{idx+1}",
                    source_a=name_a,
                    source_b=name_b,
                    swap_type="brics",
                )
                if rec:
                    results.append(rec)
            except Exception:
                continue

    return results


# --- Single-molecule self-recombination ---

def self_recombine(mol: Chem.Mol, name: str) -> List[dict]:
    """
    When only one molecule is provided, fragment it at multiple bonds
    and recombine its own pieces in every arrangement.
    """
    frags = fragmenter.fragment_by_bonds(mol)
    results = []

    for i in range(len(frags)):
        for j in range(i + 1, len(frags)):
            si1, si2, _ = frags[i]
            sj1, sj2, _ = frags[j]

            combos = [
                (si1, sj1, f"cut{i+1}_A-cut{j+1}_A"),
                (si1, sj2, f"cut{i+1}_A-cut{j+1}_B"),
                (si2, sj1, f"cut{i+1}_B-cut{j+1}_A"),
                (si2, sj2, f"cut{i+1}_B-cut{j+1}_B"),
            ]

            for fx, fy, label in combos:
                product = _join_fragments(fx, fy)
                if product:
                    rec = _make_product_record(
                        smiles=product,
                        product_id=f"{name}_self_{label}",
                        source_a=name,
                        frag_a=fx,
                        frag_b=fy,
                        swap_type=label,
                    )
                    if rec:
                        results.append(rec)

    return results
