"""
3D Conformer Generation
=======================
Generates 3D conformations for product molecules using RDKit's
ETKDG (Experimental Torsion Knowledge Distance Geometry) method.

Takes the result dicts from the recombiner/properties pipeline and
produces RDKit Mol objects with embedded 3D coordinates, ready for
SDF/MOL file export by the output module.

Workflow:
    results (from recombiner) -> generate_conformers() -> Mol objects with 3D coords
                                                       -> output.save_sdf()
"""

from typing import Dict, List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom


def generate_single_conformer(
    smiles: str,
    num_conformers: int = 1,
    max_attempts: int = 500,
    random_seed: int = 42,
    optimize: bool = True,
    force_field: str = "MMFF",
) -> Optional[Chem.Mol]:
    """
    Generate 3D conformer(s) for a single SMILES string.

    Args:
        smiles:          Canonical SMILES of the molecule.
        num_conformers:  How many conformers to embed (best is kept if 1).
        max_attempts:    Max embedding attempts for ETKDG.
        random_seed:     Seed for reproducibility.
        optimize:        Whether to energy-minimize after embedding.
        force_field:     "MMFF" (Merck) or "UFF" (Universal).

    Returns:
        RDKit Mol with 3D coordinates, or None if generation failed.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    # Add hydrogens for realistic 3D geometry
    mol = Chem.AddHs(mol)

    # ETKDG parameters
    params = rdDistGeom.ETKDGv3()
    params.randomSeed = random_seed
    params.numThreads = 0  # use all available cores
    params.useSmallRingTorsions = True

    # Embed conformer(s)
    conf_ids = rdDistGeom.EmbedMultipleConfs(mol, numConfs=num_conformers, params=params)

    if len(conf_ids) == 0:
        # Fallback: try with random coordinates
        params.useRandomCoords = True
        conf_ids = rdDistGeom.EmbedMultipleConfs(mol, numConfs=num_conformers, params=params)
        if len(conf_ids) == 0:
            return None

    # Energy minimize each conformer
    if optimize:
        energies = []
        for cid in conf_ids:
            try:
                if force_field == "MMFF":
                    props = AllChem.MMFFGetMoleculeProperties(mol)
                    if props is not None:
                        ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid)
                        if ff is not None:
                            ff.Minimize(maxIts=500)
                            energies.append((cid, ff.CalcEnergy()))
                            continue
                # Fallback to UFF
                ff = AllChem.UFFGetMoleculeForceField(mol, confId=cid)
                if ff is not None:
                    ff.Minimize(maxIts=500)
                    energies.append((cid, ff.CalcEnergy()))
                else:
                    energies.append((cid, float("inf")))
            except Exception:
                energies.append((cid, float("inf")))

        # If we only want 1 conformer, keep the lowest energy one
        if num_conformers == 1 and len(energies) > 1:
            best_cid = min(energies, key=lambda x: x[1])[0]
            # Remove all other conformers
            to_remove = [cid for cid, _ in energies if cid != best_cid]
            for cid in sorted(to_remove, reverse=True):
                mol.RemoveConformer(cid)

    return mol


def generate_conformers_batch(
    results: List[dict],
    num_conformers: int = 1,
    optimize: bool = True,
    force_field: str = "MMFF",
) -> List[Tuple[dict, Chem.Mol]]:
    """
    Generate 3D conformers for every product in the results list.

    Args:
        results:         List of result dicts (must have 'product_smiles' key).
        num_conformers:  Conformers per molecule.
        optimize:        Energy-minimize after embedding.
        force_field:     "MMFF" or "UFF".

    Returns:
        List of (result_dict, mol_with_3d) tuples.
        Molecules that failed 3D generation are skipped.
    """
    output = []
    total = len(results)
    success = 0
    failed = 0

    for idx, r in enumerate(results):
        smi = r["product_smiles"]
        mol = generate_single_conformer(
            smi,
            num_conformers=num_conformers,
            optimize=optimize,
            force_field=force_field,
        )
        if mol is not None:
            # Store metadata as mol properties for SDF output
            mol.SetProp("_Name", r.get("product_id", f"product_{idx+1}"))
            mol.SetProp("SMILES", smi)
            mol.SetProp("Source_A", r.get("source_a", ""))
            mol.SetProp("Source_B", r.get("source_b", ""))
            mol.SetProp("Fragment_A", r.get("frag_a", ""))
            mol.SetProp("Fragment_B", r.get("frag_b", ""))
            mol.SetProp("Swap_Type", r.get("swap_type", ""))

            # Numeric properties
            for key in ["mol_weight", "logp", "tpsa", "hbd", "hba",
                         "rotatable_bonds", "num_rings", "num_atoms"]:
                if key in r:
                    mol.SetProp(key, str(r[key]))

            if "lipinski_violations" in r:
                mol.SetProp("lipinski_violations", str(r["lipinski_violations"]))

            output.append((r, mol))
            success += 1
        else:
            failed += 1

        # Progress for large batches
        if total > 20 and (idx + 1) % 50 == 0:
            print(f"  [conformer] Progress: {idx+1}/{total} "
                  f"(success: {success}, failed: {failed})")

    print(f"  [conformer] 3D generation complete: {success} succeeded, {failed} failed "
          f"out of {total}")

    return output


def generate_input_conformers(
    mol_data: Dict[str, Chem.Mol],
    num_conformers: int = 1,
    optimize: bool = True,
) -> Dict[str, Chem.Mol]:
    """
    Generate 3D conformers for the original input molecules too,
    so they can be included in the SDF output for comparison.
    """
    result = {}
    for name, mol in mol_data.items():
        smi = Chem.MolToSmiles(mol)
        mol_3d = generate_single_conformer(smi, num_conformers=num_conformers, optimize=optimize)
        if mol_3d is not None:
            mol_3d.SetProp("_Name", f"INPUT_{name}")
            mol_3d.SetProp("SMILES", smi)
            mol_3d.SetProp("Role", "input_molecule")
            result[name] = mol_3d
    return result
