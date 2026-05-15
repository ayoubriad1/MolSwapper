"""
Output & Reporting
==================
Handles saving results in multiple formats:
    - CSV:  Flat table of all products with properties.
    - SDF:  3D molecular structures with embedded metadata (requires conformer module).
    - PNG:  Grid image of 2D molecule structures.

Consumes the result dicts produced by recombiner and enriched by
the properties module, plus 3D Mol objects from the conformer module.
"""

import csv
import os
from typing import Dict, List, Tuple

from rdkit import Chem
from rdkit.Chem import Draw


# --- CSV output ---

def save_csv(results: List[dict], output_path: str):
    """Save results to a CSV file."""
    if not results:
        print("No results to save.")
        return

    fieldnames = list(results[0].keys())
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"Saved {len(results)} products to {output_path}")


# --- SDF output ---

def save_sdf(
    conformer_data: List[Tuple[dict, Chem.Mol]],
    output_path: str,
    input_mols_3d: Dict[str, Chem.Mol] = None,
):
    """
    Write an SDF file containing 3D structures of all products.

    Each molecule entry includes:
        - 3D atomic coordinates
        - All metadata as SD properties (source molecules, fragments,
          swap type, molecular descriptors)

    Args:
        conformer_data:  List of (result_dict, mol_with_3d) from conformer module.
        output_path:     Path for the .sdf file.
        input_mols_3d:   Optional dict of {name: mol_3d} for original input molecules
                         to include at the top of the SDF for reference.
    """
    if not conformer_data and not input_mols_3d:
        print("No 3D molecules to save.")
        return

    writer = Chem.SDWriter(output_path)
    count = 0

    # Write input molecules first (if provided)
    if input_mols_3d:
        for name, mol in input_mols_3d.items():
            if mol is not None and mol.GetNumConformers() > 0:
                writer.write(mol)
                count += 1

    # Write product molecules
    for result_dict, mol in conformer_data:
        if mol is not None and mol.GetNumConformers() > 0:
            writer.write(mol)
            count += 1

    writer.close()
    print(f"Saved {count} 3D molecules to {output_path}")


def save_individual_mol_files(
    conformer_data: List[Tuple[dict, Chem.Mol]],
    output_dir: str,
):
    """
    Save each product as a separate .mol file in the given directory.
    File names are derived from the product_id.
    """
    os.makedirs(output_dir, exist_ok=True)
    count = 0

    for result_dict, mol in conformer_data:
        if mol is None or mol.GetNumConformers() == 0:
            continue
        pid = result_dict.get("product_id", f"product_{count+1}")
        safe_name = pid.replace("/", "_").replace("\\", "_").replace(" ", "_")
        filepath = os.path.join(output_dir, f"{safe_name}.mol")

        writer = Chem.SDWriter(filepath)
        writer.write(mol)
        writer.close()
        count += 1

    print(f"Saved {count} individual .mol files to {output_dir}")


# --- Summary & verbose ---

def print_summary(
    results: List[dict],
    mol_data: Dict[str, Chem.Mol],
    conformer_count: int = 0,
):
    """Print a human-readable summary of inputs and outputs."""
    print("\n" + "=" * 60)
    print("FRAGMENT SWAPPING SUMMARY")
    print("=" * 60)
    print(f"Input molecules: {len(mol_data)}")
    for name, mol in mol_data.items():
        print(f"  {name}: {Chem.MolToSmiles(mol)}")

    print(f"\nTotal products generated: {len(results)}")

    if results:
        weights = [r.get("mol_weight", 0) for r in results]
        atoms = [r.get("num_atoms", 0) for r in results]
        print(f"Molecular weight range: {min(weights):.1f} - {max(weights):.1f}")
        print(f"Heavy atom count range: {min(atoms)} - {max(atoms)}")

        if "logp" in results[0]:
            logps = [r["logp"] for r in results]
            print(f"LogP range: {min(logps):.2f} - {max(logps):.2f}")

        if "lipinski_violations" in results[0]:
            passing = sum(1 for r in results if r.get("lipinski_violations", 99) <= 1)
            print(f"Lipinski-compliant: {passing}/{len(results)}")

    if conformer_count > 0:
        print(f"3D conformers generated: {conformer_count}")

    print("=" * 60)


def print_verbose(results: List[dict], max_rows: int = 20):
    """Print a detailed table of the first N results."""
    for r in results[:max_rows]:
        pid = r.get("product_id", "?")
        smi = r.get("product_smiles", "?")
        mw = r.get("mol_weight", "?")
        print(f"  {pid}: {smi} (MW={mw})")
    if len(results) > max_rows:
        print(f"  ... and {len(results) - max_rows} more")


# --- Grid image ---

def save_grid_image(
    results: List[dict],
    output_path: str,
    mols_per_row: int = 4,
    max_mols: int = 24,
    img_size: tuple = (300, 300),
):
    """
    Save a grid image of the first N product molecules as PNG.
    """
    mols = []
    legends = []
    for r in results[:max_mols]:
        mol = Chem.MolFromSmiles(r["product_smiles"])
        if mol:
            mols.append(mol)
            legends.append(r.get("product_id", ""))

    if not mols:
        print("No valid molecules to draw.")
        return

    img = Draw.MolsToGridImage(
        mols,
        molsPerRow=mols_per_row,
        subImgSize=img_size,
        legends=legends,
    )
    img.save(output_path)
    print(f"Saved molecule grid image to {output_path}")
