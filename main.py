"""
MolSwapper CLI
==============
Entry point that wires together all modules:
    parser -> fragmenter -> recombiner -> properties -> conformer -> output
    combinator runs alongside to report the theoretical combination space.

Usage:
    py -3.13 -m MolSwapper --smiles "c1ccccc1CCO" "c1ccccc1CCNC"
    py -3.13 -m MolSwapper --file mols.csv --method both --lipinski --sdf results.sdf
    py -3.13 -m MolSwapper --file library.smi --deduplicate --sdf output.sdf --mol-dir mol_files/
"""

import argparse
import sys
from datetime import datetime

from rdkit import Chem

from .parser import parse_molecules
from .fragmenter import fragment_all_by_bonds, fragment_all_by_brics
from .recombiner import bond_cut_recombine, brics_recombine, self_recombine
from .properties import enrich_results, deduplicate, filter_by_lipinski, filter_by_range
from .conformer import generate_conformers_batch, generate_input_conformers
from .combinator import calculate_bond_combinations, print_combination_report
from .output import (
    save_csv, save_sdf, save_individual_mol_files,
    print_summary, print_verbose, save_grid_image,
)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Fragment and recombine molecules to generate novel chimeric structures.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic: two molecules, bond-cut
  py -3.13 -m MolSwapper --smiles "c1ccccc1CCO" "c1ccccc1CCNC"

  # Full pipeline with 3D conformers in SDF
  py -3.13 -m MolSwapper --smiles "c1ccccc1CCO" "c1ccccc1CCNC" --sdf output.sdf

  # Large library with filters and individual .mol files
  py -3.13 -m MolSwapper -f library.csv --method both --deduplicate --lipinski \\
      --sdf all_products.sdf --mol-dir mol_files/ --image grid.png

  # Just see the combination space without generating
  py -3.13 -m MolSwapper -f library.smi --count-only

  # Multiple conformers per molecule, UFF force field
  py -3.13 -m MolSwapper --smiles "c1ccccc1CCO" "c1ccccc1CCNC" \\
      --sdf output.sdf --num-conformers 5 --force-field UFF
        """,
    )

    # Input
    inp = p.add_argument_group("Input")
    inp.add_argument("--smiles", nargs="+", help="SMILES strings of input molecules")
    inp.add_argument("--file", "-f", help="Input file (.txt, .csv, .smi)")

    # Method
    meth = p.add_argument_group("Fragmentation")
    meth.add_argument("--method", choices=["bonds", "brics", "both"], default="bonds",
                      help="Fragmentation method (default: bonds)")
    meth.add_argument("--max-products", type=int, default=50,
                      help="Max BRICS products per molecule pair (default: 50)")

    # Filtering
    filt = p.add_argument_group("Filtering")
    filt.add_argument("--deduplicate", action="store_true", help="Remove duplicate products")
    filt.add_argument("--lipinski", action="store_true",
                      help="Keep only Lipinski-compliant products")
    filt.add_argument("--min-weight", type=float, help="Minimum molecular weight")
    filt.add_argument("--max-weight", type=float, help="Maximum molecular weight")
    filt.add_argument("--min-atoms", type=int, help="Minimum heavy atom count")
    filt.add_argument("--max-atoms", type=int, help="Maximum heavy atom count")
    filt.add_argument("--max-logp", type=float, help="Maximum LogP value")

    # 3D conformer generation
    conf = p.add_argument_group("3D Conformer Generation")
    conf.add_argument("--sdf", help="Output SDF file with 3D conformers")
    conf.add_argument("--mol-dir", help="Save each product as individual .mol file in this directory")
    conf.add_argument("--num-conformers", type=int, default=1,
                      help="Number of 3D conformers per molecule (default: 1)")
    conf.add_argument("--force-field", choices=["MMFF", "UFF"], default="MMFF",
                      help="Force field for energy minimization (default: MMFF)")
    conf.add_argument("--no-optimize", action="store_true",
                      help="Skip energy minimization (faster but less accurate)")
    conf.add_argument("--include-inputs", action="store_true",
                      help="Include 3D structures of input molecules in SDF output")

    # Combination analysis
    comb = p.add_argument_group("Combination Analysis")
    comb.add_argument("--count-only", action="store_true",
                      help="Only calculate and print combination space, do not generate")

    # Output
    out = p.add_argument_group("Output")
    out.add_argument("--output", "-o", help="Output CSV path (auto-generated if omitted)")
    out.add_argument("--image", help="Save a molecule grid image (PNG)")
    out.add_argument("--verbose", "-v", action="store_true", help="Print detailed results")

    return p


def main():
    parser = build_parser()
    args = parser.parse_args()

    if not args.smiles and not args.file:
        parser.error("Provide molecules via --smiles or --file")

    # =====================================================
    # Step 1: Parse
    # =====================================================
    mol_data = parse_molecules(smiles_list=args.smiles, filepath=args.file)
    if not mol_data:
        print("Error: No valid molecules parsed.")
        sys.exit(1)
    print(f"Parsed {len(mol_data)} valid molecule(s)")

    # =====================================================
    # Step 2: Combination space analysis (always shown)
    # =====================================================
    combo_report = calculate_bond_combinations(mol_data)
    print_combination_report(combo_report)

    if args.count_only:
        print("\n--count-only mode: stopping before generation.")
        sys.exit(0)

    # =====================================================
    # Step 3: Fragment & Recombine
    # =====================================================
    all_results = []

    if len(mol_data) == 1:
        name, mol = list(mol_data.items())[0]
        print(f"\nSingle molecule mode: self-recombining fragments of '{name}'")

        if args.method in ("bonds", "both"):
            results = self_recombine(mol, name)
            all_results.extend(results)
            print(f"  Bond-cut self-recombination: {len(results)} products")

        if args.method in ("brics", "both"):
            frags = fragment_all_by_brics(mol_data)
            for n, f in frags.items():
                print(f"  BRICS fragments for {n}: {f}")
    else:
        if args.method in ("bonds", "both"):
            results = bond_cut_recombine(mol_data)
            all_results.extend(results)
            print(f"\nBond-cut recombination: {len(results)} products")

        if args.method in ("brics", "both"):
            if args.verbose:
                frags = fragment_all_by_brics(mol_data)
                for n, f in frags.items():
                    print(f"  BRICS fragments for {n}: {f}")
            results = brics_recombine(mol_data, max_per_pair=args.max_products)
            all_results.extend(results)
            print(f"BRICS recombination: {len(results)} products")

    # =====================================================
    # Step 4: Enrich with properties
    # =====================================================
    all_results = enrich_results(all_results)

    # =====================================================
    # Step 5: Filter
    # =====================================================
    if args.deduplicate:
        before = len(all_results)
        all_results = deduplicate(all_results)
        print(f"Deduplicated: {before} -> {len(all_results)} unique products")

    if args.lipinski:
        before = len(all_results)
        all_results = filter_by_lipinski(all_results)
        print(f"Lipinski filter: {before} -> {len(all_results)} compliant products")

    has_range = any([args.min_weight, args.max_weight, args.min_atoms,
                     args.max_atoms, args.max_logp])
    if has_range:
        before = len(all_results)
        all_results = filter_by_range(
            all_results,
            min_weight=args.min_weight,
            max_weight=args.max_weight,
            min_atoms=args.min_atoms,
            max_atoms=args.max_atoms,
            max_logp=args.max_logp,
        )
        print(f"Range filter: {before} -> {len(all_results)} products")

    # =====================================================
    # Step 6: 3D Conformer Generation (if requested)
    # =====================================================
    conformer_data = []
    conformer_count = 0

    if args.sdf or args.mol_dir:
        print(f"\nGenerating 3D conformers ({args.num_conformers} per molecule, "
              f"force field: {args.force_field}, "
              f"optimize: {not args.no_optimize})...")

        conformer_data = generate_conformers_batch(
            all_results,
            num_conformers=args.num_conformers,
            optimize=not args.no_optimize,
            force_field=args.force_field,
        )
        conformer_count = len(conformer_data)

    # =====================================================
    # Step 7: Output
    # =====================================================
    if args.verbose:
        print_verbose(all_results)

    # CSV (always produced)
    output_path = args.output
    if output_path is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_path = f"swapped_molecules_{timestamp}.csv"
    save_csv(all_results, output_path)

    # SDF with 3D structures
    if args.sdf and conformer_data:
        input_mols_3d = None
        if args.include_inputs:
            print("Generating 3D conformers for input molecules...")
            input_mols_3d = generate_input_conformers(
                mol_data,
                num_conformers=args.num_conformers,
                optimize=not args.no_optimize,
            )
        save_sdf(conformer_data, args.sdf, input_mols_3d=input_mols_3d)

    # Individual .mol files
    if args.mol_dir and conformer_data:
        save_individual_mol_files(conformer_data, args.mol_dir)

    # Grid image
    if args.image:
        save_grid_image(all_results, args.image)

    # Final summary
    print_summary(all_results, mol_data, conformer_count=conformer_count)


if __name__ == "__main__":
    main()
