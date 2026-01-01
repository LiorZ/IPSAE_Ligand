# ipsae_ligand.py
# Script for calculating the ipSAE score for scoring protein-ligand interactions in AlphaFold3 and Boltz1 models
# Adapted from ipsae.py for protein-protein interactions
# https://www.biorxiv.org/content/10.1101/2025.02.10.637595v1

# This script treats:
#   - Chain A as the protein (or optionally specified chain)
#   - All other chains (B, C, D, etc.) as small molecule ligands
# It produces ipSAE scores for each protein-ligand pair
# 
# Note: Chains that contain standard amino acids are treated as proteins/peptides
# Chains with HETATM records (residue_seq_num = '.') are treated as ligands

# MIT license: script can be modified and redistributed for non-commercial and commercial use

# Usage:
#  For AlphaFold3:
#    python ipsae_ligand.py <path_to_af3_pae_json_file> <path_to_af3_cif_file> <pae_cutoff> <dist_cutoff> [protein_chain]
#    python ipsae_ligand.py fold_protein_ligand_full_data_0.json fold_protein_ligand_model_0.cif 10 10
#
#  For Boltz1:
#    python ipsae_ligand.py <path_to_boltz1_pae_npz_file> <path_to_boltz1_cif_file> <pae_cutoff> <dist_cutoff> [protein_chain]
#    python ipsae_ligand.py pae_protein_ligand_model_0.npz protein_ligand_model_0.cif 10 10

import sys
import os
import json
import numpy as np

np.set_printoptions(threshold=np.inf)


# Define the ptm and d0 functions
def ptm_func(x, d0):
    return 1.0 / (1 + (x / d0) ** 2.0)


ptm_func_vec = np.vectorize(ptm_func)


def calc_d0(L, pair_type='ligand'):
    """Calculate d0 value for scoring.
    
    For ligands, we use a smaller minimum d0 value since ligands are typically small.
    """
    L = float(L)
    if L < 27:
        L = 27
    min_value = 1.0
    if pair_type == 'ligand':
        min_value = 1.0  # Use same minimum for ligand interactions
    d0 = 1.24 * (L - 15) ** (1.0 / 3.0) - 1.8
    return max(min_value, d0)


def calc_d0_array(L, pair_type='ligand'):
    """Calculate d0 values for an array of L values."""
    L = np.array(L, dtype=float)
    L = np.maximum(27, L)
    min_value = 1.0
    return np.maximum(min_value, 1.24 * (L - 15) ** (1.0 / 3.0) - 1.8)


def parse_cif_atom_line(line, fielddict, include_ligands=False):
    """Parse atom line from mmCIF file.
    
    Args:
        line: ATOM/HETATM line from mmCIF file
        fielddict: Dictionary mapping field names to column indices
        include_ligands: If True, parse ligand atoms (those with '.' as residue_seq_num)
    
    Returns:
        Dictionary with atom information, or None if skipping this atom
    """
    linelist = line.split()
    atom_num = linelist[fielddict['id']]
    atom_name = linelist[fielddict['label_atom_id']]
    residue_name = linelist[fielddict['label_comp_id']]
    chain_id = linelist[fielddict['label_asym_id']]
    residue_seq_num = linelist[fielddict['label_seq_id']]
    x = linelist[fielddict['Cartn_x']]
    y = linelist[fielddict['Cartn_y']]
    z = linelist[fielddict['Cartn_z']]

    is_ligand = (residue_seq_num == ".")
    
    if is_ligand and not include_ligands:
        return None

    atom_num = int(atom_num)
    if not is_ligand:
        residue_seq_num = int(residue_seq_num)
    x = float(x)
    y = float(y)
    z = float(z)

    return {
        'atom_num': atom_num,
        'atom_name': atom_name,
        'residue_name': residue_name,
        'chain_id': chain_id,
        'residue_seq_num': residue_seq_num if not is_ligand else 1,  # Use 1 for ligand atoms
        'x': x,
        'y': y,
        'z': z,
        'is_ligand': is_ligand
    }


def get_representative_atom_name(residue_name):
    """Get the representative atom name for a residue/ligand.
    
    For proteins: CA (alpha carbon)
    For nucleic acids: C1' 
    For ligands: Try to use a central heavy atom
    """
    # Standard amino acid CA
    aa_set = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
              "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
    
    # Modified amino acids also use CA
    modified_aa = {"TPO", "SEP", "PTR", "MLY", "MSE"}
    
    # Nucleic acids use C1'
    nuc_set = {"DA", "DC", "DT", "DG", "A", "C", "U", "G"}
    
    if residue_name in aa_set or residue_name in modified_aa:
        return "CA"
    elif residue_name in nuc_set:
        return "C1'"
    else:
        # For ligands, we'll include all heavy atoms (non-H)
        return None  # Will handle all atoms


# Ensure correct usage
if len(sys.argv) < 5:
    print("Usage for AF3 and Boltz1 protein-ligand complexes:")
    print("")
    print("AlphaFold3:")
    print("   python ipsae_ligand.py <path_to_pae_json_file> <path_to_mmcif_file> <pae_cutoff> <dist_cutoff> [protein_chain]")
    print("   python ipsae_ligand.py fold_protein_ligand_full_data_0.json fold_protein_ligand_model_0.cif 10 10")
    print("")
    print("Boltz1:")
    print("   python ipsae_ligand.py <path_to_pae_npz_file> <path_to_mmcif_file> <pae_cutoff> <dist_cutoff> [protein_chain]")
    print("   python ipsae_ligand.py pae_protein_ligand_model_0.npz protein_ligand_model_0.cif 10 10")
    print("")
    print("This script calculates ipSAE scores for protein - ligand interactions")
    print("By default, uses chain A as the protein. Optionally specify the protein chain.")
    sys.exit(1)

pae_file_path = sys.argv[1]
pdb_path = sys.argv[2]
pae_cutoff = float(sys.argv[3])
dist_cutoff = float(sys.argv[4])
specified_protein_chain = sys.argv[5] if len(sys.argv) > 5 else None

pae_string = str(int(pae_cutoff))
if pae_cutoff < 10:
    pae_string = "0" + pae_string
dist_string = str(int(dist_cutoff))
if dist_cutoff < 10:
    dist_string = "0" + dist_string

# Determine file type (AF3 or Boltz1)
if ".cif" in pdb_path and pae_file_path.endswith(".json"):
    af3 = True
    boltz1 = False
elif ".cif" in pdb_path and pae_file_path.endswith(".npz"):
    af3 = False
    boltz1 = True
else:
    print("This script requires:")
    print("  - For AF3: mmCIF (.cif) structure file and JSON (.json) PAE file")
    print("  - For Boltz1: mmCIF (.cif) structure file and NPZ (.npz) PAE file")
    sys.exit(1)

pdb_stem = pdb_path.replace(".cif", "")
path_stem = f'{pdb_path.replace(".cif", "")}_ligand_{pae_string}_{dist_string}'

file_path = path_stem + ".txt"
OUT = open(file_path, 'w')

# Parse mmCIF file to get structure information
# We need to read:
# 1. Protein residues (chain A or specified chain) - use CA atoms
# 2. Ligand atoms (all other chains) - use all atoms for true ligands

protein_residues = []  # List of protein residues with CA coordinates
ligand_atoms = {}  # Dict: chain_id -> list of atoms (for true ligands)
other_chain_residues = {}  # Dict: chain_id -> list of residues (for other protein/peptide chains)
atomsitefield_num = 0
atomsitefield_dict = {}
all_atoms_by_chain = {}  # Store all atoms for each chain first

# Standard amino acids for identifying proteins
residue_set = {"ALA", "ARG", "ASN", "ASP", "CYS",
               "GLN", "GLU", "GLY", "HIS", "ILE",
               "LEU", "LYS", "MET", "PHE", "PRO",
               "SER", "THR", "TRP", "TYR", "VAL"}

# Modified amino acids
modified_aa = {"TPO", "SEP", "PTR", "MLY", "MSE"}

# Nucleic acids
nuc_residue_set = {"DA", "DC", "DT", "DG", "A", "C", "U", "G"}

# Track chain information
chain_has_residues = {}  # chain_id -> True if has sequential residues (is protein/peptide)
chain_has_ligand_atoms = {}  # chain_id -> True if has ligand atoms (. residue seq)

print(f"Reading structure from {pdb_path}...")

# First pass: collect all atoms and determine chain types
with open(pdb_path, 'r') as PDB:
    for line in PDB:
        if line.startswith("_atom_site."):
            line = line.strip()
            (atomsite, fieldname) = line.split(".")
            atomsitefield_dict[fieldname] = atomsitefield_num
            atomsitefield_num += 1

        if line.startswith("ATOM") or line.startswith("HETATM"):
            atom = parse_cif_atom_line(line, atomsitefield_dict, include_ligands=True)
            if atom is None:
                continue
            
            chain_id = atom['chain_id']
            if chain_id not in all_atoms_by_chain:
                all_atoms_by_chain[chain_id] = []
            all_atoms_by_chain[chain_id].append(atom)
            
            if atom['is_ligand']:
                chain_has_ligand_atoms[chain_id] = True
            else:
                chain_has_residues[chain_id] = True

# Determine protein chain
protein_chain = specified_protein_chain
if protein_chain is None:
    # Use the first chain that has residues (usually 'A')
    for chain_id in sorted(all_atoms_by_chain.keys()):
        if chain_has_residues.get(chain_id, False):
            protein_chain = chain_id
            break

if protein_chain is None:
    print("Error: No protein chain found in the structure")
    sys.exit(1)

# Second pass: categorize atoms into protein residues and ligands
for chain_id, atoms in all_atoms_by_chain.items():
    for atom in atoms:
        residue_name = atom['residue_name']
        
        if chain_id == protein_chain:
            # This is our protein chain
            if not atom['is_ligand']:
                if atom['atom_name'] == "CA" or "C1'" in atom['atom_name']:
                    protein_residues.append({
                        'atom_num': atom['atom_num'],
                        'coor': np.array([atom['x'], atom['y'], atom['z']]),
                        'res': residue_name,
                        'chainid': chain_id,
                        'resnum': atom['residue_seq_num'],
                        'residue': f"{residue_name:3}   {chain_id:3} {atom['residue_seq_num']:4}"
                    })
        else:
            # This is a ligand or other chain
            if atom['is_ligand']:
                # True ligand (has '.' for residue seq num)
                if chain_id not in ligand_atoms:
                    ligand_atoms[chain_id] = []
                ligand_atoms[chain_id].append({
                    'atom_num': atom['atom_num'],
                    'coor': np.array([atom['x'], atom['y'], atom['z']]),
                    'res': residue_name,
                    'chainid': chain_id,
                    'atom_name': atom['atom_name'],
                    'identifier': f"{residue_name}_{atom['atom_name']}_{chain_id}",
                    'is_true_ligand': True
                })
            elif residue_name in residue_set or residue_name in modified_aa or residue_name in nuc_residue_set:
                # This is another protein/peptide/nucleic acid chain - treat as "ligand" for our purposes
                if chain_id not in other_chain_residues:
                    other_chain_residues[chain_id] = []
                if atom['atom_name'] == "CA" or "C1'" in atom['atom_name']:
                    other_chain_residues[chain_id].append({
                        'atom_num': atom['atom_num'],
                        'coor': np.array([atom['x'], atom['y'], atom['z']]),
                        'res': residue_name,
                        'chainid': chain_id,
                        'resnum': atom['residue_seq_num'],
                        'residue': f"{residue_name:3}   {chain_id:3} {atom['residue_seq_num']:4}",
                        'is_true_ligand': False
                    })

if protein_chain is None:
    print("Error: No protein chain found in the structure")
    sys.exit(1)

print(f"Found protein chain: {protein_chain} with {len(protein_residues)} residues")
print(f"Found true ligand chains: {list(ligand_atoms.keys())}")
for chain, atoms in ligand_atoms.items():
    print(f"  Chain {chain}: {len(atoms)} atoms, residue: {atoms[0]['res']}")
print(f"Found other protein/peptide chains: {list(other_chain_residues.keys())}")
for chain, residues in other_chain_residues.items():
    print(f"  Chain {chain}: {len(residues)} residues")

# Load PAE data (AF3 or Boltz1)
print(f"\nLoading PAE data from {pae_file_path}...")

if af3:
    # Load AlphaFold3 JSON data
    if os.path.exists(pae_file_path):
        with open(pae_file_path, 'r') as file:
            data = json.load(file)
    else:
        print(f"AF3 PAE file does not exist: {pae_file_path}")
        sys.exit()

    # Get token information from AF3
    token_chain_ids = data['token_chain_ids']
    token_res_ids = data['token_res_ids']
    pae_matrix_full = np.array(data['pae'])
    atom_plddts = np.array(data['atom_plddts'])

elif boltz1:
    # Load Boltz1 NPZ data
    # Boltz1 filenames pattern:
    # structure_model_0.cif
    # pae_structure_model_0.npz
    # plddt_structure_model_0.npz
    # confidence_structure_model_0.json
    
    if not os.path.exists(pae_file_path):
        print(f"Boltz1 PAE file does not exist: {pae_file_path}")
        sys.exit()
    
    # Load PAE matrix
    data_pae = np.load(pae_file_path)
    pae_matrix_full = np.array(data_pae['pae'])
    
    # Load pLDDT values
    plddt_file_path = pae_file_path.replace("pae", "plddt")
    if os.path.exists(plddt_file_path):
        data_plddt = np.load(plddt_file_path)
        atom_plddts = np.array(100.0 * data_plddt['plddt'])
    else:
        print(f"Warning: Boltz1 pLDDT file not found: {plddt_file_path}")
        atom_plddts = np.zeros(pae_matrix_full.shape[0])
    
    # For Boltz1, we need to build token_chain_ids and token_res_ids from the structure
    # Each token corresponds to a residue or atom in order
    token_chain_ids = []
    token_res_ids = []
    
    # Build token list from protein residues first (sorted by chain, then by atom_num)
    all_tokens = []
    for res in protein_residues:
        all_tokens.append((res['atom_num'], res['chainid'], res['resnum']))
    
    # Add ligand atoms
    for chain_id, atoms in ligand_atoms.items():
        for atom in atoms:
            all_tokens.append((atom['atom_num'], atom['chainid'], 1))  # Use 1 for ligand res_id
    
    # Add other chain residues
    for chain_id, residues in other_chain_residues.items():
        for res in residues:
            all_tokens.append((res['atom_num'], res['chainid'], res['resnum']))
    
    # Sort by atom number to get correct order
    all_tokens.sort(key=lambda x: x[0])
    
    # Extract chain_ids and res_ids
    for atom_num, chain_id, res_id in all_tokens:
        token_chain_ids.append(chain_id)
        token_res_ids.append(res_id)

print(f"PAE matrix shape: {pae_matrix_full.shape}")
print(f"Number of tokens: {len(token_chain_ids)}")

# Build mapping from tokens to our structure data
# For proteins: token corresponds to residue
# For ligands: token corresponds to atom

# Get protein token indices (chain A residues)
protein_token_indices = []
for i, (chain, res_id) in enumerate(zip(token_chain_ids, token_res_ids)):
    if chain == protein_chain:
        protein_token_indices.append(i)

print(f"Protein tokens: {len(protein_token_indices)} (indices {protein_token_indices[0]}-{protein_token_indices[-1]})")

# Get ligand token indices for each ligand chain
ligand_token_indices = {}
for chain in ligand_atoms.keys():
    ligand_token_indices[chain] = []
    for i, (tok_chain, res_id) in enumerate(zip(token_chain_ids, token_res_ids)):
        if tok_chain == chain:
            ligand_token_indices[chain].append(i)
    if ligand_token_indices[chain]:
        print(f"Ligand chain {chain} tokens: {len(ligand_token_indices[chain])}")

# Get token indices for other protein/peptide chains
other_chain_token_indices = {}
for chain in other_chain_residues.keys():
    other_chain_token_indices[chain] = []
    for i, (tok_chain, res_id) in enumerate(zip(token_chain_ids, token_res_ids)):
        if tok_chain == chain:
            other_chain_token_indices[chain].append(i)
    if other_chain_token_indices[chain]:
        print(f"Other chain {chain} tokens: {len(other_chain_token_indices[chain])}")

# Calculate ipSAE scores for each protein-ligand pair
print("\n" + "="*80)
print("Calculating ipSAE scores for protein-ligand interactions")
print("="*80)

# Header for output
OUT.write("# IPSAE Ligand Scores - Protein-Ligand Interactions\n")
OUT.write(f"# Protein chain: {protein_chain}\n")
OUT.write(f"# PAE cutoff: {pae_cutoff}, Distance cutoff: {dist_cutoff}\n")
OUT.write("#\n")
OUT.write("Protein  Ligand  LigandName  Type        NumProtRes  NumLigAtoms  ipSAE      ipSAE_max  MeanPAE    MinPAE     NumGoodPAE  MeanPlddt  Model\n")

results = []

# Process true ligands (small molecules)
for lig_chain, lig_atoms_list in ligand_atoms.items():
    if lig_chain not in ligand_token_indices or len(ligand_token_indices[lig_chain]) == 0:
        print(f"Warning: No tokens found for ligand chain {lig_chain}, skipping")
        continue
    
    lig_name = lig_atoms_list[0]['res']  # Get ligand name from first atom
    lig_tok_idx = ligand_token_indices[lig_chain]
    
    # Extract PAE submatrix for protein-ligand interaction
    # Rows: protein tokens, Columns: ligand tokens
    pae_submatrix = pae_matrix_full[np.ix_(protein_token_indices, lig_tok_idx)]
    
    # Calculate statistics
    mean_pae = np.mean(pae_submatrix)
    min_pae = np.min(pae_submatrix)
    
    # Count pairs with PAE < cutoff
    num_good_pae = np.sum(pae_submatrix < pae_cutoff)
    
    # Get plddt values for ligand atoms (with bounds checking)
    try:
        valid_indices = [lig_atoms_list[i]['atom_num'] - 1 
                        for i in range(min(len(lig_atoms_list), len(lig_tok_idx)))
                        if 0 <= lig_atoms_list[i]['atom_num'] - 1 < len(atom_plddts)]
        lig_plddt = atom_plddts[valid_indices] if valid_indices else np.array([])
        mean_plddt = np.mean(lig_plddt) if len(lig_plddt) > 0 else 0.0
    except (IndexError, KeyError):
        mean_plddt = 0.0
    
    # Calculate ipSAE score
    num_prot_res = len(protein_token_indices)
    num_lig_atoms = len(lig_tok_idx)
    
    # Calculate number of "good" interactions for d0
    good_lig_atoms = np.sum(np.any(pae_submatrix < pae_cutoff, axis=0))
    good_prot_res = np.sum(np.any(pae_submatrix < pae_cutoff, axis=1))
    
    n0 = good_prot_res + good_lig_atoms
    if n0 < 1:
        n0 = 1
    
    d0 = calc_d0(n0, 'ligand')
    
    # Calculate ipSAE: for each protein residue, calculate PTM score against ligand atoms
    ipsae_by_res = np.zeros(num_prot_res)
    ipsae_max = 0.0
    
    for i in range(num_prot_res):
        pae_row = pae_submatrix[i, :]
        good_mask = pae_row < pae_cutoff
        
        if np.any(good_mask):
            good_paes = pae_row[good_mask]
            n0_res = np.sum(good_mask)
            d0_res = calc_d0(n0_res, 'ligand')
            ptm_scores = ptm_func_vec(good_paes, d0_res)
            ipsae_by_res[i] = np.mean(ptm_scores)
            if ipsae_by_res[i] > ipsae_max:
                ipsae_max = ipsae_by_res[i]
    
    # Calculate ipSAE using all good PAE values
    good_paes_all = pae_submatrix[pae_submatrix < pae_cutoff]
    if len(good_paes_all) > 0:
        ipsae_all = np.mean(ptm_func_vec(good_paes_all, d0))
    else:
        ipsae_all = 0.0
    
    results.append({
        'protein_chain': protein_chain,
        'ligand_chain': lig_chain,
        'ligand_name': lig_name,
        'ligand_type': 'ligand',
        'num_prot_res': num_prot_res,
        'num_lig_atoms': num_lig_atoms,
        'ipsae': ipsae_all,
        'ipsae_max': ipsae_max,
        'mean_pae': mean_pae,
        'min_pae': min_pae,
        'num_good_pae': num_good_pae,
        'mean_plddt': mean_plddt
    })
    
    # Write to output file
    outstring = (f"{protein_chain:8} {lig_chain:7} {lig_name:11} {'ligand':11} "
                 f"{num_prot_res:10d}  {num_lig_atoms:11d}  "
                 f"{ipsae_all:8.6f}   {ipsae_max:8.6f}   "
                 f"{mean_pae:8.4f}   {min_pae:8.4f}   "
                 f"{num_good_pae:10d}  {mean_plddt:9.2f}  {pdb_stem}\n")
    OUT.write(outstring)
    
    print(f"\nProtein {protein_chain} - Ligand {lig_chain} ({lig_name}):")
    print(f"  Type: Small molecule ligand")
    print(f"  Number of protein residues: {num_prot_res}")
    print(f"  Number of ligand atoms: {num_lig_atoms}")
    print(f"  ipSAE: {ipsae_all:.6f}")
    print(f"  ipSAE_max (per residue): {ipsae_max:.6f}")
    print(f"  Mean PAE: {mean_pae:.4f}")
    print(f"  Min PAE: {min_pae:.4f}")
    print(f"  Number of PAE < {pae_cutoff}: {num_good_pae}")
    print(f"  Mean ligand pLDDT: {mean_plddt:.2f}")

# Process other protein/peptide chains (treat as ligands)
for other_chain, other_residues in other_chain_residues.items():
    if other_chain not in other_chain_token_indices or len(other_chain_token_indices[other_chain]) == 0:
        print(f"Warning: No tokens found for chain {other_chain}, skipping")
        continue
    
    other_tok_idx = other_chain_token_indices[other_chain]
    
    # Extract PAE submatrix for protein-other chain interaction
    pae_submatrix = pae_matrix_full[np.ix_(protein_token_indices, other_tok_idx)]
    
    # Calculate statistics
    mean_pae = np.mean(pae_submatrix)
    min_pae = np.min(pae_submatrix)
    num_good_pae = np.sum(pae_submatrix < pae_cutoff)
    
    # Get plddt values for the other chain residues (with bounds checking)
    try:
        valid_indices = [other_residues[i]['atom_num'] - 1 
                        for i in range(min(len(other_residues), len(other_tok_idx)))
                        if 0 <= other_residues[i]['atom_num'] - 1 < len(atom_plddts)]
        other_plddt = atom_plddts[valid_indices] if valid_indices else np.array([])
        mean_plddt = np.mean(other_plddt) if len(other_plddt) > 0 else 0.0
    except (IndexError, KeyError):
        mean_plddt = 0.0
    
    num_prot_res = len(protein_token_indices)
    num_other_res = len(other_tok_idx)
    
    # Calculate number of "good" interactions for d0
    good_other_res = np.sum(np.any(pae_submatrix < pae_cutoff, axis=0))
    good_prot_res = np.sum(np.any(pae_submatrix < pae_cutoff, axis=1))
    
    n0 = good_prot_res + good_other_res
    if n0 < 1:
        n0 = 1
    
    d0 = calc_d0(n0, 'protein')
    
    # Calculate ipSAE
    ipsae_by_res = np.zeros(num_prot_res)
    ipsae_max = 0.0
    
    for i in range(num_prot_res):
        pae_row = pae_submatrix[i, :]
        good_mask = pae_row < pae_cutoff
        
        if np.any(good_mask):
            good_paes = pae_row[good_mask]
            n0_res = np.sum(good_mask)
            d0_res = calc_d0(n0_res, 'protein')
            ptm_scores = ptm_func_vec(good_paes, d0_res)
            ipsae_by_res[i] = np.mean(ptm_scores)
            if ipsae_by_res[i] > ipsae_max:
                ipsae_max = ipsae_by_res[i]
    
    # Calculate ipSAE using all good PAE values
    good_paes_all = pae_submatrix[pae_submatrix < pae_cutoff]
    if len(good_paes_all) > 0:
        ipsae_all = np.mean(ptm_func_vec(good_paes_all, d0))
    else:
        ipsae_all = 0.0
    
    results.append({
        'protein_chain': protein_chain,
        'ligand_chain': other_chain,
        'ligand_name': 'peptide',
        'ligand_type': 'peptide',
        'num_prot_res': num_prot_res,
        'num_lig_atoms': num_other_res,
        'ipsae': ipsae_all,
        'ipsae_max': ipsae_max,
        'mean_pae': mean_pae,
        'min_pae': min_pae,
        'num_good_pae': num_good_pae,
        'mean_plddt': mean_plddt
    })
    
    # Write to output file
    outstring = (f"{protein_chain:8} {other_chain:7} {'peptide':11} {'peptide':11} "
                 f"{num_prot_res:10d}  {num_other_res:11d}  "
                 f"{ipsae_all:8.6f}   {ipsae_max:8.6f}   "
                 f"{mean_pae:8.4f}   {min_pae:8.4f}   "
                 f"{num_good_pae:10d}  {mean_plddt:9.2f}  {pdb_stem}\n")
    OUT.write(outstring)
    
    print(f"\nProtein {protein_chain} - Chain {other_chain} (peptide/protein):")
    print(f"  Type: Peptide/protein")
    print(f"  Number of protein residues: {num_prot_res}")
    print(f"  Number of other chain residues: {num_other_res}")
    print(f"  ipSAE: {ipsae_all:.6f}")
    print(f"  ipSAE_max (per residue): {ipsae_max:.6f}")
    print(f"  Mean PAE: {mean_pae:.4f}")
    print(f"  Min PAE: {min_pae:.4f}")
    print(f"  Number of PAE < {pae_cutoff}: {num_good_pae}")
    print(f"  Mean chain pLDDT: {mean_plddt:.2f}")

OUT.close()

print(f"\n{'='*80}")
print(f"Results written to: {file_path}")
print(f"{'='*80}")
