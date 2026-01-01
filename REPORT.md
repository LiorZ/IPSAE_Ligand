# IPSAE Ligand Adaptation Report

## Overview

This report describes the adaptation of the IPSAE (interprotein Structural Alignment Error) scoring function for small molecule ligand-protein complexes in AlphaFold3 models.

## Background

The original `ipsae.py` script was designed to calculate ipSAE scores for pairwise protein-protein interactions in AlphaFold2 and AlphaFold3 models. The scoring is based on the Predicted Aligned Error (PAE) matrix from AlphaFold predictions, which measures the confidence in the relative positions of residue pairs.

## Problem Statement

The goal was to adapt the IPSAE scoring for protein-ligand interactions where:
- Chain A represents the protein
- Chains B, C, D, etc. represent small molecule ligands
- An ipSAE score is calculated for each protein-ligand pair

## Implementation

### New Script: `ipsae_ligand.py`

A new Python script `ipsae_ligand.py` was created specifically for protein-ligand interactions in AlphaFold3 models.

### Key Features

1. **Automatic Chain Classification**
   - Chains with sequential residue numbers (standard amino acids, modified residues) are identified as protein chains
   - Chains with `.` as residue sequence number (HETATM records in mmCIF) are identified as true ligand chains
   - Other protein/peptide chains are also identified and scored

2. **Flexible Protein Chain Selection**
   - By default, uses the first protein chain found (typically chain A)
   - Optional command-line argument to specify which chain is the protein

3. **Comprehensive Output**
   - ipSAE score for each protein-ligand pair
   - ipSAE_max (maximum per-residue score)
   - Mean and minimum PAE values
   - Count of good PAE pairs (below cutoff)
   - Mean pLDDT for ligand atoms
   - Distinction between true ligands and peptide chains

### Technical Implementation

#### Chain Detection
The script performs two passes over the mmCIF file:
1. First pass: Collect all atoms and determine chain types
2. Second pass: Categorize atoms into protein residues and ligands

#### PAE Matrix Processing
- Extracts the PAE submatrix for protein-ligand interactions
- Uses protein token indices (rows) vs ligand token indices (columns)
- Handles the token-based structure of AF3 PAE matrices

#### ipSAE Calculation
For each protein-ligand pair:
1. Count protein residues and ligand atoms with good PAE values
2. Calculate d0 parameter based on the number of interacting elements
3. Compute PTM-like scores for each PAE value below cutoff
4. Report both overall ipSAE (mean of all good scores) and ipSAE_max (best per-residue score)

## Usage

```bash
python ipsae_ligand.py <pae_json_file> <mmcif_file> <pae_cutoff> <dist_cutoff> [protein_chain]

# Examples:
python ipsae_ligand.py fold_aurka_0_tpx2_0_full_data_0.json fold_aurka_0_tpx2_0_model_0.cif 10 10
python ipsae_ligand.py fold_aurka_0_tpx2_0_full_data_0.json fold_aurka_0_tpx2_0_model_0.cif 10 10 A
```

## Example Output

Testing with the AURKA/TPX2/ATP/Mg complex:

```
Reading structure from Example/fold_aurka_0_tpx2_0_model_0.cif...
Found protein chain: A with 276 residues
Found true ligand chains: ['C', 'D', 'E']
  Chain C: 31 atoms, residue: ATP
  Chain D: 1 atoms, residue: MG
  Chain E: 1 atoms, residue: MG
Found other protein/peptide chains: ['B']
  Chain B: 43 residues

Protein A - Ligand C (ATP):
  Type: Small molecule ligand
  ipSAE: 0.881437
  ipSAE_max (per residue): 0.533034
  Mean PAE: 3.4472

Protein A - Ligand D (MG):
  Type: Small molecule ligand
  ipSAE: 0.887526
  ipSAE_max (per residue): 0.519068
  Mean PAE: 3.2966

Protein A - Ligand E (MG):
  Type: Small molecule ligand
  ipSAE: 0.875538
  ipSAE_max (per residue): 0.389736
  Mean PAE: 3.4534

Protein A - Chain B (peptide/protein):
  Type: Peptide/protein
  ipSAE: 0.727848
  ipSAE_max (per residue): 0.448952
  Mean PAE: 5.7617
```

## Output File Format

The output file (`*_ligand_XX_XX.txt`) contains:
- Header with protein chain and cutoff information
- One line per protein-ligand pair with the following columns:
  - Protein: Protein chain identifier
  - Ligand: Ligand chain identifier
  - LigandName: Residue name (ATP, MG, etc.) or "peptide"
  - Type: "ligand" or "peptide"
  - NumProtRes: Number of protein tokens
  - NumLigAtoms: Number of ligand tokens
  - ipSAE: Overall ipSAE score
  - ipSAE_max: Maximum per-residue ipSAE
  - MeanPAE: Mean PAE value
  - MinPAE: Minimum PAE value
  - NumGoodPAE: Count of PAE values below cutoff
  - MeanPlddt: Mean pLDDT for ligand
  - Model: Input filename

## Interpretation

- **ipSAE scores** range from 0 to 1, with higher values indicating better predicted interaction
- **ipSAE > 0.8** typically indicates high confidence in the protein-ligand interaction
- **Low mean PAE** (e.g., < 5 Å) suggests accurate relative positioning
- **High pLDDT** for ligand atoms indicates confident structure prediction

## Files Modified/Created

1. **ipsae_ligand.py** - New script for protein-ligand ipSAE scoring
2. **README.md** - Updated with documentation for the new script
3. **REPORT.md** - This report explaining the implementation

## Limitations

1. Currently only supports AlphaFold3 mmCIF format (not AF2 PDB or Boltz1)
2. Distance cutoff parameter is not currently used (future enhancement)
3. By-residue output file not yet implemented for ligand interactions

## Future Enhancements

1. Add support for distance-based filtering
2. Add by-residue output for detailed analysis
3. Support for multiple protein chains
4. PyMOL script generation for visualization
