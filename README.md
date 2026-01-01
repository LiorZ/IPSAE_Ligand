# IPSAE
Scoring function for interprotein interactions in AlphaFold2 and AlphaFold3

# Installation
Simply download the Python script ipsae.py. It may be necessary to install the Python numpy package with:

    pip install numpy

# Usage:                                                                                                                                                                                                     
AlphaFold2:

     python ipsae.py <path_to_af2_json_file> <path_to_af2_pdb_file> <pae_cutoff> <dist_cutoff>   
     python ipsae.py RAF1_KSR1_scores_rank_001_alphafold2_multimer_v3_model_4_seed_003.json RAF1_KSR1_unrelaxed_rank_001_alphafold2_multimer_v3_model_4_seed_003.pdb 15 15

AlphaFold3:

     python ipsae.py <path_to_af3_json_file> <path_to_af3_cif_file> <pae_cutoff> <dist_cutoff>                    
     python ipsae.py fold_aurka_tpx2_full_data_0.json fold_aurka_tpx2_model_0.cif 10 10

Boltz1:

     python ipsae.py <path_to_boltz1_pae_npz_file> <path_to_boltz1_cif_file> <pae_cutoff> <dist_cutoff>
     python ipsae.py pae_AURKA_TPX2_model_0.npz  AURKA_TPX2_model_0.cif 10 10


# IPSAE Ligand - Protein-Small Molecule Interactions

For scoring protein-ligand interactions in AlphaFold3 models, use the `ipsae_ligand.py` script:

## Usage

    python ipsae_ligand.py <path_to_af3_json_file> <path_to_af3_cif_file> <pae_cutoff> <dist_cutoff> [protein_chain]
    python ipsae_ligand.py fold_aurka_0_tpx2_0_full_data_0.json fold_aurka_0_tpx2_0_model_0.cif 10 10
    python ipsae_ligand.py fold_aurka_0_tpx2_0_full_data_0.json fold_aurka_0_tpx2_0_model_0.cif 10 10 A

This script:
- Treats chain A (or specified chain) as the protein
- Treats all other chains as ligands (small molecules) or peptides
- Calculates ipSAE scores for each protein-ligand pair

## Output ligand score file

Example output for AURKA (protein chain A) with ATP, two Mg ions, and TPX2 peptide:

    Protein  Ligand  LigandName  Type        NumProtRes  NumLigAtoms  ipSAE      ipSAE_max  MeanPAE    MinPAE     NumGoodPAE  MeanPlddt  Model
    A        C       ATP         ligand             296           31  0.881437   0.533034     3.4472     0.8000         8719      91.93  fold_aurka_0_tpx2_0_model_0
    A        D       MG          ligand             296            1  0.887526   0.519068     3.2966     1.0000          281      91.47  fold_aurka_0_tpx2_0_model_0
    A        E       MG          ligand             296            1  0.875538   0.389736     3.4534     1.3000          281      85.96  fold_aurka_0_tpx2_0_model_0
    A        B       peptide     peptide            296           43  0.727848   0.448952     5.7617     0.9000        11236      78.29  fold_aurka_0_tpx2_0_model_0

## Output columns

**Protein** = Protein chain identifier

**Ligand** = Ligand/other chain identifier

**LigandName** = Residue name of the ligand (e.g., ATP, MG) or "peptide" for protein chains

**Type** = "ligand" for small molecules or "peptide" for protein/peptide chains

**NumProtRes** = Number of protein residues/tokens

**NumLigAtoms** = Number of ligand atoms/tokens

**ipSAE** = ipSAE score calculated using all PAE values below the cutoff

**ipSAE_max** = Maximum per-residue ipSAE score

**MeanPAE** = Mean PAE value for the protein-ligand interaction

**MinPAE** = Minimum PAE value for the protein-ligand interaction  

**NumGoodPAE** = Number of PAE values below the cutoff

**MeanPlddt** = Mean pLDDT value for the ligand atoms

**Model** = Input model filename


# Output chain-chain score file (original ipsae.py)

AlphaFold2 complex of RAF1 (chain A), KSR1 (chain B), and MEK1 (chain C). RAF1 and KSR1 kinase domains make pseudosymmetric dimer, similar to BRAF homodimer. In this model, MEK1 kinase domain interacts directly with RAF1 kinase domain (AF3 places MEK1 on KSR1; both complexes exist in PDB). 

    Chn1 Chn2  PAE Dist  Type   ipSAE    ipSAE_d0chn ipSAE_d0dom  ipTM_af  ipTM_d0chn     pDockQ     pDockQ2    LIS       n0res  n0chn  n0dom   d0res   d0chn   d0dom  nres1   nres2   dist1   dist2  Model
    A    B     15   15   asym  0.647603    0.868775    0.753413    0.570    0.378501      0.2342     0.0142     0.4128     290   1571    553    6.26   12.57    8.29    261     292      92      93   RAF1_KSR1_MEK1_9f755_rank001_iptm0.57_ptm0.54_nrx_4_000
    B    A     15   15   asym  0.641387    0.867977    0.749378    0.570    0.482265      0.2342     0.0143     0.3782     285   1571    543    6.21   12.57    8.22    255     288      95      94   RAF1_KSR1_MEK1_9f755_rank001_iptm0.57_ptm0.54_nrx_4_000
    A    B     15   15   max   0.647603    0.868775    0.753413    0.570    0.482265      0.2342     0.0143     0.3955     290   1571    553    6.26   12.57    8.29    288     292      94      95   RAF1_KSR1_MEK1_9f755_rank001_iptm0.57_ptm0.54_nrx_4_000
    
    A    C     15   15   asym  0.666048    0.810562    0.750777    0.570    0.765550      0.3016     0.1158     0.3642     364   1041    635    6.93   10.71    8.77    270     365     114     105   RAF1_KSR1_MEK1_9f755_rank001_iptm0.57_ptm0.54_nrx_4_000
    C    A     15   15   asym  0.626076    0.818914    0.749619    0.570    0.438184      0.3016     0.0912     0.3470     284   1041    597    6.20   10.71    8.55    313     284      96     112   RAF1_KSR1_MEK1_9f755_rank001_iptm0.57_ptm0.54_nrx_4_000
    A    C     15   15   max   0.666048    0.818914    0.750777    0.570    0.765550      0.3016     0.1158     0.3556     364   1041    635    6.93   10.71    8.77    284     365     114     105   RAF1_KSR1_MEK1_9f755_rank001_iptm0.57_ptm0.54_nrx_4_000
    
    B    C     15   15   asym  0.323178    0.576574    0.413315    0.570    0.542858      0.0000     0.0000     0.1329     352   1316    563    6.83   11.74    8.35    207     356       0       0   RAF1_KSR1_MEK1_9f755_rank001_iptm0.57_ptm0.54_nrx_4_000
    C    B     15   15   asym  0.323743    0.619089    0.457644    0.570    0.287348      0.0000     0.0000     0.1335     288   1316    566    6.24   11.74    8.37    277     289       0       0   RAF1_KSR1_MEK1_9f755_rank001_iptm0.57_ptm0.54_nrx_4_000
    B    C     15   15   max   0.323743    0.619089    0.457644    0.570    0.542858      0.0000     0.0000     0.1332     288   1316    566    6.24   11.74    8.37    289     356       0       0   RAF1_KSR1_MEK1_9f755_rank001_iptm0.57_ptm0.54_nrx_4_000



**Chn1**=first chain

**Chn2**=second chain

**PAE**=PAE cutoff

**Dist**=Distance cutoff for CA-CA contacts

**Type**="asym" or "max"; asym means asymmetric ipTM/ipSAE values; max is maximum of asym values

**ipSAE**=ipSAE value for given PAE cutoff and d0 determined by number of residues in 2nd chain with PAE<cutoff 

**ipSAE_d0chn**=ipSAE calculated with PAE cutoff and d0 = sum of chain lengths

**ipSAE_d0dom**=ipSAE calculated with PAE cutoff and d0 = total number of residues in both chains with any interchain PAE<cutoff

**ipTM_af**=AlphaFold ipTM values. For AF2, this is for whole complex from json file. For AF3, this is symmetric pairwise value from summary json file.   

**ipTM_d0chn**=ipTM (no PAE cutoff) calculated from PAE matrix and d0 = sum of chain lengths 

**pDockQ**=score from pLDDTs from Bryant, Pozotti, and Eloffson. https://www.nature.com/articles/s41467-022-28865-w

**pDockQ2**=score based on PAE, calculated pairwise: from Zhu, Shenoy, Kundrotas, Elofsson. https://academic.oup.com/bioinformatics/article/39/7/btad424/7219714

**LIS**=Local Interaction Score based on transform of PAEs from Kim, Hu, Comjean, Rodiger, Mohr, Perrimon. https://www.biorxiv.org/content/10.1101/2024.02.19.580970v1

**n0res**=number of residues for d0 in ipSAE calculation

**n0chn**=number of residues in d0 in ipSAE_d0chn calculation

**n0dom**=number of residues in d0 in ipSAE_d0dom calculation

**d0res**=d0 for ipSAE

**d0chn**=d0 for ipSAE_d0chn

**d0dom**=d0 for ipSAE_d0dom

**nres1**=number of residues in chain1 with PAE<cutoff with residues in chain2

**nres2**=number of residues in chain2 with PAE<cutoff with residues in chain1

**dist1**=number of residues in chain 1 with PAE<cutoff and dist<cutoff from chain2

**dist2**=number of residues in chain 2 with PAE<cutoff and dist<cutoff from chain1

**pdb**=AlphaFold filename

# Output by-residue score file

    i   AlignChn ScoredChain  AlignResNum  AlignResType  AlignRespLDDT      n0chn  n0dom  n0res    d0chn     d0dom     d0res    ptm_pae   psae_d0chn  psae_d0dom  psae_d0res 
    1       A         B            1           MET           16.32          1571    541      0    12.569     8.210     1.001     0.1362      0.0000      0.0000      0.0000
    2       A         B            2           GLU           22.12          1571    541      0    12.569     8.210     1.001     0.1362      0.0000      0.0000      0.0000
    3       A         B            3           HIS           24.73          1571    541      0    12.569     8.210     1.001     0.1362      0.0000      0.0000      0.0000
    4       A         B            4           ILE           25.78          1571    541      0    12.569     8.210     1.001     0.1362      0.0000      0.0000      0.0000
    5       A         B            5           GLN           25.27          1571    541      0    12.569     8.210     1.001     0.1362      0.0000      0.0000      0.0000
    6       A         B            6           GLY           28.12          1571    541      0    12.569     8.210     1.001     0.1363      0.0000      0.0000      0.0000


**i**=residue in model (from 1 to total number of residues in model)

**AlignChn**=chainid of aligned residue (in PAE calculation)

**ScoredChn**= chainid of scored residues (with PAE less than cutoff)

**AlignResNum**=residue number of aligned residue

**AlignResType**=residue type of aligned residue (three letter code)

**AlignRespLDDT**=plDDT of aligned residue

**n0chn**=number of residues in d0 in ipSAE_d0chn calculation

**n0dom**=number of residues in d0 in ipSAE_d0dom calculation

**n0res**=number of residues for d0 in ipSAE calculation

**d0chn**=d0 for ipSAE_d0chn

**d0dom**=d0 for ipSAE_d0dom

**d0res**=d0 for ipSAE

**pTM_pae**=ipTM calculated from PAE matrix and d0 = sum of chain lengths 

**pSAE_d0chn**=residue-specific ipSAE calculated with PAE cutoff and d0 = sum of chain lengths (n0chn)

**pSAE_d0dom**=residue-specific ipSAE calculated with PAE cutoff and d0 = total number of residues in both chains with any interchain PAE<cutoff (n0dom)

**pSAE**=residue-specific ipSAE value for given PAE cutoff and d0 determined by number of residues in 2nd chain with PAE<cutoff (n0res)
