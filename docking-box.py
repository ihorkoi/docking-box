import oddt
from oddt import interactions
import pandas as pd
import numpy as np


def get_hbonds(protein, ligand):
    prot_atoms, lig_atoms, strict = oddt.interactions.hbonds(
        protein, ligand, cutoff=3.4)

    # Get the protein atom involved in the first H-bond
    resis = []
    resns = []
    chains = []
    atomns = []

    for p_at in lig_atoms:
        resis.append(p_at['resnum'])
        resns.append(p_at['resname'])
        # Map back to the atom in the protein object
        main_struct_at = protein.atoms[p_at['id']]
        # Get the underlying rdkit Atom object:
        main_struct_at_rdk = main_struct_at.Atom
        # Get the chain ID and PDB-format atom name
        chains.append(main_struct_at_rdk.GetPDBResidueInfo().GetChainId())
        atomns.append(main_struct_at_rdk.GetPDBResidueInfo().GetName())

    colnames = ['residue_number', 'residue_name', 'chain', 'atom']
    prot_df = pd.DataFrame(data=np.array(
        [resis, resns, chains, atomns]).T, columns=colnames)
    return prot_df


def get_hydrophobic(protein, ligand):
    prot_atoms = oddt.interactions.pi_stacking(
        protein, ligand, cutoff=3.5)

    resis = []
    resns = []
    chains = []
    atomns = []
    return prot_atoms


''' for p_at in lig_atoms:
     resis.append(p_at['resnum'])
     resns.append(p_at['resname'])
     # Map back to the atom in the protein object
     main_struct_at = protein.atoms[p_at['id']]
     # Get the underlying rdkit Atom object:
     main_struct_at_rdk = main_struct_at.Atom
     # Get the chain ID and PDB-format atom name
     chains.append(main_struct_at_rdk.GetPDBResidueInfo().GetChainId())
     atomns.append(main_struct_at_rdk.GetPDBResidueInfo().GetName())
 colnames = ['residue_number', 'residue_name', 'chain', 'atom']
 prot_df = pd.DataFrame(data=np.array(
     [resis, resns, chains, atomns]).T, columns=colnames)'''


def create_3D_box(pocket):

    for mol in oddt.toolkit.readfile('mol2', f'{pocket}'):
        coord_matrix = []
    for atom in mol.atoms:
        coord_matrix.append([coord for coord in atom.coords])

    x = []
    y = []
    z = []

    for atom in coord_matrix:
        x.append(atom[0])
        y.append(atom[1])
        z.append(atom[2])

    min_x = min(x)
    max_x = max(x)

    min_y = min(y)
    max_y = max(y)

    min_z = min(z)
    max_z = max(z)

    return [[min_x, min_y, min_z], [min_x, min_y, max_z], [min_x, max_y, min_z], [min_x, max_y, max_z], [max_x, min_y, min_z], [max_x, min_y, max_z], [max_x, max_y, min_z], [max_x, max_y, max_z]]


def measure_molecule(matrix):
    x = []
    y = []
    z = []

    for atom in matrix:
        x.append(atom[0])
        y.append(atom[1])
        z.append(atom[2])

    min_x = min(x)
    max_x = max(x)

    min_y = min(y)
    max_y = max(y)

    min_z = min(z)
    max_z = max(z)

    return min_x, max_x, min_y, max_y, min_z, max_z


def is_atom_in_box(atom, box):
    if box[0][0] <= atom[0] <= box[-1][0] and box[0][1] <= atom[1] <= box[-1][1] and box[0][2] <= atom[2] <= box[-1][2]:
        return True
    return False


box = create_3D_box(
    '/home/receptor/Work/Projects/CR3/production/CD11b-input/3/pocket0.mol2')

protein = next(oddt.toolkit.readfile(
    'pdb', '/home/receptor/Work/Projects/CR3/production/CD11b-input/6/CD11b-f6.pdb'))
protein.protein = True
for mol in oddt.toolkit.readfile('sdf', '/home/receptor/test.sdf'):
    coord_matrix = []
    # for atom in mol.atoms:
    # coord_matrix.append([coord for coord in atom.coords])
    # print(coord_matrix)

    # print(is_atom_in_box([coord for coord in atom.coords], box))
    # print(get_hbonds(protein=protein, ligand=mol))
    print(get_hbonds(protein=protein, ligand=mol))
