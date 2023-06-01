import prolif as plf
import MDAnalysis as mda

prot = mda.Universe(
    '/home/receptor/Work/Projects/CR3/production/CD11b-input/6/CD11b-f6.pdb')
protein = plf.Molecule.from_mda(prot)
ligand = plf.sdf_supplier('/home/receptor/test.sdf')


fp = plf.Fingerprint()
fp.run_from_iterable(ligand, protein)
df = fp.to_dataframe()
print(df)
