import os
location = os.getcwd()

x = os.listdir("../ligandy")
standard_aa = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
y = os.listdir("../ligandy_only_atoms")


def make_files_with_ATOM_lines_only(ligands_list):
    for ligand in ligands_list:
        in_file = '../ligandy/'+ligand
        out_file = '../ligandy_only_atoms/only_atoms_'+ligand
        with open(in_file) as in_pdb_file:
            with open(out_file, "w") as out_pdb_file:
                for line in in_pdb_file:
                    if line[:4] == 'ATOM':
                        out_pdb_file.write(line)
        in_pdb_file.close()
        out_pdb_file.close()

def standard_aa_checking(ligands_list):
    from biopandas.pdb import PandasPdb
    for ligand in ligands_list:
        in_file = '../ligandy_only_atoms/'+ligand
        p = PandasPdb().read_pdb(in_file)
        data_atom = p.df['ATOM']
        res_names = data_atom['residue_name']
        L = list(res_names)
        for i in L:
            if i not in standard_aa:
                print(i)

def main(x,y):
    make_files_with_ATOM_lines_only(x)
    standard_aa_checking(y)

if __name__ == '__main__':
    main(x,y)