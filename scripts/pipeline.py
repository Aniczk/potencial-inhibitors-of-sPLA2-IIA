import os
import re
import pandas as pd
import matplotlib.pyplot as plt


def create_folder(dir):
    """Create a folder called <dir> if it doesn't exist."""
    if not os.path.exists(dir):
        os.mkdir(dir)

def get_receptor_seq(line):
    """Return the receptor sequence"""
    if line_as_list[0]=='receptor_sequence':
        return line_as_list[2]

def get_index(pattern, sequence):
    """Return the amino acid index in the sequence corresponding to the beginning of the entered pattern."""
    p = re.compile(pattern)
    p2 = p.findall(sequence) 
    return sequence.index(p2[0])

def read_readme_file(file):
    """Parse the README file.
    Return the amino acid indexes in the sequence corresponding to the beginning of the sites:
    AALSYGFYG, CG*GG and D*CC**HD."""
    index,ligand = None,None
    with open(file) as fp:
        for line in fp:
            line_as_list = line.strip().split(" ")
            if line_as_list[0]=='receptor_sequence':
                seq = line_as_list[2]
                index_ca_loop = get_index("CG\wGG", seq)
                index_active_site = get_index("D\wCC\w\wHD", seq) 
                if 'AALSYGFYG' in seq:
                    index_aalsygfyg = seq.index('AALSYGFYG')
            if line_as_list[0]=='ligand_sequence':
                ligand = line_as_list[2]
    return index_aalsygfyg, index_ca_loop, index_active_site, ligand


def get_charts(index, length):
    """Return the list of characters corresponding to the specified indices"""
    L = []
    for i in range(index+1,index+1+length):
        L.append('A'+str(i))
    return L

def create_dataframe(cols):
    """Create a dataframe with columns <cols>"""
    data = {}
    for i in cols:
        data[i] = []
    return pd.DataFrame(data)

def plot_contacts(df,t, peptide):
    """Save the plot in .jpg format showing the frequency of peptide contacts with AALSYGFYG site of phospholipase."""
    df.plot(kind='scatter', x='receptor', y=peptide, c='contact', title=t)
    directory = os.path.join("..","plots","AALSYGFYG_contacts_scatter")
    create_folder(directory)
    plt.savefig(os.path.join("..","plots", "AALSYGFYG_contacts_scatter",t+'.jpg'))
    plt.close()

def plot_number_of_interactions(df):
    """Save the plot in .jpg format showing how many models of phospholipase the given peptide was in contact with."""
    df = df.sort_values('interactions', ascending=False)
    df.plot.bar(legend=False)
    plt.ylabel("Number of interactions")
    plt.xlabel("Ligands")
    directory = os.path.join("..","plots")
    create_folder(directory)
    plt.savefig(os.path.join("..","plots",'AALSYGFYG_number_of_interactions.jpg'),bbox_inches='tight')
    plt.close()

def read_map_file(file,ind_1,ind_2,ind_3,peptide, L2,L3,L4):
    """Parse .map file.
    Return dataframes and lists of amino acid names and numbers that correspond to sequence fragments:
    AALSYGFYG - L1
    CG*GG - L2
    D*CC**HD - L3"""
    df = create_dataframe(['receptor', peptide, 'contact'])
    df2 = create_dataframe(['receptor_position', peptide])
    df3 = create_dataframe(['receptor_position', peptide])
    df4 = create_dataframe(['receptor_position', peptide])
    
    AALSYGFYG_ch = get_charts(ind_1, len('AALSYGFYG'))
    ca_ch = get_charts(ind_2, 5)
    active_site_ch = get_charts(ind_3, 8)
    
    if AALSYGFYG_ch not in L2: L2.append(AALSYGFYG_ch)
    if ca_ch not in L3: L3.append(ca_ch)
    if active_site_ch not in L4: L4.append(active_site_ch)

    with open(file) as fp:
        for line in fp:
            line_as_list = line.strip().split(" ")

            # inhibitor and AALSYGFYG contact
            if line_as_list[0] in AALSYGFYG_ch:
                row = {'receptor': str(line_as_list[0]),
                peptide: str(line_as_list[1]),
                'contact': float(line_as_list[2])}

                row2 = {'receptor_position': str(line_as_list[0]),
                peptide: float(line_as_list[2])}

                df = df.append(row, ignore_index=True)
                df2 = df2.append(row2, ignore_index=True)

            # inhibitor and Ca binding site contact
            if line_as_list[0] in ca_ch:
                row3 = {'receptor_position': str(line_as_list[0]),
                peptide: float(line_as_list[2])}

                df3 = df3.append(row3, ignore_index=True)

            # inhibitor and active site contact
            if line_as_list[0] in active_site_ch:
                row4 = {'receptor_position': str(line_as_list[0]),
                peptide: float(line_as_list[2])}

                df4 = df4.append(row4, ignore_index=True)

    return df, df2, df3, df4, L2, L3, L4

def plot_mean(df, filename):
    """Save the plot in .jpg format showing mean contact for each inhibitor."""
    df.plot(kind="bar")
    plt.ylabel("Mean contact")
    plt.xlabel("Ligands")
    create_folder(os.path.join("..","plots"))
    plt.savefig(os.path.join("..","plots",filename),bbox_inches='tight')
    plt.close()

def region_as_list_of_symbols(region):
    """Return list of symbols corresponding to the <region>."""
    symbols = []
    for i in region:
        if i not in symbols:
            symbols.append(i)
        else:
            count = symbols.count(i)
            symbols.append(i+'\''*count)
    return symbols


def dataframe_with_contacts(df, L, symbols):
    """Return dataframe with contacts."""
    df_new = create_dataframe([])
    for i, j, k in zip(L[0], L[1], symbols):
        v = (df.loc[i] + df.loc[j])/2
        df_new[k] = v
    return df_new

def dataframe_with_Ca_binding_site_contact(df):
    """"Return dataframe with Ca binding site contact."""
    df_new = create_dataframe([])
    for i, j, k in zip(['A28', 'A29', 'A30', 'A31', 'A32'],['A48', 'A49', 'A50', 'A51', 'A52'],['C','G','V','G\'','G\'\'']):
        v = (df.loc[i] + df.loc[j])/2
        df_new[k] = v
    return df_new


def plot_positions_contact(datas, sequence, L, region):
    """Save the plot in .jpg format showing mean contact at each inhibitor's position."""
    symbols = region_as_list_of_symbols(sequence)

    result_1 = pd.concat(datas, axis=0)
    df_groupby = result_1.groupby(['receptor_position']).mean()

    df_new = dataframe_with_contacts(df_groupby, L, symbols)    
    plot_mean(df_new, region+'_mean_contact_each_pos.jpg')

    df_mean_rows = df_new.mean(axis = 1)
    plot_mean(df_mean_rows, region+'_mean_contact.jpg')



folder = '..'

sub_folders = [name for name in os.listdir(folder) if os.path.isdir(os.path.join(folder, name))]

def main():

    datas, datas2, datas3, datas4, titles, ligands, L2, L3, L4 = [], [], [], [], [], [], [], [], []
    counter = 0
    
    for i in sub_folders:
        path_map = os.path.join("..",i,"MAPS","cluster_1.map")
        path_readme = os.path.join("..",i,"README")
        if os.path.isfile(path_map) and os.path.isfile(path_readme):
            x = i.split('_')
            ligand_name = x[-1]    
            index_aalsygfyg, index_ca, index_active_site, ligand = read_readme_file(path_readme)
            d1, d2, d3, d4, L2, L3, L4 = read_map_file(path_map,index_aalsygfyg, index_ca, index_active_site, ligand_name, L2, L3, L4)
            datas.append(d1)
            datas2.append(d2)
            datas3.append(d3)
            datas4.append(d4)
            titles.append(i)
            ligands.append(ligand)
    df_ligands = create_dataframe(['interactions'])
    df_positions = create_dataframe(['receptor_position'])

    for i,t in enumerate(titles):
        x = t.split('_')
        ligand_name = x[-1]
        if not (datas[i]['contact']==0).all():
            plot_contacts(datas[i],t,ligand_name)

            if ligands[i] not in list(df_ligands.index):
                df_ligands.loc[ligands[i]] = 1
            else:
                v = int(df_ligands.loc[ligands[i]])
                df_ligands.loc[ligands[i]] = v+1
    plot_number_of_interactions(df_ligands)
    print(datas2)

    plot_positions_contact(datas2, 'AALSYGFYG', L2, 'oligomers_reg')
    plot_positions_contact(datas3, 'CGVGG', L3, 'Ca_binding_site')
    plot_positions_contact(datas4, 'DRCCVTHD', L4, 'active_site')


if __name__ == '__main__':
    main()
