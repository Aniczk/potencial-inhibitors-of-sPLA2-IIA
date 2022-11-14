import os
import re
import pandas as pd
import matplotlib.pyplot as plt


def create_folder(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

def get_receptor_seq(line):
	if line_as_list[0]=='receptor_sequence':
		return line_as_list[2]

def read_readme_file(file):
	index,ligand = None,None
	with open(file) as fp:
		for line in fp:
			line_as_list = line.strip().split(" ")
			if line_as_list[0]=='receptor_sequence':
				seq = line_as_list[2]
				p = re.compile("CG\wGG")
				ca_loop = p.findall(seq) 
				index_ca_loop = seq.index(ca_loop[0])
				if 'AALSYGFYG' in seq:
					index_aalsygfyg = seq.index('AALSYGFYG')
			if line_as_list[0]=='ligand_sequence':
				ligand = line_as_list[2]
	return index_aalsygfyg, index_ca_loop, ligand


def get_charts(index, length):
	L = []
	for i in range(index+1,index+1+length):
		L.append('A'+str(i))
	return L

def create_dataframe(cols):
	data = {}
	for i in cols:
		data[i] = []
	return pd.DataFrame(data)

def plot_contacts(df,t, peptide):
	df.plot(kind='scatter', x='receptor', y=peptide, c='contact', title=t)
	directory = os.path.join("..","plots","AALSYGFYG_contacts_scatter")
	create_folder(directory)
	plt.savefig(os.path.join("..","plots", "AALSYGFYG_contacts_scatter",t+'.jpg'))
	plt.close()

def plot_number_of_interactions(df):
	df = df.sort_values('interactions', ascending=False)
	df.plot.bar(legend=False)
	plt.ylabel("Number of interactions")
	plt.xlabel("Ligands")
	directory = os.path.join("..","plots")
	create_folder(directory)
	plt.savefig(os.path.join("..","plots",'AALSYGFYG_number_of_interactions.jpg'),bbox_inches='tight')
	plt.close()

def read_map_file(file,ind_1,ind_2,peptide):
	df = create_dataframe(['receptor', peptide, 'contact'])
	df2 = create_dataframe(['receptor_position', peptide])
	df3 = create_dataframe(['receptor_position', peptide])
	AALSYGFYG_ch = get_charts(ind_1, len('AALSYGFYG'))
	ca_ch = get_charts(ind_2, 5)
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

	return df, df2, df3

def plot_mean(df, filename):
	df.plot(kind="bar")
	plt.ylabel("Mean contact")
	plt.xlabel("Ligands")
	create_folder(os.path.join("..","plots"))
	plt.savefig(os.path.join("..","plots",filename),bbox_inches='tight')
	plt.close()

def dataframe_with_AALSYGFYG_contact(df):
	df_new = create_dataframe([])
	for i, j, k in zip(['A17', 'A18', 'A19', 'A20', 'A21', 'A22', 'A23', 'A24', 'A25'],['A37', 'A38', 'A39', 'A40', 'A41', 'A42', 'A43', 'A44', 'A45'],['A','A\'','L','S','Y','G','F','Y\'','G\'']):
		v = (df.loc[i] + df.loc[j])/2
		df_new[k] = v
	return df_new

def dataframe_with_Ca_binding_site_contact(df):
	df_new = create_dataframe([])
	for i, j, k in zip(['A28', 'A29', 'A30', 'A31', 'A32'],['A48', 'A49', 'A50', 'A51', 'A52'],['C','G','V','G\'','G\'\'']):
		v = (df.loc[i] + df.loc[j])/2
		df_new[k] = v
	return df_new

def plot_positions_contact(datas, region):
	result_1 = pd.concat(datas, axis=0)
	df_groupby = result_1.groupby(['receptor_position']).mean()

	if region == "AALSYGFYG":
		df_new = dataframe_with_AALSYGFYG_contact(df_groupby)
	else: # region == Ca_binding_site
		df_new = dataframe_with_Ca_binding_site_contact(df_groupby)
		
	plot_mean(df_new, region+'_mean_contact_each_pos.jpg')

	df_mean_rows = df_new.mean(axis = 1)
	plot_mean(df_mean_rows, region+'_mean_contact.jpg')




folder = '..'

sub_folders = [name for name in os.listdir(folder) if os.path.isdir(os.path.join(folder, name))]

def main():

	datas, datas2, datas3, titles, ligands = [], [], [], [], []
	counter = 0
	
	for i in sub_folders:
		path_map = os.path.join("..",i,"MAPS","cluster_1.map")
		path_readme = os.path.join("..",i,"README")
		if os.path.isfile(path_map) and os.path.isfile(path_readme):
			x = i.split('_')
			ligand_name = x[-1]	
			index_aalsygfyg, index_ca, ligand = read_readme_file(path_readme)
			d1,d2,d3 = read_map_file(path_map,index_aalsygfyg, index_ca, ligand_name)
			datas.append(d1)
			datas2.append(d2)
			datas3.append(d3)
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

	plot_positions_contact(datas2, 'AALSYGFYG')
	plot_positions_contact(datas3, 'Ca')

if __name__ == '__main__':
	main()
