import pandas as pd
df = pd.read_csv("data/interim/intact_pairs_with_pdb_contacts_topk_diverse.tsv", sep="\t")

# 上位の“過密”候補（結晶接触や巨大会合体由来の可能性）
print(df.sort_values("n_atom_contacts", ascending=False)
        .head(10)[["uniprot_a","uniprot_b","pdb_id","chain_id_a","chain_id_b","n_atom_contacts","min_atom_distance"]])

# “点接触”っぽい（距離は短いのに原子数が少ない）
sus = df[(df["min_atom_distance"]<1.5) & (df["n_atom_contacts"]<20)]
print("suspicious small-contact hits:", len(sus))

# ペアごとのPDB多様性（1.0に近いほど良い）
g = df.groupby(["uniprot_a","uniprot_b"])
print("mean distinct PDB count per pair:", (g["pdb_id"].nunique()/g.size()).mean())
