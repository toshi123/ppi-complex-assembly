from ppi_complex.db import HumanProteomeDB
db = HumanProteomeDB('data/processed/human_proteome.sqlite', 'data/processed/human_residues.parquet')

# result = db.search_by_gene('LIPA')
# uniprot_id = result[0]['uniprot_id']
# uniprot_id = 'P01033'
uniprot_id = 'P08100'

protein = db.get_protein(uniprot_id, include_go=True, include_locations=True)
print('=== 基本情報 ===')
for name, value in protein.items():
    print(f"  {name}: {value}")
# print(f"UniProt ID: {protein['uniprot_id']}")
# print(f"Description: {protein['description']}")
# print(f"Length: {protein['length']} aa")
# print(f"Surface: {protein['n_surface']}, Buried: {protein['n_buried']}")
# for go in protein["go"]:
#     print(f'GO terms: {go["category"]}: {go["go_term"]}')

print('\n=== 遺伝子名 ===')
for g in db.get_gene_names(uniprot_id):
    print(f"  {g['gene_name']} ({g['name_type']})")

print('\n=== 局在（分類済み） ===')
print(protein.get('locations', []))

print('\n=== 詳細局在（UniProt） ===')
for loc in db.get_uniprot_subcellular_locations(uniprot_id):
    print(f"  {loc}")

print('\n=== シグナル/トランジットペプチド ===')
for r in db.get_uniprot_regions(uniprot_id):
    if 'peptide' in r['region_type']:
        print(f"  {r['region_type']}: {r['start']}-{r['end']}")

print('\n=== 発現パターン ===')
print(db.get_uniprot_expression(uniprot_id))

print('\n=== PPIパートナー ===')
interfaces = db.get_interfaces(uniprot_id)
for partner in interfaces:
    print(f"  {partner}:")
    print(f"    pdb_id: {interfaces[partner]['pdb_id']}")
    print(f"    chains: {interfaces[partner]['chains']}")
    print(f"    residues: {interfaces[partner]['residues']}")
    print(f"    n_residues: {interfaces[partner]['n_residues']}")
    print(f"    bsa_total: {interfaces[partner]['bsa_total']}")
    print(f"    n_atom_contacts: {interfaces[partner]['n_atom_contacts']}")
    print(f"    min_atom_distance: {interfaces[partner]['min_atom_distance']}")
    print(f"    stats:")
    if interfaces[partner]['stats']:
        for key, value in interfaces[partner]['stats'].items():
            print(f"      {key}: {value}")
    else:
        print("      None")

print('\n=== 膜貫通トポロジー ===')
transmembrane = db.get_transmembrane_topology(uniprot_id)
for tm in transmembrane['transmembrane']:
    print(f"  Helical region: {tm['start']}-{tm['end']}")
for topo in transmembrane['topology']:
    print(f"  Topological domain: {topo['start']}-{topo['end']}: {topo['side']}: {topo['description']}")

print('\n=== ジスルフィド結合 ===')
for bond in db.get_uniprot_disulfide_bonds(uniprot_id):
    print(f"  {bond['position1']}-{bond['position2']}")

print('\n=== 領域特徴 ===')
for feature in db.get_uniprot_regions(uniprot_id):
    print(f"  {feature['region_type']}: {feature['start']}-{feature['end']}")

print('\n=== 残基特徴 ===')
for feature in db.get_uniprot_residue_features(uniprot_id):
    print(feature)
    print(f"  {feature['feature_type']}: {feature['position']}: {feature['description']}")

print('\n=== ディスオーダ ===')
for residue in db.get_disordered_residues(uniprot_id):
    print(residue)

print('\n=== 全残基 ===')
for residue in db.get_residues(uniprot_id):
    print(residue)