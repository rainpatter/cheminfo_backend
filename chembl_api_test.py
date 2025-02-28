# https://github.com/chembl/notebooks/blob/main/ChEMBL_webresource_client_examples.ipynb


from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt



similarity = new_client.similarity
results = similarity.filter(smiles="CO[C@@H](CCC#C\C=C/CCCC(C)CCCCC=C)C(=O)[O-]", similarity=70).only(['molecule_chembl_id', 'pref_name', 'similarity', 'molecule_structures'])

mol = results[0]['molecule_structures']['molfile']
m = Chem.MolFromMolBlock(mol)

img = Draw.MolToImage(m)

print(img)

imgplot = plt.imshow(img)
plt.show()


