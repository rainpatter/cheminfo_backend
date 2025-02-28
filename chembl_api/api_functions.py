from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
import io
import base64


def retrieve_mol_data_from_smiles(smiles):
    similarity = new_client.similarity
    results = list(similarity.filter(smiles=smiles, similarity=70).only(['molecule_chembl_id', 'pref_name', 'similarity', 'molecule_structures']))
    # for i in range(len(results)):
    #     mol = results[i]['molecule_structures']['molfile']
    #     m = Chem.MolFromMolBlock(mol)
    #     img = Draw.MolToImage(m)
    #     print(img)
    #     results[i] = dict(results[i], image=img)
    # print([result.keys() for result in results])
    return results

retrieve_mol_data_from_smiles('C=O')