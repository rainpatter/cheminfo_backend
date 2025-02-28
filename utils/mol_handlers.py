import re

def smiles_from_query_string(query_string):
    cleaned_smiles = re.sub('%3','=', query_string)
    return cleaned_smiles