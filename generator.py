# generates image set from SMILES

from rdkit import Chem
import pandas as pd
import numpy as np


def get_mol_image(smiles_string, img_filename, save_dir):
    # TODO: needs assertion that the string is good   
    mol = Chem.MolFromSmiles(smiles_string)
    mol_img = Chem.Draw.MolToImage(mol)
    mol_img.save(f'{save_dir}/{img_filename}.jpg')
    return


def main():
    df = pd.read_excel('data/SMILES_ID.xlsx', index_col=0)
    for i, sml in enumerate(df['SMILES']):
        try:
            get_mol_image(sml, df['ID'][i], 'imgs')
        except:
            continue


if __name__ == '__main__':
    main()