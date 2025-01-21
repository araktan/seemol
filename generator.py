import sys
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
import numpy as np

def get_x_y(smiles_string, image_dims=(60, 60)):
    # TODO: needs assertion that the string is good   
    mol = Chem.MolFromSmiles(smiles_string)
    img = Draw.MolToImage(mol, size=image_dims)
    thresh = 200
    fn = lambda x : 0 if x > thresh else 255
    r = img.convert('L').point(fn, mode='1')
    #r = img.convert('L')
    adj = Chem.GetAdjacencyMatrix(mol)
    padded_adj = np.zeros(image_dims)
    for index, bond_val in np.ndenumerate(adj):
        padded_adj[index] += bond_val
    return np.expand_dims(np.array(r), -1), np.expand_dims(padded_adj, -1)

def main():
    cmd_args = sys.argv
    if len(cmd_args) > 1:   
        print('Setting shape: ', cmd_args[1:])
        output_shape = (int(cmd_args[1]), int(cmd_args[2]))
    else:
        print('Default shape: ', [60, 60])
        output_shape = (60, 60)


    # generte data
    df = pd.read_csv('data/SMILES.csv')

    x_data, y_data = [],[]
    for SMILES in df['SMILES']:
        try:
            img, adj_matrix = get_x_y(SMILES, image_dims=output_shape)
            x_data.append(img)
            y_data.append(adj_matrix)
        except:
            continue

    x_data, y_data = np.array(x_data), np.array(y_data)
    
    with open('test_np_X.npy', 'wb') as f:
        np.save(f, x_data, allow_pickle=True)
    
    with open('test_np_Y.npy', 'wb') as f:
        np.save(f, y_data, allow_pickle=True)
    
    return

if __name__ == '__main__':
    main()

