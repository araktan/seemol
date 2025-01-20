I am doing this not because it's easy, but because I thought it would be easy.

The goal of the project is to make **image recognition** for **molecular structures**: i.e. supply an image of molecular structure and getting a **SMILES string**. 

The first difficult bit: the challenge appears to be to infer the molecular graph from an image. The most straightforward and least rewarding approach was to try to have a bw-image and produce an adjacency matrix. This would then be input and the target. However adjacency matrices are permutable and are of varying size. Zero padding is a problem since most of the data would be artificialy added zeros.  Variable input data is yet another can of worms.

First win would be getting a similar graph from predicted adj matrix.