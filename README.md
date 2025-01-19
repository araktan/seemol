I am doing this not because it's easy, but because I thought it would be easy.

The goal of the project is to make **image recognition** for **molecular structures**: i.e. supply an image of molecular structure and getting a **SMILES string**. 

The first difficult bit: the challenge appears to be to infer the molecular graph from an image. The most straightforward and least rewarding approach was to try to have a bw-image and produce an adjacency matrix. This would then be input and the target. However adjacency matrices are permutable and are of varying size. Zero padding is a problem since most of the data would be artificialy added zeros.  Variable input data is yet another can of worms.

Now the next bit is managing the layer output dims and using an appropriate sequence of information. Next is identifying ways of encoding a graph in a more consistent way.

Training: Switching to iterables is another thing to do to save on memory. It seems to be possible to run pytorch backend on my gpu however there are some memory issues to debug.