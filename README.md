# ClonalTree

**Reconstructing the evolutionary history of a BCR lineage with minimum spanning tree and clonotype abundances**

ClonalTree is a new algorithm to reconstruct BCR lineage trees that incorporates genotype abundance into a minimum spanning tree to infer maximum parsimony trees.

**CONTACT**  
  E-mail: 
  juliana.silva_bernardes@sorbonne-universite.fr 
  
## Inputs
 
  * The alignment of clonally related BCR heavy chain sequences in FASTA format. The naive B cell heavy chain sequence must be included in the alignment, named as "naive". There are ways to provide the genotype abundance for each sequence :
    * By the sequence ID repetition in the alignment, or
    * By integrating the abundance of each sequence in its ID, for instance, the sequence SeqX with an abundance of 12 will have the following ID:
      >SeqX@12
  * See [example input files](https://github.com/julibinho/ClonalTree/Example/Input)

## Outputs

  * ClonalTree returns:

    - [repertoire_name].nk : the reconstructed BCR lineage tree in [newick format](https://en.wikipedia.org/wiki/Newick_format) 

    - [repertoire_name].nk.csv :  a table in csv format, containing the parent relationship and cost.
  * See [example output files](https://github.com/julibinho/ClonalTree/Example/Output)

     
      
## Requirements 

  * numpy :
      ```
      $ conda install numpy
      ```
      or 
      ```
      $ pip install numpy
      ```

  * Biopython
      ```
      $ pip install biopython
      ```

  * ete3 :
      ```
      $ pip install ete3
      ```

  * networkx (for the evaluation) :
      ```
      $ pip install networkx
      ```


## Using ClonalTree 
   The command line for launching the ClonalTree is:

  ```
  $ python clonalTree.py  -i [seq_alignment_file] -o [seq_alignment_file] [...options]

  ```
### Arguments

  * -a 0, if considering abundance, otherwise , -a 1
  * -r 0, if performing revision, otherwise -r 1
  * -t 0, if performing trimming tree, otherwise -t 1


  For instance the following command can be run in the src/ folder:
  ```
  $ python clonalTree.py  -i ../Example/Input/simulation200.fasta -o ../Example/Output/clonalTree.abRT.nk -a 0 -r 0 -t 0
  ```
                      
  Output files will be placed as such:
  ```
  ~Example/Output/[seq_alignment_file].nk
                  [seq_alignment_file].nk.csv
 ```
 [seq_alignment_file] is the alignement of the clonaly related B cells receptor heavly chain seuqneces.

## License, Patches, and Ongoing Developements

  * The program is distributed under the ....
  * [Feature requests and open issues](https://github.com/julibinho/ClonalTree/issues).
