
### The Correctness Of Ancestral Reconstruction (COAR)

```
      $ python COAR.py -a [groundTruth].nk  -b [inferredTree].nk  -f [seq_alignment_file].fasta -o [output_file] 
```



### The Most Recent Common Ancestor (MRCA) 

```
      $ python MRCA.py -a [groundTruth].nk -b [inferredTree].nk -f [seq_alignment_file].fasta -o [output_file] 
```

### Graph Editing Distance (GED)
```
      $ python compareNewickTrees.py [groundTruth].naive.nk [inferredTree].nk
```

### Examples

```
      # GED - we provide a toy example, real data should be executed in a grid cluster
      $ python compareNewickTrees.py  ../Examples/trees/A.nk ../Examples/trees/B.nk
      
      #COAR
      $ python COAR.py ../Examples/trees/30_1.GT.naive.nk ../Examples/trees/30_1.clonalTree.nk ../Examples/input/30_1.clonalTree.nk
      
      #MRCA
      $ python MRCA.py ../Examples/trees/30_1.GT.naive.nk ../Examples/trees/30_1.clonalTree.nk ../Examples/input/30_1.clonalTree.nk
```
