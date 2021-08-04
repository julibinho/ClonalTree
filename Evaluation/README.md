
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
      # GED - we provide a toy example
      $ python compareNewickTrees.py  ../Example/trees/A.nk ../Example/trees/B.nk
```
