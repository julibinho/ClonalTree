
import numpy as np
from optparse import OptionParser
from scipy.spatial import distance
from MSTree import *
from BasicTree import *
from BasicSeq import *
import sys 




#===================================================================================
#						Main
#===================================================================================
def main():
	usage = usage = "python clonalTree.py -i <fastaFile> -r <revision> -o <outputFile> \n"
	parser = OptionParser(usage)
	parser.add_option("-i", "--fastaFile", dest="fastaFile",  help="sequences in fasta format")
	parser.add_option("-o", "--outputFile", dest="outputFile",  help="output file")
	parser.add_option("-a", "--useAbundance", dest="useAbundance",  help="if 0 it uses abundance")
	parser.add_option("-r", "--revision", dest="revision",  help="if 0 it performs revision")
	parser.add_option("-t", "--trim", dest="trim",  help="if 0 it performs trimming tree")
	
	(options, args) = parser.parse_args()
	if len(sys.argv) < 5:
		parser.error("incorrect number of arguments")
	
	fastaFile = options.fastaFile
	outputFile = options.outputFile
	useAbundance = options.useAbundance
	revision = options.revision
	trim =  options.trim

	useAbundance = makeBoolean(useAbundance)
	revision = makeBoolean(revision)
	trim =  makeBoolean(trim)

	#labels, root, arraySeqs, abundance  = readFastaRepeat(fastaFile)

	labels, root, arraySeqs, abundance, dico =  readFastaAbundance(fastaFile)#; print (labels)
	
	adjMatrix = createAdjMatrix(arraySeqs) #;print(adjMatrix)
	
	tree, infoTree = kruskalMST(adjMatrix, root, labels, abundance, useAbundance) #;print (infoTree)
	    
	#print (tree.get_ascii(show_internal=True))

	if trim:
		tree = trimming(tree, labels, adjMatrix) #;print (tree.get_ascii(show_internal=True))
	
	if revision:
		tree = editTree(tree, adjMatrix, labels)
	tree.write(format=1, outfile=outputFile)
	f = open(outputFile +'.csv', 'w')
	f.write(infoTree); f.close()
	print ('done')
	

#===================================================================================
if __name__ == "__main__":
	main()

