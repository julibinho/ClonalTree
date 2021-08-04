from Bio import SeqIO
import numpy as np
from optparse import OptionParser
import random
from ete3 import Tree
import sys 
sys.path.insert(1, '../clonalTree/')

from BasicTree import *
from BasicSeq  import *

def makeBoolean(var):
	if var == '0':
		return True
	else:
		return False

			
#===================================================================================
def MRCA(ancLeaves1, ancLeaves2, dicSeq):
	sc = 0
	for k,v in ancLeaves1.items():
		seqs = k.split('-')
		key2 = seqs[1]+'-'+seqs[0]
		v2 = ''
		if k in ancLeaves2.keys():
			v2 = ancLeaves2[k]
		elif key2 in ancLeaves2.keys():
			v2 = ancLeaves2[key2]
		if v2 == '':
			print ("ERROR: the trees have diferent nodes")
			print (v, v2)
		   
		elif v != v2:
			print (k, v, v2)
			s = hamming_distance(dicSeq[v], dicSeq[v2])/len(dicSeq[v])
			#print (hamming_distance(dicSeq[v], dicSeq[v2]), s)
			sc = sc + s
	return (sc/len(ancLeaves1))

#===================================================================================
#						Main
#===================================================================================
def main():
	usage = usage = "python MRCA.py -a <nkTree1> -R <R> -b <nkTree2> -L <L> -f <fasta>-o <outputFile> \n"
	parser = OptionParser(usage)
	parser.add_option("-a", "--nkTree1", dest="nkTree1", help="nk file for gcTree")
	parser.add_option("-R", "--nk1Rooted", dest="R", help="is it a rooted Tree (0=True, 1=False)")
	parser.add_option("-b", "--nkTree2", dest="nkTree2", help="nk file for clonalTree")
	parser.add_option("-L", "--nk2Rooted", dest="L", help="is it a rooted Tree (0=True, 1=False)")
	parser.add_option("-f", "--fasta", dest="fastaFile", help="fasta file")
	parser.add_option("-o", "--outputFile", dest="outputFile", help="output file")
	
	(options, args) = parser.parse_args()
	if len(sys.argv) < 9:
		parser.error("incorrect number of arguments")
	
	
	nkTree1 = options.nkTree1
	nkTree2 = options.nkTree2
	fastaFile = options.fastaFile
	outputFile = options.outputFile
	aR = options.R
	bR = options.L

	aRooted = makeBoolean(aR)
	bRooted = makeBoolean(bR)

	
	tree1 = readNKTree(nkTree1, aRooted)
	#print (tree1.get_ascii(show_internal=True))

	tree2 = readNKTree(nkTree2, bRooted)
	#print (tree2.get_ascii(show_internal=True))

	labels, dicSeq  = readFasta2(fastaFile)
	
	cpleaves1 = findCommonAncestorLeaves(tree1, labels)
	#print (len(cpleaves1))

	cpleaves2 = findCommonAncestorLeaves(tree2, labels)
	#print (len(cpleaves2))
	
	if len(cpleaves1) != len(cpleaves2):
		print('ERROR Trees with different amount of nodes', len(cpleaves1), len(cpleaves2))

	score = MRCA(cpleaves1, cpleaves2, dicSeq)
	
	print ("MRCA= ", score)


#===================================================================================
if __name__ == "__main__":
	main()

