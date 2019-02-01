#!/usr/bin/python

import sys
import getopt
import os

def main():
	params = parseArgs()
	kmers = list()
	if params.k and params.s:
		kmers=get_kmers(params.k, params.s)
	elif params.l:
		kmers=read_kmers(params.l)
	else:
		params.display_help("No input provided.")
	
	#calculate directed overlap graph
	graph = KmerOverlapGraph()
	for mer in kmers:
		graph.addNode(mer)
	graph.setEdges()
	
	deBrujin = deBrujinGraph(graph)


#de Brujin graph implementation
class deBrujinGraph:
	def __init__(self, overlaps):
		self.edgeList = dict()
		self.numEdges = 0
		self.construct(overlaps)
		
	def construct(self, overlaps):
		for node in overlaps:
			if node.hasNeighbors():
				newEdge=Edge(node.getId())
				self.edgeList[node.getId()] = newEdge
				for e in node.getNeighbors():
					if e.getId() not in self.edgeList.keys():
						rightNeighbor = Edge(e.getId())
						self.edgeList[e.getId()] = rightNeighbor
						rightNeighbor.addLeft(newEdge)
					newEdge.addRight(rightNeighbor)

class Edge:
	def __init__(self, key):
		self.key = key
		self.leftNeighbors = dict()
		self.rightNeighbors = dict()
	
	def addLeft(self, other):
		self.leftNeighbors[other.getId()] = other
		
	def addRight(self, other):
		self.rightNeighbors[other.getId()] = other

#Directed graph of kmers (as nodes), with edges as k-1 suffix overlaps
class KmerOverlapGraph:
	def __init__(self):
		self.nodeList = dict()
		self.numNodes = 0
		
	def addNode(self, key):
		if key in self.nodeList.keys():
			return(self.nodeList[key])
		self.numNodes = self.numNodes + 1
		newNode = Node(key)
		self.nodeList[key] = newNode
		this = newNode.getId()
		#self.setEdges()
		return(newNode)
	
	def setEdges(self):
		#calculate edges for new node
		for node1 in self.nodeList.values():
			for node2 in self.nodeList.values():
				left = node1.getId()
				right = node2.getId()
				#print(left[1:])
				#print(right[:-1])
				if left[1:] == right[:-1]:
					#print("k-1 overlap:",left, right)
					if right not in node1.getNeighborIds():
						self.nodeList[left].addNeighbor(node2)

	def getNode(self, key):
		if key in self.nodeList.keys():
			return(self.nodeList[key])
		else:
			return(None)

	def __contains__(self, name):
		return(key in self.nodeList.keys())
	
	def getNodes(self):
		return(self.nodeList)
	
	def getNodeIds(self):
		return(self.nodeList.keys())
	
	def __iter__(self):
		for n in self.nodeList.values():
			yield(n)

#node class for the KmerOverlapGraph
class Node:
	def __init__(self, key):
		self.id = str(key)
		self.neighbors = {}
		self.numNeighbors=0
	
	def addNeighbor(self, nbr):
		self.numNeighbors += 1
		self.neighbors[nbr.getId()] = nbr #nbr should be a node object
	
	def __str__(self):
		ret = str(self.id) + " -> ["
		com=False
		for x in self.getNeighbors():
			if com == True:
				ret = ret + ","
			ret = ret + str(x)
			com = True

		ret = ret + "]"
		return (ret)
	
	def getNeighbors(self):
		return(self.neighbors.values())
		
	def getNumNeighbors(self):
		return(self.numNeighbors)
	
	def hasNeighbors(self):
		return(self.numNeighbors > 0)
	 
	def getNeighborIds(self):
		return(self.neighbors.keys())

	def __iter__(self):
		for n in self.neighbors.values():
			yield(n)
	
	def getId(self):
		return(self.id)
	
	def getWeight(self, nbr):
		return(self.neighbors[nbr])
	
	def __eq__(self, other):
		return(self.id == other.getId())

#reads kmers from txt file, one on each line 
def read_kmers(l):
	if os.path.exists(l):
		with open(l, 'r') as fh:
			try:
				ret = list()
				for line in fh:
					line = line.strip()
					if not line:
						continue
					ret.append(line)
				return(ret)
			except IOError:
				print("Could not read file ",l)
				sys.exit(1)
			finally:
				fh.close()
	else:
		raise FileNotFoundError("File %s not found!"%l)
		
#naive assembly of kmers, assuming they are in order and assuming k-1 overlap
#input is a list of kmers
#output is a contiguous string
def naive_assembly(kmers):
	contig=str()
	for mer in kmers:
		if len(contig) == 0:
			contig = mer
		else:
			contig = contig + mer[-1]
	return(contig)

#function to return k-mers, in order, from a string s
def get_kmers(k, s):
	mers=list()
	keep_going=True
	start=0
	end=start+k
	while end <= len(s):
		mers.append(s[start:end])
		start+=1
		end+=1
	return(mers)

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'k:l:s:h')
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.k=None
		self.s=None
		self.l=None

		#First pass to see if help menu was called
		for o, a in options:
			if o in ("-h"):
				self.display_help("Exiting because help menu was called.")

		#Second pass to set all args.
		for opt, arg_raw in options:
			arg = arg_raw.replace(" ","")
			arg = arg.strip()
			opt = opt.replace("-","")
			#print(opt,arg)
			if opt =="k":
				self.k = int(arg)
			elif opt == "s":
				self.s=str(arg)
			elif opt == "l":
				self.l=str(arg)
			elif opt in ('h', 'help'):
				pass
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if self.k and not self.s:
			self.display_help("Need to specify -s with -k")
		elif self.s and not self.k:
			self.display_help("Need to specify -k with -s")
		elif not self.s and not self.k and not self.l:
			self.display_help("Need to specify either -l, or -k and -s")


	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\ncontigFromReads.py\n")
		print ("Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu")
		print ("\nUsage: ", sys.argv[0], "[-k <kmer> -s <sequence] || [-l list] \n")

		print("""
	Input options:
		-k		: K-mer length 
		-s		: Sequence to sample kmers from 
	- or -
		-l		: .TXT file containing k-mers to build graph from
	Other:
		-h,--help	: Displays help menu""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
