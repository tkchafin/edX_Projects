#!/usr/bin/python

import sys

def main():
	# if not sys.argv[1] and sys.argv[2]:
	# 	sys.exit(1)
	# 	print("Usage: readsFromContigs.py <k> <string>")

	# k=int(sys.argv[1])
	# s=str(sys.argv[2])

	# kmers=get_kmers(k, s)
	
	
	
	#calculate directed overlap graph
	graph = KmerOverlapGraph(k)
	for mer in kmers:
		print("adding:", mer)
		graph.addNode(mer)
		
	for node in graph:
		print(node)
	
#Directed graph of kmers (as nodes), with edges as k-1 suffix overlaps
class KmerOverlapGraph:
	def __init__(self, k=2):
		self.nodeList = {}
		self.numNodes = 0
		self.k = int(k)
		
	def addNode(self, key):
		self.numNodes = self.numNodes + 1
		newNode = Node(key)
		self.nodeList[key] = newNode
		this = newNode.getId()
		self.setEdges()
		return(newNode)
	
	def setEdges(self):
		#calculate edges for new node
		for node1 in self.nodeList.values():
			for node2 in self.nodeList.values():
				left = node1.getId()
				right = node2.getId()
				if left != right:
					print(left[-(self.k-1):])
					print(right[:(self.k-1)])
					if left[-(self.k-1):] == right[:(self.k-1)]:
						print("k-1 overlap:",left, right)
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
		return(self.nodeList.keys())
	
	def __iter__(self):
		return iter(self.nodeList.values())

#node class for the KmerOverlapGraph
class Node:
	def __init__(self, key):
		self.id = str(key)
		self.neighbors = {}
	
	def addNeighbor(self, nbr):
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
		return(self.neighbors)
	 
	def getNeighborIds(self):
		return(self.neighbors.keys())
	
	def getId(self):
		return(self.id)
	
	def getWeight(self, nbr):
		return(self.neighbors[nbr])
	
	def __eq__(self, other):
		return(self.id == other.getId())

		
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

#Call main function
if __name__ == '__main__':
    main()
