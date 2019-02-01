#!/usr/bin/python

import sys

def main():
	if not sys.argv[1] and sys.argv[2]:
		sys.exit(1)
		print("Usage: readsFromContigs.py <k> <string>")

	k=int(sys.argv[1])
	s=str(sys.argv[2])

	kmers=get_kmers(k, s)
	
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
		
		#calculate edges for new node
		print(key, "checking for overlaps")
		for node in self.nodeList.values():
			other = node.getId()
			if other != this:
				print(this[-(self.k-1):])
				print(other[:(self.k-1)])
				if this[-(self.k-1):] == other[:(self.k-1)]:
					print("k-1 overlap:",key, other)
					self.nodeList[this].addNeighbor[node]
		return(newNode)

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
	
	def addNeighbor(self, nbr, weight=0):
		self.neighbors[nbr] = weight #nbr should be a node object
	
	def __str__(self):
		return str(self.id) + " -> " + str([x.id + "," for x in self.neighbors])
	
	def getNeighbors(self):
		return(self.neighbors)
	
	def getId(self):
		return(self.id)
	
	def getWeight(self, nbr):
		return(self.neighbors[nbr])
	
	def __hash(self):
		return hash(self.id)
	
	def __eq__(self, other):
		return(self.id == other.getId())
	
	def __ne__(self, other):
		return not(self == other)
		
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
