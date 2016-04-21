class CodeGraph:

	def __init__(self, gnm, numnode):
		self.graphname = gnm
		self.numnodes = numnode
		self.nodes = [None] * numnode
		self.edges = []

	# Method in CodeGraph class starts with 'g'
	# set graph name
	def gsetname(self, gnm):
		self.graphname = gnm
	
	# get graph name
	def ggetname(self):
		return self.graphname

	# set graph number of nodes
	def gsetnumofnodes(self, nds):
		self.numnodes = nds

	# get graph number of nodes
	def ggetnumofnodes(self):
		return self.numnodes

	# set graph node
	def gsetnode(self, nd):
		n = Nodes(nd)
		self.nodes[nd-1] = n

	# get graph all nodes
	def ggetnodes(self):
		return self.nodes

	# get graph node
	def ggetnode(self, nd):
		return self.nodes[nd-1]
	
	# set node adjacent
	def gsetadjnodes(self, nd, adj):
		#import pdb; pdb.set_trace()
		self.nodes[nd-1].setadjnodes(adj)
	
	# get node adjacent
	def ggetadjnodes(self, nd):
		return self.nodes[nd-1].getadjnodes()

	# get edges of graph
	def ggetedges(self):
		return self.edges

	# get colours used in graph
	def ggetgraphcols(self):
		cols = []
		for nd in self.ggetnodes():
			if nd.getcolour() not in cols:
				cols.append(nd.getcolour())
		return cols

	# called from assignrandomcols method
	def methodradcol(self, cols):
		partition = [' '] * len(cols) 
		ct = 0
		for node in self.ggetnodes():
			if ct >= len(cols):
				ct = 0
			import random
			node.resetcolour()
			node.setcolour(cols[ct])

			tmp = partition[(-node.getcolour())-1]
			tmp = tmp + str(node.getnodenm()) + " "
			partition[(-node.getcolour())-1] = tmp
			ct = ct + 1

		return partition

	# assign random colours to nodes 
	# repeat until graph is not coloured with given colours
	def assignrandomcols(self, cols):
		part = self.methodradcol(cols)
		while len(self.ggetgraphcols()) != len(cols):
			part = self.methodradcol(cols)
		return part

	# add adjacency for remaining nodes which only 
	# appeared as adjacent to some node
	def addremadjacent(self, remadj):
		if None in self.nodes:
			var = self.nodes.index(None) + 1
			Flag = True
		else:
			Flag = False

		if Flag:
			self.gsetnode(var)
			for adj in remadj:
				if str(var) in adj:
					for i in adj.split():
						if '-' in i:
							self.gsetadjnodes(var, -int(i))
							break
			self.addremadjacent(remadj)
	
	# fix if node a appears to be adjacent of b 
	# but b does not have a in its adjacency list.
	def correctadjacent(self):
		for i in self.nodes:
			mn = i.vernm
			for j in i.adjver:
				if mn not in self.nodes[j-1].adjver:
					self.gsetadjnodes(j, mn)
	
	# store edges of graph
	def generateedges(self):
		for nd in self.ggetnodes():
			for adj in nd.getadjnodes():
				if nd.getnodenm() > adj:
					tmp = str(adj) + " - " + str(nd.getnodenm())
				elif nd.getnodenm() < adj:
					tmp = str(nd.getnodenm()) + " - " + str(adj)
				# store in node data structure 
				nd.setedges(tmp)
				if tmp not in self.edges:
					# store in graph data structure
					self.edges.append(tmp)

	#  logic from assignment number1 Graphs & Algorithms
	def checkbipartite(self):
		Q = [None] * (self.numnodes + 1)
		X = [None] * (self.numnodes + 1)
		Q[1] = 1
		M = 1
		X[1] = 1
		k = 1
		while k <= M:
			u = Q[k]
			for v in self.ggetadjnodes(u):
				if X[v] == None:
					M = M + 1
					Q[M] = v
					X[v] = -X[u]
				elif X[v] == X[u]:
					return {'chk': False}
			k = k + 1
			if k == self.numnodes:
				break
		return {'chk': True, 'chnum': 2, 'coloured': True, 'sets': X}

	# Check if a graph is complete graph
	def checkcomplete(self):
		n = self.numnodes
		deg = n-1
		edges = (n*(n-1))/2
		alledges = []
		X = [None] * (self.numnodes + 1)

		for u in self.nodes:
			X[u.vernm] = -u.vernm
			if u.getdegree() != deg:
				return {'chk' : False}

		for u in self.nodes:
			newedges = u.getedges()
			alledges = alledges + newedges

		alledges = set(alledges)
		numedges = len(alledges)

		if numedges != edges:
			return {'chk': False}
		else:
			return {'chk': True, 'chnum': deg+1, 'coloured': True, 'sets': X}
	# Check if a graph is ring of odd number of edges	
	def checkoddring(self):
		X = [None] * (self.numnodes + 1)
		deg = 2
		for u in self.nodes:
			if u.getdegree() != deg:
				return {'chk' : False}

	# Get chromatic number of input graph
	def getchromnumber(self):
		# This method returns a chk, chnum, coloured, sets 

		# check if bipartite graph
		val1 = self.checkbipartite()
		if val1['chk']: return val1

		# check if complete graph
		val2 = self.checkcomplete()
		if val2['chk']: return val2

		# check if a ring with odd numbers 
		val3 = self.checkoddring()
		if val3['chk']: return val3

		chnumb = 0
		for v in self.ggetnodes():
			deg = v.getdegree()
			if deg > chnumb:
				chnumb = deg
		return {'chk': False, 'chnum': chnumb , 'coloured': False, 'sets': None}

# Node class
class Nodes:
	def __init__(self, nm):
		self.vernm = nm
		self.adjver = []
		# check if coloured of not
		self.coloured = False
		# what colour is coloured 
		self.colour = None
		self.sdeg = 0
		self.edges = []

	# get node name
	def getnodenm(self):
		return self.vernm

	# set node adjacent
	def setadjnodes(self, adj):
		self.adjver.append(adj)

	# get node adjacent
	def getadjnodes(self):
		return self.adjver

	# set node edge
	def setedges(self, edge):
		if edge not in self.edges:
			self.edges.append(edge)

	# get node edges
	def getedges(self):
		return self.edges

	# get degree of node
	def getdegree(self):
		return len(self.adjver)

	# get edges connected to node
	def getedges(self):
		edges = []
		u = self.vernm
		for v in self.adjver:
			if u < v:
				edges.append(str(u)+'-'+str(v))
			elif v < u:
				edges.append(str(v)+'-'+str(u))
		return edges

	# check if node is coloured or not.
	def iscoloured(self):
		return self.coloured 

	# set colour of node
	def setcolour(self, val):
		self.coloured = True
		self.colour = val
		return True

	# reset colour of node
	def resetcolour(self):
		self.coloured = False
		self.colour = None

	# get colour of node
	def getcolour(self):
		return self.colour

	# set saturated degree of node
	def setsatdeg(self, val):
		self.sdeg = val 

	# get saturated degree of node
	def satdeg(self):
		return self.sdeg

# Tabus class will store actual restricted moves 
class Tabus:
	def __init__(self, node):
		self.node = node
		self.tabus = []

# TabuList will store objects of TAbus class
class TabuList:
	def __init__(self, limit):
		self.limit = limit
		self.node = []
		#maintain objects of tabus
		self.nodetabu = {}
	
	# add node and it partition as tabu
	def addtabu(self,node, part):
		length = len(self.node)
		if length == self.limit:
			nd = self.node.pop(0)
			tbs = self.nodetabu[nd]
			if len(tbs.tabus) > 1:
				tbs.tabus.pop(0)
			else:
				del self.nodetabu[nd]

		if node in self.nodetabu:
			# No need to create tabus object
			self.node.append(node)
			self.nodetabu[node].tabus.append(part)
		elif node not in self.nodetabu:
			# Create new tabus object
			newnd = Tabus(node)
			self.nodetabu[node] = newnd
			newnd.tabus.append(part)
			self.node.append(node)

	# remove a node from tabu list	
	def removetabu(self, node, part):
		self.node.remove(node)
		tbs = self.nodetabu[node]
		if len(tbs.tabus) > 1:
			tbs.tabus.pop(0)
		else:
			del self.nodetabu[nd]

	# check if tabu, if yes return partition
	def istabu(self, node):
		tmp = node in self.node
		if tmp:
			return True, self.nodetabu[node].tabus
		else:
			return False, []
	
def main():

	# Read the input graph.
	fnm = raw_input("Please enter name of input file: ")
	print "\nYou entered:", fnm
	rd = open(fnm, 'r')
	data = ''
	rdgraph = []
	
	# This code is reading all data from first line as provided
	# in assignment programming question.
	
	while 1:
		line = rd.readline()
		if not line:
			break
		else: 
			data = line

	rd.close()
	gfnm = '' 	# graph name
	tmp = 0		# tmp var
	ct = 1		# counter
	nodes = 0	# number of vertices
	remainingadj = []
	val = ''
	flag = 0
	
	for detail in data.split():
		if ct == 1:
			gfnm = detail
			ct = ct + 1
		elif ct == 2:
			graph = CodeGraph(gfnm, int(detail)) 
			ct = ct + 1
		elif ct == 3:
			if int(detail) == 0:
				remainingadj.append(val)
				break
			if int(detail) < 0:
				if flag == 1:
					remainingadj.append(val)
					val = ''
	
				graph.gsetnode(-int(detail))
				# tmp is used to store u and connect all v to it.
				tmp = -int(detail)
				flag = 1
			else:
				graph.gsetadjnodes(tmp, int(detail))
				val = val + " " + str(-tmp) + " " + detail
	
	# correction method
	graph.addremadjacent(remainingadj)

	# correction method
	graph.correctadjacent()
	
	# populate edges of graph
	graph.generateedges()

	# print node and its adjacency	
	print "\n************ Node and its adjacency *************\n"
	for i in graph.ggetnodes():
		print "\t", i.vernm, "\t",i.adjver
	print "\n*************************************************\n"

	
	algo = raw_input("Which algorithm to use? \n\n\t1. Sequential\n\t2. DSatur\n\t3. TabuCol\n\t4. TabuCol Multi-Threaded\n\t5. DFS\n\t6. BFS\nNote: If a graph is Bipartite, Complete or Ring, none of the above algo's will be used to return the result.\n\nChoice: ")
	print "\nYou entered:", algo
	algo = int(algo)

	algolist = [1,2,3,4,5,6]
	if algo not in algolist:
		algo = algolist[0]
		print "\nUsing default [DSatur Sequential] as the option doesn't exist!\n"

	# Methods used in algorithms **********************

	# Sorting nodes w.r.t degrees 
	def HSort(sequence):

		# Start Heap Sort
		def createheap(sequence):
			start = (len(sequence) - 2) / 2
			while start >= 0:
				arrangeelement(sequence, start, len(sequence) - 1)
				start -= 1
			
		# Shift element downward
		def arrangeelement(sequence, start, end):
			root = start
			while root * 2 + 1 <= end:
				child = root * 2 + 1
				if child + 1 <= end and sequence[child].getdegree() < sequence[child+1].getdegree():
					child += 1

				if child <= end and sequence[root].getdegree() < sequence[child].getdegree():
					tmp = sequence[root]
					sequence[root] = sequence[child]
					sequence[child] = tmp
					root = child
				else:
					return
		createheap(sequence)
		end = len(sequence) - 1
		while end > 0:
			tmp = sequence[end]
			sequence[end] = sequence[0]
			sequence[0] = tmp

			arrangeelement(sequence, 0, end - 1)
			end -= 1

		# Return sorted nodes
		return sequence

	# Sorting nodes w.r.t degrees 
	def HSortSDeg(sequence):

		# Start Heap Sort
		def createheap(sequence):
			start = (len(sequence) - 2) / 2

			while start >= 0:
				arrangeelement(sequence, start, len(sequence) - 1)
				start -= 1
			
		# Shift element downward
		def arrangeelement(sequence, start, end):
			root = start

			while root * 2 + 1 <= end:
				child = root * 2 + 1

				if child + 1 <= end and sequence[child].satdeg() < sequence[child + 1].satdeg():
					child += 1
				elif child +1 <= end and sequence[child].getdegree() < sequence[child + 1].getdegree():
					child += 1

				if child <= end and sequence[root].satdeg() < sequence[child].satdeg():
					tmp = sequence[root]
					sequence[root] = sequence[child]
					sequence[child] = tmp
					root = child
				elif child <= end and sequence[root].getdegree() < sequence[child].getdegree():
					tmp = sequence[root]
					sequence[root] = sequence[child]
					sequence[child] = tmp
					root = child
				else:
					return
		createheap(sequence)
		end = len(sequence) - 1

		while end > 0:
			tmp = sequence[end]
			sequence[end] = sequence[0]
			sequence[0] = tmp

			arrangeelement(sequence, 0, end - 1)
			end -= 1

		# Return sorted nodes
		return sequence

	# Colour the vertex with lowest available colour
	def sdcolourvertices(node, colours, usedcol, sdgraph):
		excludecol = []
		adjnds = []
		breakout = False

		for adj in node.getadjnodes():
			nd = sdgraph.ggetnode(int(adj))
			adjnds.append(nd)
			if nd.iscoloured() == True:
				excludecol.append(nd.getcolour())
		colournow = list(set(colours).difference(excludecol))
		colournow.sort(reverse=True)
		# [-1, -2, -3, -4, -5, -6, -7,.......,-cromnum]

		if len(colournow) == 0:
			print "No colour left to colour the node!"
			# Colouring failed!
			breakout = True
		else:
			# Assign colour to node
			done = node.setcolour(colournow[0])
			if done:
				# Increment saturation degree of adj nodes by one.
				for ad in adjnds:
					ad.setsatdeg(ad.satdeg() + 1)

			if colournow[0] not in usedcol:
				usedcol.append(colournow[0])

		return breakout, usedcol

	def colourlow(node, colours, sdgraph):
		excludecol = []
		adjnds = []
		retcol = None

		for adj in node.getadjnodes():
			nd = sdgraph.ggetnode(int(adj))
			adjnds.append(nd)
			if nd.iscoloured() == True:
				excludecol.append(nd.getcolour())

		colournow = list(set(colours).difference(excludecol))
		colournow.sort(reverse=True)

		if len(colournow) != 0:
			retcol = colournow[0]

		return retcol


	# Calculate the funtion value based on partition
	# Count 1 when an edge has both end points in same
	# partition, 0 otherwise.
	def calcfunval(G, partition):
		count = 0
		for edge in G.ggetedges():
			tmp = edge
			endpts = tmp.split()
			endpts.remove('-')
			for eachpart in partition:
				check1 = " "+str(endpts[0])+" " in eachpart
				check2 = " "+str(endpts[1])+" " in eachpart
				if check1 and check2:
					count = count + 1
		return count


	# check if a given edge has same colour endpoints 
	# if the ends points are of same colour return index
	# of partition [colour] and the endpoints.
	def samecolends(G, partition):
		stuckedges = []

		for edge in G.ggetedges():
			tmp = edge
			endpts = tmp.split()
			endpts.remove('-')
			for eachpart in partition:
				check1 = " "+str(endpts[0])+" " in eachpart
				check2 = " "+str(endpts[1])+" " in eachpart
				if check1 and check2:
					stuckedges.append(edge)
		if len(stuckedges) > 0:
			import random; tmp = random.choice(stuckedges)
			endpts = tmp.split()
			endpts.remove('-')
			for eachpart in partition:
				check1 = " "+str(endpts[0])+" " in eachpart
				check2 = " "+str(endpts[1])+" " in eachpart
				if check1 and check2:
					indx = partition.index(eachpart)
			return endpts, indx 
		else:
			return False, None
	
	# get a random number in range
	def getrandomnum(stop):
		import random; randnum = random.randrange(stop)
		return randnum

	# colour the other endpoint of a same coloured endpoint
	# edge with some other colour [move it to another partition]
	def makemovemaketabu(partition, partindx, y, tbu, funval, graphtabu):
		graphtabu.ggetnode(int(y.strip())).setcolour(-(partindx+1))
		partition[partindx] = partition[partindx].replace(y, " ")

		# get the value of move other than tabus	
		move = partindx
		val, restrictmoves = tbu.istabu(int(y.strip()))
		Flag = True

		# get the restricted moves from tabu list 
		restrictmoves.append(partindx)

		if val:
			while Flag:
				move = getrandomnum(len(partition))
				if move not in restrictmoves:
					Flag = False
		else:
			while Flag:
				move = getrandomnum(len(partition))
				if move != partindx:
					Flag = False
		part = partition[:]
		arr = []

		for i in range(len(part)):
			if i in restrictmoves:
				continue
			arr.append(-(i+1))
		tmp = colourlow(graphtabu.ggetnode(int(y.strip())), arr, graphtabu)

		if tmp != None:
			move = -(int(tmp)+1)

		partition[int(move)] = partition[int(move)] + y.strip()+" "

		return partition

	# Tabucolour algorithm runs till 1500 iteration max
	def tabucol(graph, numofcols):
		import copy
		tabugraph = copy.deepcopy(graph)
		alltabucolours = []
		tabulimit = 7
		tabu = TabuList(tabulimit)

		for i in range(numofcols):
			alltabucolours.append(-i-1)

		s = tabugraph.assignrandomcols(alltabucolours)
		fval = calcfunval(tabugraph, s)
		stmp = s
		count = 0
		plot = ""
		plot1 = ""

		while fval != 0 and count < 500:
			xy, partindx = samecolends(tabugraph, stmp)
			if len(xy) != 0:
				snew = makemovemaketabu(stmp[:], partindx, " "+str(xy[1])+" ", tabu, fval, tabugraph)
				snewfval = calcfunval(tabugraph,snew)
				plot = plot + "\t" + str(snewfval)

				if snewfval < fval:
					plot1 = plot1 + "\t" + str(snewfval)

					# insert y in tabulist for partindx
					tabu.addtabu(int(xy[1].strip()), partindx)
					
					# colour the graph with new colouring
					for each in snew:
						for nde in each.strip().split():
							ndcol = tabugraph.ggetnode(int(nde))
							ndcol.setcolour(-(snew.index(each)+1))

					# make fval and stmp as new vals
					fval = snewfval
					stmp = snew[:]
			else:
				print "No endpoints have same colours!"
				break
			count = count + 1
		
		result = None

		if fval !=0 and count == 500:
			result = False
		elif fval == 0:
			result = True
		else:
			result = False

		return tabugraph, count, result

	# Tabucolour algorithm method which is called by 
	# multiple threads at same time
	def tabumultithreaded(graph, numcols):
		tabugraph, count, result = tabucol(graph, numcols)

		if result:
			return "Colouring found for given graph in "+str(count)+" iterations and "+str(numcols)+" colours."
		else:
			return "\nCan't find colouring for given graph in "+str(count)+" iterations and "+str(numcols)+" colours.\n"

	# get the adjacent non visited vertex of a given node
	def getnonvisitednext(G, V, u):
		nd = G.ggetnode(u)
		nds = nd.getadjnodes()
		ndstake = []

		for i in nds:
			if i not in V:
				ndstake.append(i)

		if len(ndstake) == 0:
			return False
		else:
			import random;
			return random.choice(ndstake)

	# Get the non visited adjacent vertices of a node. 
	def getnonvisitedbfsnds(G, V, u):
		nd = G.ggetnode(u)
		nds = nd.getadjnodes()
		nxtnds = []

		for v in nds:
			if v not in V:
				nxtnds.append(v)

		if len(nxtnds) == 0:
			return False
		else:
			return nxtnds

	# This method is called in TabuCol instead of 
	# randomly colouring the graph at beginning.
	def colourwithsequential(tmpgraph, chromnum):
		allcolours = []
		usedcol = []

		for i in range(chromnum):
			allcolours.append(-i-1)

		for j in tmpgraph.ggetnodes():
			breakout, usedcol = sdcolourvertices(j, allcolours, usedcol, tmpgraph)
			if breakout:
				break

		return len(usedcol)

	# End all methods **********************

	chnum = graph.getchromnumber()
	print "\nChromatic number of graph is: ", chnum['chnum']
	print "\n***************************************\n"

	if chnum['coloured']:
		print "This output is not from the selected option as the graph is either bipartite, complete or ring!\n"
		for i in range(len(chnum['sets'])):
			if chnum['sets'][i] != None:
				print "       ",i, "is coloured ", chnum['sets'][i]
	else:
		chromnum = chnum['chnum']

		# Using Sequential Algorithm
		if algo == 1: 

			print "\nUsing Sequential algorithm!\n"
			import time
			starttm = time.time()

			# Copy graph variable to method variable
			# Make a deepcopy not a reference copy
			import copy
			seqgraph = copy.deepcopy(graph)

			allcolours = []
			usedcol = []

			for i in range(chromnum):
				allcolours.append(-i-1)

			# Store nodes in increasing order of their degree
			seq = [] 
			flag = 1
			for i in seqgraph.ggetnodes():
				seq.append(i)

			sortedseq = seq[:]

			# Start colouring from the vertex with highest degree
			for j in sortedseq:
				breakout, usedcol = sdcolourvertices(j, allcolours, usedcol, seqgraph)
				if breakout:
					break

			print "\nUsed colours were: ", usedcol
			print "\nTime of execution is: ", time.time() - starttm
			print "\nNumber of colours used: ", len(usedcol)

		# Using DSatur algorithm
		elif algo ==2:
			print "\nUsing DSatur algorithm!\n"

			import time
			starttm = time.time()

			# Copy graph variable to method variable
			# Make a deepcopy not a reference copy
			import copy
			sdgraph = copy.deepcopy(graph)

			allcolours = []
			usedcol = []

			for i in range(chromnum):
				allcolours.append(-i-1)

			# Store nodes in increasing order of their degree
			seq = [] 
			flag = 1

			for i in sdgraph.ggetnodes():
				seq.append(i)

			sortedseq = seq[:]

			while len(sortedseq) != 0:
				# Sort nodes according to their satdeg or degree
				sortedseq = HSortSDeg(sortedseq[:])
				sortedseq.reverse()
				breakout, usedcol = sdcolourvertices(sortedseq[0], allcolours, usedcol, sdgraph)
				sortedseq.remove(sortedseq[0])
				if breakout:
					break

			print "\nUsed colours were: ", usedcol
			print "\nTime of execution is: ", time.time() - starttm
			print "\nNumber of colours used: ", len(usedcol)

		# Using TabuCol Algorithm
		elif algo == 3: 			

			print "\nUsing TabuCol algorithm!\n"

			import copy
			tmpgraph = copy.deepcopy(graph)
			numcols = colourwithsequential(tmpgraph, chromnum) 

			import time
			starttm = time.time()

			print "colours are: ", numcols + int(round(((chromnum/100.0)*10)-numcols/8))
			numcols = numcols + int(round(((chromnum/100.0)*10)-numcols/8))

			tabugraph, count, result = tabucol(graph, numcols)
			
			while not result:
				print "\nCan't find colouring for given graph in", count, "iterations and", numcols, "colours.\n"
				numcols = numcols + 1
				tabugraph, count, result = tabucol(graph, numcols)
				if result:
					break
				starttm = time.time()
			if result:
				print "Colouring found for given graph in", count, "iterations and", numcols, "colours."
			else:
				print "\nCan't find colouring for given graph in", count, "iterations and", numcols, "colours.\n"

			print "Time of execution is: ", time.time() - starttm

		# Using TabuCol Multi-Threaded Algorithm
		elif algo == 4: 

			print "\nUsing TabuCol Multi-threaded algorithm!\n"

			import copy
			tmpgraph = copy.deepcopy(graph)
			numcols = colourwithsequential(tmpgraph, chromnum)

			import time
			from multiprocessing.pool import ThreadPool
			starttm = time.time()

			p1 = ThreadPool(processes=4)

			r1 = p1.apply_async(tabumultithreaded, (graph, numcols,))
			r2 = p1.apply_async(tabumultithreaded, (graph, numcols + 2,))
			r3 = p1.apply_async(tabumultithreaded, (graph, numcols + int(round((chromnum/100.0)*10)),))
			r4 = p1.apply_async(tabumultithreaded, (graph, numcols + int(round((chromnum/100.0)*10)) + 1,))

			print r1.get()
			print r2.get()
			print r3.get()
			print r4.get()

			print "Time of execution is: ", time.time() - starttm

		# Using DFS to colour graph
		elif algo == 5:

			print "Using DFS to colour the graph!\n"
			import time
			starttm = time.time()

			import copy
			dfsgraph= copy.deepcopy(graph)

			visited = []
			stack = []
			#u = raw_input("Please enter the value of start node: ")
			#print "\nYou entered:", u, '\n'
			#u = int(u)
			u = 1
			stack.append(u)
			visited.append(u)
			allcolours = []
			usedcol = []

			for i in range(chromnum):
				allcolours.append(-i-1)

			while len(stack) != 0:
				breakout = False
				vnext = getnonvisitednext(dfsgraph, visited, u)	

				if not vnext:
					# When there is no node to proceed
					colnd = stack.pop()
					#print "pop node: ", colnd, "\n\n"
					sendnd = dfsgraph.ggetnode(colnd)
					breakout, usedcol = sdcolourvertices(sendnd, allcolours, usedcol, dfsgraph)
					if len(stack) - 1 == -1:
						breakout = True
					else:
						u = stack[len(stack)-1]
				else:
					# When there is next node to proceed
					stack.append(vnext)
					visited.append(vnext)
					u = vnext

				if breakout:
					break

			print "Time of execution is: ", time.time() - starttm
			print "\nUsed colours: ", usedcol, "\n"
			print "Number of colours used: ", len(usedcol)

		# Using BFS to colour the graph
		elif algo == 6:

			print "Using BFS to colour the graph!\n"
			import time
			starttm = time.time()

			import copy
			bfsgraph = copy.deepcopy(graph)

			visited = []
			queue = []
			#u = raw_input("Please enter the value of start node: ")
			#print "\nYou entered:", u, '\n'
			#u = int(u)
			u = 1
			queue.append(u)
			visited.append(u)
			allcolours = []
			usedcol = []
			
			for i in range(chromnum):
				allcolours.append(-i-1)

			while len(queue) != 0:
				u = queue[0]
				nxtnds = getnonvisitedbfsnds(bfsgraph, visited, u)

				if not nxtnds:
					# When there is no node to proceed
					queue.remove(u)
					sendnd = bfsgraph.ggetnode(u)
					breakout, usedcol = sdcolourvertices(sendnd, allcolours, usedcol, bfsgraph)
				else:
					# When there is next node to proceed
					for ech in nxtnds:
						queue.append(ech)
						visited.append(ech)
					queue.remove(u)
					sendnd = bfsgraph.ggetnode(u)
					breakout, usedcol = sdcolourvertices(sendnd, allcolours, usedcol, bfsgraph)

				if breakout:
					break

			print "Time of execution is: ", time.time() - starttm
			print "\nUsed colours: ", usedcol, "\n"
			print "Number of colours used: ", len(usedcol)

	print "\n*****************End*******************\n"

if __name__ == '__main__':
	main()
