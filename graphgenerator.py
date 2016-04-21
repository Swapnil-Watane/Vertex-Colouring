def GenerateGraph(nds, prob):
	graph = [' '] * nds 
	for i in range(1, nds+1):
		for j in range(i+1, nds+1):
			import random;
			r = random.uniform(0.0, 1.0)
			# Check condition and add edge
			if r <= prob:
				graph[i-1] = graph[i-1] + str(j)+" "

	filename = str(nds) + "nds-" + str(prob) + ".txt"
	graphnm = "Graph" + str(nds)
	# Open File
	target = open(filename, 'w')
	# Clear all
	target.truncate()
	# Write Graph Name
	target.write(graphnm)
	# Write \r like newline
	target.write("\r")
	# Write Number of Nodes
	target.write(str(nds))
	# Write \r like newline
	target.write("\r")

	for line in graph:
		# Write -nodes and its adjacency
		target.write("-" + str(graph.index(line)+1) + line.rstrip())
		# Write \r like newline
		target.write("\r")
	
	# Write end of graph 0
	target.write("0")
	# File close
	target.close()

if __name__ == '__main__':
	# Number of nodes
	nodes = 800
	# Probability
	probability = 0.7
	GenerateGraph(nodes, probability)

