"# Vertex-Colouring"

Algorithms implemented are as below:

1.	Sequential Colouring
2.	Degree Saturation Algorithm
3.	Tabu Colouring
4.	Tabu Colouring Multithreading 
5.	DFS Colouring 
6.	BFS Colouring 

Time required for execution order:
Seq, BFS, DFS, DSat, Tabu Mtt, Tabu.

Lowest colouring found order:
DSat, BFS, DFS, Seq, Tabu Mtt, Tabu.



Methodology:

Sequential: Simple 'for' loop. Starts colouing sequentially from 1 to n.

Degree Saturation: There are two type of variables here: Degree and Saturation Degree. Degree is number of neighbours of current node and saturation degree is the number of adjacent coloured nodes of current node. The nodes are sorted accourding to first saturation degree and if two node have same saturation degree the they are sorted according to there degree. The node with highest degree is coloured first. The process repeats itself that saturation degree of nodes keep on changing. 