"""main_cea_is.py

all-in-one file containing classes and function definitions for cluster exact
approximation (cea) of Ising systems (is) for arbitrary bond configurations and
neighbor-relations

Contains class defintions for Edge and UndirectedMultiGraph. Notation
consistent with 

        "Some Basic Definitions in Graph Theory",
        John W. Essam and Michael E. Fisher,
        Review of Modern Physics, 42 (1970) 271

Implements algorithm (running time O(G.v^3)) computing approximate groundstate
spin configurations using a cluster exact approximation scheme as discussed in
AKH1996

        "Cluster-exact approximation of spin glass groundstates",
        A.K. Hartmann,
        Physica A 224 (1996) 480-488

and HR2001

        "Optimization algorithms in Physics",
        A.K. Hartmann, H. Rieger (Eds.),
        Wiley-VCH

For the computation of the cluster groundstates it uses the max-flow min-cut
theorem, i.e., the capacity of a min-capacity cut equals the flow value of a
max-flow in network. For the actual computation of the min-cut, the preflow
push algorithm (running time of O(N^2 sqrt(M)) for N vertices and M edges).
Therefore, the code requires the Networkx graph lib (see
http://networkx.github.io).


Author : Oliver Melchert 
Date   : 04/24/2015 
"""
import sys
import copy
import random
import networkx as nx


# GRAPH DATA STRUCTURES {{{1

class Edge(object):
        """unordered pair of vertices

        EF1970 1.5: unordered pair of distinct vertices from  
        vertex set V. [i,j] is said to be incident with vertices
        i,j and connects them.
        """
        def __init__(self,i=None,j=None):
                """default constructor for a new instance of class Edge

                \prop i vertex from vertex set V
                \prop j vertex from vertex set V
                """
                self.i = i              # unordered pair of vertices
                self.j = j
                self._id   = None       # integer id for edge bookkeeping

        @property
        def id(self):
                return self._id

        @property
        def isLoop(self):
                return True if self.i==self.j else False

        def otherEnd(self,this):
                if this in [self.i,self.j]:
                        return self.i if this==self.j else self.j
                else: 
                        return None

        def __repr__(self):
                """string representation of edge instance"""
                return "[%d,%d].%s"%(self.i,self.j,str(self.id) if self.id!=None else 'X')


class UndirectedMultiGraph(object):
        '''undirected multi graph G = (V,E) without loops

        Contains class defintions for Edge and UndirectedMultiGraph. Notation
        consistent with 

        "Some Basic Definitions in Graph Theory",
        John W. Essam and Michael E. Fisher,
        Review of Modern Physics, 42 (1970) 271

        EF1970 1.12: undirected graph G = (V,E) is an vertex set V, having at
        least one member, together with an associated edge set E. The term 
        multigraph may be used to emphasize that multiedges are allowed.

        Author : Oliver Melchert
        Date   : Feb. 05, 2015 
        '''
	def __init__(self):
		"""default constructor for a new instance of class myGraph"""
		self._nVertices = 0	# number of vertices 
		self._nEdges    = 0  	# number of edges 
                self._eCtr      = 0     # counter for integer edgeIds
		self._incEdgId  = {}    # unordered set of node:(edgeId,otherEnd) pairs 
                self._E         = {}    # unordered set of edgeId:edge pairs 

        @property	
	def V(self): 
		"""returns vertex set of the graph, see EF1970 1.1"""
		return self._incEdgId.keys()

        @property
        def E(self):
                """returns set of individual edges of the graph, see EF1970 1.7"""
                return self._E.values()

        @property
	def v(self): 
                """returns number of vertices in graph, see EF1970 1.1"""
		return self._nVertices

        @property
	def e(self): 
                """returns number of edges in graph, see EF1970 1.7"""
		return self._nEdges

	def adjVertices(self,i):
                """returns set of vertices j adjacent to vertex j, see EF1970 1.31"""
		return set([otherEnd for (eId,otherEnd) in self._incEdgId[i].iteritems()]) 

	def incEdges(self,i):
                """returns set of edges [i,j] incident to vertex i, see EF1970 1.5"""
		return set([self._E[eId] for (eId,otherEnd) in self._incEdgId[i].iteritems()]) 

	def deg(self,i):
                """degree of a vertex, see EF1970 1.19"""
		return len(self._incEdgId[i])

	def addVertex(self,i):
		if i not in self.V: 
			self._incEdgId[i]={}
			self._nVertices += 1

	def delVertex(self,i):
                for eij in self.incEdges(i):
                        self.delEdge(eij)
                del self._incEdgId[i]
                self._nVertices -= 1

	def addEdge(self,e):
                if not e.isLoop:
                   self.addVertex(e.i)
                   self.addVertex(e.j)
                   e._id = self._eCtr
                   self._eCtr+=1
                   self._E[e._id] = e
                   self._incEdgId[e.j][e._id]=e.i
                   self._incEdgId[e.i][e._id]=e.j
                   self._nEdges += 1

	def delEdge(self,e):
                del self._incEdgId[e.i][e._id]
                del self._incEdgId[e.j][e._id]
                del self._E[e._id]
                self._nEdges -= 1

        def contains(self,e):
                """True if graph contains edge.
                
                Comparison is made based on the terminal vertices, not the edge id
                """
                return True if (e.j in self._incEdgId[e.i].values()) else False

        def __str__(self):
                """return string with graph in lemon graph format"""
                string = '@nodes\nlabel deg\n'
                for i in self.V:
                        string += '%3d %3d\n'%(i,self.deg(i))

                string += '@edges\n     label\n'
                for e in self.E:
                        string += '%3d  %3d  %3d\n'%(e.i,e.j,e.id)

                return string

# 1}}}

# READ CFGs 1{{{

def readGraph_bondList_csg(fName):
        """read bond configuration in "Colone Spin glass server" format 
        """
        G = UndirectedMultiGraph()
        for line in open(fName):
               c = line.split()
               if len(c)==3 and c[0]!="c":
                       eij     = Edge(int(c[0])-1, int(c[1])-1)
                       eij.wgt = int(c[2])
                       G.addEdge(eij)
        return G 

# 1}}}

# OBSERVABLES {{{1

def magnetization(spin):
        """magnetization per spin
        
        \param[in] spin Array containing spin orientation for each vertex  
        """
        return float(sum(spin))/len(spin)


def energy(G,spin):
        """energy per spin

        NOTE:
        -# edges need to have a ".wgt" attribute

        \param[in] G Undirected weighted graph
        \param[in] spin Spin orientation for energy computation
        """
        locErg = [-spin[e.i]*spin[e.j]*e.wgt for e in G.E]
        return float(sum(locErg))/G.v

def energy_f(G,h,spin):
        """energy for given bond, spin, and field setup
        
        NOTE:
        -# edges need to have a ".wgt" attribute

        INPUT:
        \param[in] G Undirected weighted graph containing bond setup
        \param[in] h Local fields at spin locations 
        \param[in] spin Spin configuration 
        
        RETURNS: (E)
        \param[ret] E Energy of spin configuration
        """
        locErg  = [spin[e.i]*spin[e.j]*e.wgt for e in G.E]
        locPiv  = [h[i]*spin[i] for i in range(G.v)]
        return -sum(locErg)-sum(locPiv)


def nBroken(G,spin):
        """determine number of broken bonds

        NOTE:
        -# edges need to have a ".wgt" attribute

        \param[in] G Undirected weighted graph
        \param[in] spin Spin orientation for computation of broken bonds 
        """
        locErg = [-spin[e.i]*spin[e.j]*e.wgt for e in G.E]
        return len(filter(lambda x: x>0, locErg))

#1}}}

# CLUSTER-EXACT APPROX. {{{1

def labelClusterVertices(G):
        """construct subgraph that exhibits no frustration

        NOTE:
        -# edges need to have a ".wgt" attribute
        -# implemented according to AKH1996:
           "Cluster-exact approximation of spin glass groundstates",
           A. K. Hartmann, 
           Physica A 224 (1996)

        \param[in] G Undirected weighted graph
        \param[out] t Gauge variables that yield ti*tj*Jij >= 0 
        """
        d = {i:0 for i in G.V}  # flag indicating if vertex has been marked 
        t = {i:0 for i in G.V}  # gauge variable to yield Jij*ti*tj >= 0

        # seed vertex for cluster construction is chosen at random
        i0    = random.choice(d.keys())
        S     = [i0]
        d[i0] = 1

        while S:
            i = S.pop(random.randint(0,len(S)-1)) 

            A = [e for e in G.incEdges(i) if t[e.otherEnd(i)]!=0]
            if not A:
                # no neighboring gauge variable set, yet
                t[i] = 1
            else:
                # check if vertex-neighborhood requires unique gauge 
                # variable at vertex i and extend cluster if possible 
                aList = [1 if e.wgt*t[e.otherEnd(i)] > 0 else -1 for e in A]
                if len(set(aList))==1:
                   t[i] = aList[0] 
                else:
                   t[i]=0
            
            # if vertex i was added to cluster, put its still unmarked
            # neighbors on stack for later processing
            if t[i]!=0:
                for j in G.adjVertices(i):
                   if d[j]==0:
                       S.append(j)
                       d[j]=1

        return t


def constructClusterSubgraph(G,s,t):
        """slice cluster subgraph from original graph

        Construct subgraph from orignial graph, comprising only 
        (relabeled) cluster vertices and their corresponding bonds. 
        Interaction of cluster spins with non-cluster spins is
        accounted for by means of local fields attached to vertices
        at cluster surface.

        NOTE:
        -# vertices of Subgraph are relabeled 0...nCluster

        \param[in]  G Undirected weighted graph
        \param[in]  s Spin configuration used to compute local fields
                      summarizing the interaction with non-cluster spins
        \param[in]  t Gauge variables that yield ti*tj*Jij >= 0 
        \param[out] A Subgraph comprising cluster vertices only
        \param[out] h Local fields at vertices summarizing the interaction
                      with non-cluster spins
        \param[out] i2c Dictionary translating vertex labels from
                      subgraph-id to original graph-id
        """
        A   = UndirectedMultiGraph() 
        h   = {}
        C   = set([i for i,ti in t.items() if ti]) 
        c2i = {ci:i for (i,ci) in enumerate(C)}
        i2c = {i:ci for (ci,i) in c2i.items()} 

        for ci in C:
             h_ci = 0
             for e in G.incEdges(ci):
                 cj = e.otherEnd(ci)
                 if cj in C:
                    if ci<cj: 
                       eNew     = Edge(c2i[ci],c2i[cj])
                       eNew.wgt = abs(e.wgt)
                       A.addEdge(eNew)
                 else:
                       h_ci += t[ci]*e.wgt*s[cj]
             h[c2i[ci]] = h_ci 

        return A,h,i2c


def constructEquivalentNetwork(G,h):
        """build equivalent network for computation of max-flow/min-cut

        NOTE:
        -# requires Networkx graph lib (see http://networkx.github.io)
        -# equivalent network constructed according to chapter 6 on
           "Maximum Flow networks" of 
           "Optimization Algorithms in Physics",
           A.K. Hartmann, H. Rieger

        \param[in] G Undirected weighted graph
        \param[in] h Local fields at vertices 
        \param[out] EN Directed equivalent network  
        \param[out] c0N  
        """
        EN = nx.DiGraph()

        # inner vertices, edges and edge-capacities
        for e in G.E:
           EN.add_edge(e.i+1,e.j+1,capacity=4*e.wgt)

        # outer source/sink-vertices, edges and edge-capacities
        w_sum = 0
        for i in G.V:

            w_ip = -2*h[i]
            for e in G.incEdges(i):
                w_ip += 2*e.wgt if i > e.otherEnd(i) else -2*e.wgt

            (W1,W2) = (0,w_ip) if w_ip > 0 else (-w_ip,0)
            EN.add_edge(0,i+1,capacity = W1)
            EN.add_edge(i+1,G.v+1,capacity = W2)
            w_sum+=- (W1+W2)

        # leave (0,N)-edge out, since it would be crossed by any cut,
        # anyways. After the capacity of the min-cut is found, c0N
        # has to be added in order to obtain the proper energy
        c0N = -sum([e.wgt for e in G.E])-w_sum/2

        return EN,c0N


def minCutVertexPartition(EN):
        """obtain vertex partition corresponding to minimum cut

        Computes the value and the vertex partition of a minimum (0,N)-cut
        in directed input network. Thereby it uses the max-flow min-cut 
        theorem, i.e., the capacity of a min-capacity cut equals the flow 
        value of a max-flow in network.

        NOTE:
        -# requires Networkx graph lib (see http://networkx.github.io)
        -# uses preflow-push algorithm (running time of O(N^2 sqrt(M)) 
           for N vertices and M edges) to compute max-flow 

        \param[in] EN Directed equivalent network  
        \param[out] (S,Sbar) Vertex partition induced by min-cut in EN
        """
        cut_value,(S,Sbar) = nx.minimum_cut(EN, 0, EN.number_of_nodes()-1)
        return (S,Sbar)


def clusterGroundstateCfg(s,t,i2c,(S,Sbar)):
        """spin configuration containing cluster groundstate
        
        compute new spin configuration for full ising system, containing
        the groundstate spin configuration of the non-frustrated cluster

        
        \param[in] s Initial spin configuration
        \param[in] t Gauge variables that yield ti*tj*Jij >= 0 
        \param[in] i2c Dictionary translating vertex labels from
                     subgraph-id to original graph-id
        \param[in] (S,Sbar) Vertex partition induced by min-cut  
        \param[out] sNew Spin configuration including cluster-groundstate 
        """
        sNew      = copy.deepcopy(s)
        upSpins   = [ ( 1, i2c[si-1]) for si in S-{0}]
        downSpins = [ (-1, i2c[si-1]) for si in Sbar-{max(Sbar)}]
        
        for (sVal,v) in upSpins + downSpins:
                sNew[v] = sVal*t[v]

        return sNew


def CEA(G,s):
        """cluster exact approximation of groundstates for Ising system

        algorithm (running time O(G.v^3)) computing approximate groundstate 
        spin configurations for Ising systems, using a cluster exact 
        approximation scheme as discussed in AKH1996

           "Cluster-exact approximation of spin glass groundstates",
           A.K. Hartmann,
           Physica A 224 (1996) 480-488

        and HR2001

           "Optimization algorithms in Physics",
           A.K. Hartmann, H. Rieger (Eds.),
           Wiley-VCH

        NOTE:
        -# Routine cannot increase configurational energy
        -# might get stuck in local minima of the energy landscape, as 
           a remedy, AKH1996 advices to use several independent runs
           per disorder instance and to use the one that ends on the
           lowest energy level (AKH1996 uses 3 runs per instance)
        -# procedure can be improved by combining it with a genetic
           algorithm that maintains a pool of several configurations
           see HR2001
        
        \param[in] G Undirected weighted graph representing Ising system
        \param[in] s Initial spin configuration
        \param[out] sNew Spin configuration including cluster-groundstate 
        """
        # label all vertices that can be attributed to a non-frustrated
        # cluster of spins and their bonds
        t        = labelClusterVertices(G)

        # construct a subgraph of G that contains only the (relabeled)
        # vertices that belong to the cluster and compute a set of 
        # local fields that encode the interaction of the cluster-spins
        # with the non-cluster spins
        A,h,i2c  = constructClusterSubgraph(G,s,t) 

        # construct the equivalent capacitated network that corresponds
        # to the cluster for which the groundstate spin configuration 
        # can be computed exactly using combinatorial optimization methods
        RN,c0N   = constructEquivalentNetwork(A,h)

        # find partition of vertices into up/down spins by using the 
        # max-flow/min-cut theorem and the preflow-push algorithm to comute
        # a minimum cut (there are possibly more than one min-cut)
        (S,Sbar) = minCutVertexPartition(RN)

        # compute new spin configuration for full ising system, containing
        # the groundstate spin configuration of the non-frustrated cluster
        sNew     = clusterGroundstateCfg(s,t,i2c,(S,Sbar))

        return sNew


# 1}}}

# UTILITY FUNCTIONS {{{1

def print_graph(self,s):
        """return string with graph in lemon graph format"""
        string = '@nodes\nlabel deg\n'
        for i in self.V:
                string += '%3d %3d\n'%(i,s[i])

        string += '@edges\n     label\n'
        for e in self.E:
                string += '%3d  %3d  %6d\n'%(e.i,e.j,e.wgt)

        return string


def print_auxGraph(self,h,s,t,i2c):
        """return string with graph in lemon graph format"""
        string = '@nodes\nlabel deg\n'
        string+= '# i hi si ti\n'
        for i in self.V:
                string += '%3d %3d %3d %3d\n'%(i2c[i],h[i],s[i2c[i]],t[i2c[i]])

        string += '@edges\n     label\n'
        for e in self.E:
                string += '%3d  %3d  %6d\n'%(i2c[e.i],i2c[e.j],e.wgt)

        return string


def print_cfg(s):
        return "".join(["+" if si>0 else "-" for si in s])

# }}}1


def main_CEA():

        bondFile="./2d_gauss_L8.cfg"
        G = readGraph_bondList_csg(bondFile)

        s = [random.choice([-1,1]) for i in range(G.v)]
        print "# CLUSTER EXACT APPROXIMATION (CEA) OF SPIN GLASS GROUNDSTATES"
        print "# SIMPLE ITERATIVE APPLICATION OF THE CEA SCHEME TO ISING SYSTEM"
        print "# SPECIFIED IN FILE: ",bondFile
        print "# (iter) (erg) (cfg)"
        print 0,energy(G,s), print_cfg(s) 

        for i in range(0,1000):
                s = CEA(G,s)
                print i+1, energy(G,s), print_cfg(s) 


main_CEA()
# EOF: main_cea_is.py
