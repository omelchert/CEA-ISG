main_cea_is.py

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
