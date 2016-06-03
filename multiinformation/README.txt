This code takes a directed acyclic graph (DAG) as input and outputs a list of inequalities defining the DAG associahedron obtained from the multiinformation. The input format is a file where each line corresponds to an edge (from i to j) and weight:

    i j weight

Examples are included (4vertex, 5vertex, etc.). The output is a .poly file consisting of inequalities, one per line. The inequalities are of the form 

    c_0 + c_1x_1 + ... + c_nx_n >= 0

and are formatted for loading in polymake for polytope computations.

The python script is m-inf.py. To run the python script, invoke the command: 

    python m-inf.py DAGfile
    
(for example, python m-inf.py 4vertex will use the 4-vertex example and output 4vertex.poly). python 2.7.6 and numpy are required to run the code. If you would like randomly selected weights, the file m-inf-rand-wts.py will do this.

-------------------------------------------------------------------------------

The python code is not very numerically friendly, so for larger graphs the output inequalities are not always correct. Tuning is possible via selection of weights and diagonal matrix in the code to obtain the correct polytope.

--------------------------------------------------------------------------

The file incidence.m2 contains Macaulay2 commands for dealing with the polymake output.  In particular, we can find the incidence matrix of from approximate inequalities and approximate vertices.  Although we have not proven the correctness of this heuristics, it works well for DAG associahedra of DAGs with up 7 nodes.
