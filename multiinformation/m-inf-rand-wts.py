import sys
import os
import itertools
from numpy import matlib,multiply,diagflat,array,newaxis,copy,diag,sqrt,log
from numpy.linalg import solve, det
import random

def getCorrelationMatrix(DAG):
    f=open(str(DAG),'r')
    edges=[int(i) for i in f.read().strip().split()] # read edges (and weights)
    f.close()
    temp = [[edges[2*i],edges[2*i+1]] for i in xrange(len(edges)/2)]
    temp = list(itertools.chain(*temp))
    n = reduce(max,temp) # number of vertices
    checkConnected(temp, n)

    adj_list = [[edges[2*i+j] for j in xrange(2)] for i in xrange(len(edges)/2)]
    for i in xrange(len(adj_list)):
        adj_list[i].append(random.randrange(25,101)/100.0)
    print adj_list
    I = matlib.identity(n)
    A = matlib.zeros((n,n))
    for edge in adj_list:
        A[edge[0]-1,edge[1]-1]=edge[2]
    product=(I-A)*(I-A).getT()
    return product,n

def checkConnected(E,V):
    all_vertices = set([i for i in xrange(1,V+1)])
    all_vertices.difference_update(set(E))
    if all_vertices:
        print "Graph not connected: vertices (" + str([i for i in all_vertices]) + ") unused."

def logDets(DAG):
    f=open(DAG+".poly",'w')
    f.write("INEQUALITIES\n")
    (K,n)=getCorrelationMatrix(DAG)
    norm_diag=diagflat(diag(K))
    vertices = [i for i in xrange(n)]
    subsets = [findsubsets(vertices,i) for i in xrange(1,n+1)]
    subsets = list(itertools.chain(*subsets))
    for subset in subsets:
        vertices=set([i for i in xrange(n)])
        vertices.difference_update(subset)
        vertices=sorted(list(vertices))
        temp_K=copy(K)
        temp_norm=copy(norm_diag)
        for x in reversed(vertices):
            temp_K=minor(temp_K,x,x)
            temp_norm=minor(temp_norm,x,x)
        vertices=sorted(list(subset))
        vector = [0]*(n+1)
        vector[0]=(.5)*log((det(temp_K)+0.0)/det(temp_norm))
        for i in list(subset):
            vector[i+1]=-1
        for i in vector:
            f.write(str(i)+" ")
        f.write('\n')
        if len(list(subset))==n:
            for i in vector:
                f.write(str(-i)+" ")
            f.write('\n')
        f.flush()
    f.close()
    os.system("polymake "+DAG+".poly F_VECTOR")

def minor(arr,i,j):
    # ith row, jth column removed
    return arr[array(range(i)+range(i+1,arr.shape[0]))[:,newaxis],
               array(range(j)+range(j+1,arr.shape[1]))]

def findsubsets(S,m):
    return list(itertools.combinations(S,m))

logDets(sys.argv[1])
