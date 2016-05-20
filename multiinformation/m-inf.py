import sys
import itertools
from numpy import matlib,multiply,diagflat,array,newaxis,copy,diag,sqrt,log
from numpy.linalg import solve, det

def getCorrelationMatrix(DAG,default=False, default_weight=1):

    #compute correlation matrix from DAG
    f=open(str(DAG),'r')
    edges=[float(i) for i in f.read().strip().split()] # read in edges (and weights)
    f.close()

    num = 3 #using specified weights
    if default == "True":
        num = 2 #using default weight

    for i in xrange(len(edges)/num):
        edges[num*i]=int(edges[num*i])
        edges[num*i+1]=int(edges[num*i+1])

    temp = [[edges[num*i],edges[num*i+1]] for i in xrange(len(edges)/num)] #get a list of all the vertices
    temp = list(itertools.chain(*temp))
    n = reduce(max,temp) #compute number of vertices
    checkConnected(temp, n) #sanity check:make sure every vertex is used

    adj_list = [[edges[num*i+j] for j in xrange(num)] for i in xrange(len(edges)/num)] #create adjacency list
    if default == "True":
        for i in xrange(len(adj_list)):
            adj_list[i].append(default_weight) #add weights to adjacency list

    I = matlib.identity(n)
    A = matlib.zeros((n,n))
    for edge in adj_list:
        A[edge[0]-1,edge[1]-1]=edge[2] #create adjacency matrix from adjacency list using weights
    
    return (I-A)*(I-A).getT(),n #(I-A)D(I-A)^T, here D=I, but in general can choose D differently
    #product_diagonal=1./sqrt(diag(product)) #get diagonal matrix to normalize diagonal of (I-A)(I-A)^T
    #return diagflat(product_diagonal)*product*diagflat(product_diagonal),n

def logDets(DAG, default=False):
    #compute multiinformation for all subsets

    f=open(DAG+".poly",'w') # file for writing inequalities
    f.write("INEQUALITIES\n")
    (K,n)=getCorrelationMatrix(DAG,default) # (I-A)(I-A)^T normalized so the diagonal is all 1s
    
    vertices = [i for i in xrange(n)]
    subsets = [findsubsets(vertices,i) for i in xrange(1,n+1)] # list all subsets of vertices
    subsets = list(itertools.chain(*subsets))

    for subset in subsets:
        vertices=set([i for i in xrange(n)])
        vertices.difference_update(subset)
        vertices=sorted(list(vertices)) #find which rows/columns to remove

        temp=copy(K) # get principal minors
        for x in reversed(vertices):
            temp=minor(temp,x,x)

        vertices=sorted(list(subset))
        vector = [0]*(n+1)
        vector[0]=log(det(temp))
        
        for i in list(subset): #find which of the x_i participate in this minor
            vector[i+1]=-1
        for i in vector:
            f.write(str(i)+" ")
        f.write('\n')
        
        if len(list(subset))==n: #sum_{i=1}^n x_i=m([n])
            for i in vector:
                f.write(str(-i)+" ")
            f.write('\n')
        f.flush()
    f.close()

def minor(arr,i,j):
    # ith row, jth column removed
    return arr[array(range(i)+range(i+1,arr.shape[0]))[:,newaxis],
               array(range(j)+range(j+1,arr.shape[1]))]

def findsubsets(S,m):
    # enumerate all the size m subsets of S
    return list(itertools.combinations(S,m))

def checkConnected(E,V):
    #check whether every vertex appears in some edge
    all_vertices = set([i for i in xrange(1,V+1)])
    all_vertices.difference_update(set(E))
    if all_vertices:
        print "Graph not connected: vertices (" + str([i for i in all_vertices]) + ") unused."

if len(sys.argv)==2:
    logDets(sys.argv[1])
else:
    logDets(sys.argv[1],sys.argv[2])

#print getCorrelationMatrix(sys.argv[1],sys.argv[2])
