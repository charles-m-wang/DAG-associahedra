import sys
import itertools
import math
from sage.matroids.constructor import Matroid

def connected(n):
    f=open("new-dependences-"+str(n),'w')
    ground_set=[str(i) for i in xrange(1,n+1)]
    matroid_ground=reduce(lambda x,y:x+y,ground_set)
    power_set=[findsubsets(ground_set,i) for i in xrange(1,n+1)]
    for rank in xrange(1,n+1):
        #power_set elts are lists, power_set[i]=rank i+1 bases
        for num_bases in xrange(1,len(power_set[rank-1])+1):
            #list all possible bases
            possible_bases=findsubsets(power_set[rank-1],num_bases)
            for bases in possible_bases:
                bases_list=[]
                #prase bases into form for Matroid constructor
                for base in bases:
                        bases_list.append(reduce(lambda x,y:x+y,base))
                Mat=Matroid(groundset=matroid_ground,bases=bases_list)
                if Mat.is_valid() and Mat.is_connected():
                    dependence_stmts=dependences(Mat,n)
                    #f.write(matroid_ground+' ')
		    output='{{' # start two lists, inner list is list of bases
                    for base in bases_list:
			base_string='{' # start list of elements in base
		    	for element in list(base):
                            base_string+=element+',' # add elt and sep comma
			base_string=base_string[:-1]+'}'+',' #remove trailing comma, close elt list, add separator
                    	output+=base_string
		    output=output[:-1]+'},{' # remove trailing comma, close inner list of bases, add separator, begin list of dependences
		    dependence='{'
                    for (i,j,K) in dependence_stmts:
                        if len(K)>0:
                            K=reduce(lambda x,y:x+y,K)
                        else:
                            K=''
			condition='{ '
			for k in list(K):
			    condition+=k+','
			condition=condition[:-1]+'}'
                        dependence+="{"+i+", "+j+"},"+condition+"},{"
                    dependence=dependence[:-2]+'}'
		    output+=dependence+'}'
                    f.write(output)
		    f.write('\n')
                    f.flush()
    f.close()
    
def dependences(M, n):
    dependences=[]
    ground=M.groundset()
    r=M.full_rank()
    
    temp=list(ground)
    for i in xrange(1,len(temp)+1):
        for j in xrange(i+1,len(temp)+1):
            x=set(temp)
            x.remove(str(i))
            x.remove(str(j))
            power_set=[findsubsets(x,a) for a in xrange(len(x)+1)]
            power_set=list(itertools.chain(*power_set))
            for K in power_set:
                if M.rank(list(K))+1 == M.rank(list(K+(str(i),))) and M.rank(list(K+(str(i),)))==M.rank(list(K+(str(j),))) and M.rank(list(K+(str(i),)))==M.rank(list(K+(str(i),str(j)))):
                    dependences.append((str(i),str(j),K))
    return dependences
    
        
def findsubsets(S,m):
    #return all m element subsets of S
    return list(itertools.combinations(S,m))
    
def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)
    
connected(int(sys.argv[1]))
