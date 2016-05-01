import sys
import itertools
from itertools import combinations as comb
import math
from sage.matroids.constructor import Matroid

#n is num elts, s is size of subsets, generate in reverse lexicographic order
def revLexBases(n,s):
    return [temp[::-1] for temp in sorted([sub[::-1] for sub in list(comb([str(x) for x in xrange(1,n+1)],s))])]
    
def complement(subset, whole_set):
    return sorted(list(set(whole_set).difference(set(subset))))

def connected_isom(n):
    f=open("connected-dependences-"+str(n),'w')
    ranks = [open("n"+str(n)+"rank"+str(i),'r').read().split() for i in xrange(1,n/2+1)]
    
    bases_lists = [revLexBases(n,i) for i in xrange(1,n/2+1)]

    for i in xrange(n/2):
        for matroid_rep in ranks[i]:
            matroid = []
            dual_matroid=[]
            for j in xrange(len(matroid_rep)):
                if matroid_rep[j] == '*':
                    matroid.append(reduce(lambda x,y:x+y,bases_lists[i][j]))
                    dual_matroid.append(reduce(lambda x,y:x+y, complement(bases_lists[i][j],tuple([str(z) for z in xrange(1,n+1)]))))
            #matroid has rank i
            #dual matroid has rank 7-i
            dual_matroid = dual_matroid[::-1]
            Mi = Matroid(groundset = reduce(lambda x,y:x+y,[str(z) for z in xrange(1,n+1)]), bases = matroid)
            M_nminusi = Matroid(groundset = reduce(lambda x,y:x+y,[str(z) for z in xrange(1,n+1)]), bases = dual_matroid)
            if Mi.is_connected():
                Mi_dep = dependences(Mi)
                f.write(parse_dependences(matroid,Mi_dep))
                f.write('\n')
            if M_nminusi.is_connected():
                M_nminusi_dep = dependences(M_nminusi)
                f.write(parse_dependences(dual_matroid,M_nminusi_dep))
                f.write('\n')
    f.flush()
    f.close()
    
def parse_dependences(bases_list, dep_list):
    output='{ { ' # start two lists, inner list is list of bases
    for base in bases_list:
        output += "{"+reduce(lambda x,y:x+','+y,list(base),'')[1:]+"},"
    output = output[:-1]+" }, { "
    for (i,j,K) in dep_list:
        output += "{{"+str(i)+","+str(j)+"},{"+reduce(lambda x,y:x+','+y,sorted(K),'')[1:]+"}},"
    output = output[:-1]+" } }"
    return output

    
def dependences(M):
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
    return list(comb(S,m))
    
connected_isom(int(sys.argv[-1]))
