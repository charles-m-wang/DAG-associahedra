import itertools
import sys
import math
from sage.matroids.constructor import Matroid

#r is the rank function to test, n is the number of elements
def fn_to_CIs(r,n):
    CIs=set()
    for i in xrange(1,n+1):
        for j in xrange(i+1,n+1):
            K=(reduce(lambda x,y:x+y, [str(x) for x in xrange(1,n+1)]).replace(str(i),'')).replace(str(j),'')
            subsets=[findsubsets(K,x) for x in xrange(0,n-1)]
            for size in subsets:
                for subset in size:
                    K=reduce(lambda x,y:x+y,subset,'')
                    if r(K+str(i))+r(K+str(j))==r(K+str(i)+str(j))+r(K):
                        CIs.add((i,j,K))			  
    return CIs

def checkGaussoids(n):
    num_gaussoid=0
    conn_mats=[]
    #f=open("new-dependences-"+str(n),'w')
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
                #parse bases into form for Matroid constructor
                for base in bases:
                    bases_list.append(reduce(lambda x,y:x+y,base))
                Mat=Matroid(groundset=matroid_ground,bases=bases_list)
                if Mat.is_valid() and Mat.is_connected():
                    #found connected matroid
                    CIs=fn_to_CIs(Mat.rank,n)
                    conn_mats.append(CIs)
                    print [reduce(lambda x,y:x+y,base,'') for base in Mat.bases()]
                    gaussoid= checkall(CIs,True)
                    
                    if gaussoid:
                        num_gaussoid+=1
                        
    print num_gaussoid
    return conn_mats


#ij|L, ik|jL \implies ik|L, ij|kL
def checkGaussoid7(CIs,output):
    pairCIs=findsubsets(CIs,2)
    for pair in pairCIs:
        (i1,j1,L1)=pair[0]
        (i2,k2,L2)=pair[1]
        (L1,L2)=(set(L1),set(L2)) 
        if i1==i2:
            #(i1,j1,L1) is the first CI statement
            if (not L2=='') and (L1.add(str(j1))==L2):
                #L = L1, need to remove j after checking
                L1.remove(str(j1))
                L=reduce(lambda x,y:x+y,sorted(reduce(lambda x,y:x+y, L1,'')),'')
                #kL = k\cup L1
                L1.add(str(k2))
                kL=reduce(lambda x,y:x+y,sorted(reduce(lambda x,y:x+y, L1,'')),'')
                L1.remove(str(k2))
                if (not (i1,k2,L) in CIs) or (not (i1,j1,kL) in CIs):
                    if output:
                        print "failed 7"
                    return False
                    
            #(i2,k2,L2) is the "first" CI statement
            if (not L1=='') and (L2.add(str(k2))==L1):
                #L = L2, need to remove k after checking
                L2.remove(str(k2))
                L=reduce(lambda x,y:x+y,sorted(reduce(lambda x,y:x+y, L2,'')),'')
                #kL = k\cup L2
                L2.add(str(j1))
                kL=reduce(lambda x,y:x+y,sorted(reduce(lambda x,y:x+y, L2,'')),'')
                L2.remove(str(j1))
                if (not (i1,j1,L) in CIs) or (not (i1,k2,kL) in CIs):
                    if output:
                        print "failed 7"
                    return False
    return True

#ij|kL, ik|jL  \implies ij|L, ik|L
def checkGaussoid8(CIs,output):
    pairCIs=findsubsets(CIs,2)
    for pair in pairCIs:
        (i1,j1,L1)=pair[0]
        (i2,k2,L2)=pair[1]
        (L1,L2)=(set(L1),set(L2))
        if (not L1=='') and (not L2=='') and (i1==i2) and (str(j1) in L2) and (str(k2) in L1) and (L2.remove(str(j1))==L1.remove(str(k2))):
            #ij|kL and ik|jL
            #should have ij|L and ik|L
            L=reduce(lambda x,y:x+y,sorted(reduce(lambda x,y:x+y, L2,'')),'')
            if (not (i1,j1,L) in CIs) or (not (i1,k2,L) in CIs):
                if output:
                    print "failed 8"
                return False
    return True

#ij|L, ik|L \implies ij|kL, ik|jL
def checkGaussoid9(CIs,output):
    pairCIs=findsubsets(CIs,2)
    for pair in pairCIs:
        (i1,j1,L1)=pair[0]
        (i2,k2,L2)=pair[1]
        (L1,L2)=(set(L1),set(L2))
        if (i1==i2) and (L1==L2):
            L1.add(str(j1))
            jL=reduce(lambda x,y:x+y,sorted(reduce(lambda x,y:x+y, L1,'')),'')
            L1.remove(str(j1))
            L1.add(str(k2))
            kL=reduce(lambda x,y:x+y,sorted(reduce(lambda x,y:x+y, L1,'')),'')
            if (not (i1,j1,kL) in CIs) or (not (i2,k2,jL) in CIs):
                if output:
                    print "failed 9"
                return False
    return True

# ij|L, ij|kL \implies ik|L or jk|L
def checkGaussoid10(CIs,output):
    pairCIs=findsubsets(CIs,2)
    for pair in pairCIs:
        (i1,j1,L1)=pair[0]
        (i2,k2,L2)=pair[1]
        (L1,L2)=(set(L1),set(L2))
        if (i1==i2) and (j1==k2):
            sym_diff=L1.symmetric_difference(L2)
            if len(sym_diff)==1:
                L=reduce(lambda x,y:x+y,sorted(reduce(lambda x,y:x+y, L1.intersection(L2),'')),'')
                k=int(sym_diff.pop())
                if (not (i1,k,L) in CIs) and (not (j1,k,L) in CIs) and (not (k,j1,L) in CIs):
                    if output:
                        print "failed 10"
                    return False
    return True

def checkall(CIs,output):
    if output:
        c7=checkGaussoid7(CIs,output)
        c8=checkGaussoid8(CIs,output)
        c9=checkGaussoid9(CIs,output)
        c10=checkGaussoid10(CIs,output)
        return c7 and c8 and c9 and c10
    return checkGaussoid7(CIs,output) and checkGaussoid8(CIs,output) and  checkGaussoid9(CIs,output) and checkGaussoid10(CIs,output)
   
def findsubsets(S,m):
    #return all m element subsets of S
    return list(itertools.combinations(S,m))
    
def subset_to_fn(subset):
    K=set(subset)
    def f(I):
        return min(1,len(K.intersection(set(I))))
    return f
    
def setfns_to_fn(setfns):
    def g(H):
        return reduce(lambda x,y:x+y, [w_ki(set(H)) for w_ki in setfns])
    return g
    
subsets = list(itertools.chain.from_iterable(itertools.combinations(['1','2','3','4'], r) for r in range(1,5)))

mss_list = list(itertools.chain.from_iterable(itertools.combinations(subsets, r) for r in range(1,16)))

#uncomment for MSS
"""total=0
gaussoids = 0
temp = set()
for mss in mss_list:
    w_K = setfns_to_fn([subset_to_fn(set(x)) for x in mss])
    CIs = fn_to_CIs(w_K,4)
    if temp.issubset(CIs) and temp.issuperset(CIs):
        continue
    else:
        total +=1
        if checkall(CIs,False):
            gaussoids+=1
        temp = CIs
print gaussoids, total"""


#uncomment for all intersections
"""allCIs=list(itertools.chain.from_iterable([checkGaussoids(i) for i in xrange(2,int(sys.argv[-1])+1)]))
longlist=list(itertools.chain.from_iterable([findsubsets(allCIs,i) for i in xrange(1,len(allCIs)+1)]))
num_gaussoid=0
total=0
for MSMP in longlist:
    total+=1
    unionCIs=MSMP[0]
    for matroid in MSMP:
        unionCIs=unionCIs.intersection(matroid)
    if checkall(unionCIs,False):
        num_gaussoid+=1
print num_gaussoid,total"""

checkGaussoids(int(sys.argv[-1]))
