import itertools

# rank functions from slide 17 of http://www.ece.drexel.edu/walsh/Yunshu_CIR.pdf
# should be non-matroidal, extreme rays of submodular cone

#for matroids
def f(i):
    def rank(X):
    	if X==frozenset([i]):
	   return 2
    	return min(len(X),2)
    return rank

def g(i):
    def rank(X):
    	if i in X:
           return min(3,len(X)+1)
        return len(X)
    return rank

def h(i,j):
    def rank(X):
    	if len(X)==2 and not X==frozenset([i,j]):
           return 3
        return min(4,2*len(X))
    return rank

def checkMatroidal():
    for i in xrange(1,5):
    	M=Matroid('1234',rank_function=f(str(i)))
	print M.is_valid()
	M=Matroid('1234',rank_function=g(str(i)))
	print M.is_valid()
    temp=['12','13','14','23','24','34']
    for i in temp:
    	M=Matroid('1234',rank_function=h(i[0],i[1]))
	print M.is_valid()

#simpler to work with

def F(i):
    def rank(X):
    	if X==i:
	   return 2
    	return min(len(X),2)
    return rank

def G(i):
    def rank(X):
    	if i in X:
           return min(3,len(X)+1)
        return len(X)
    return rank

def H(i,j):
    def rank(X):
    	if len(X)==2 and not X==i+j and not X==j+i:
           return 3
        return min(4,2*len(X))
    return rank

#f(Ki)+f(Kj) >= f(Kij)+f(K)

#if it is = then we have i \independent j | K

def fn_to_CIs():
    #f=open("CIs-fij",'w')
    f=open("CIs-gi2_readable",'w')
    #for temp in ['12','13','14','23','24','34']:
    	#r=H(temp[0],temp[1])
    for i in xrange(1,5):
        r=F(str(i))
        #f.write("rank function:f_"+temp+"\n")
	f.write("rank_function:g_"+str(i)+"\n")
        #f.write("{")
	#output=''
	for i in xrange(1,5):
	    for j in xrange(i+1,5):
	    	K=('1234'.replace(str(i),'')).replace(str(j),'')
		subsets=[findsubsets(K,x) for x in xrange(0,3)]
		for size in subsets:
		    for subset in size:
		    	K=reduce(lambda x,y:x+y,subset,'')
			if r(K+str(i))+r(K+str(j))==r(K+str(i)+str(j))+r(K):
			    f.write(str(i)+" "+str(j)+" "+str(K)+"\n")
                            #output+="{{"+str(i)+","+str(j)+"},{"+K+"}},"
	#f.write(output[:-1]+"}\n")
        f.write('\n')
        f.flush()
    f.close()


def findsubsets(S,m):
    #return all m element subsets of S
    return list(itertools.combinations(S,m))

#checkMatroidal()

fn_to_CIs()
