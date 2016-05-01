
needsPackage "GraphicalModels";

-- Read in files containing lists of connected matroids.
-- Each matroid is encoded by the list of its bases 
-- and the list of conditional dependence statements
M_2 = (lines get "new-dependences-2")/value;
M_3 = (lines get "new-dependences-3")/value;
M_4 = (lines get "new-dependences-4")/value;
M_5 = (lines get "new-dependences-5")/value;
M_6 = (lines get "new-dependences-6")/value;


-- find all "connected" matroids, possibly with loops
genMatroids = n -> flatten apply(2..n, i -> 
    (
	S := subsets(1..n,i);
	flatten apply(S, s-> 
	    -- we will transport a matroid on i elements into subset s
	    apply(M_i, m-> -- m is a matroid
		(   
		    newBases := apply(m_0, l-> l/(j->s_(j-1))); 
		    loops := set(1..n)-s;
		    newCDs := flatten table(m_1, subsets loops, (l,L)->  
			-- l is a CD statement i \indep j | K
			-- L is extra elements on condition on
				{
			    	    l_0/(j->s_(j-1)), 
			    	    sort(toList(L)|l_1/(j->s_(j-1)))
			    } );
		    {newBases, newCDs}
		    )
	    	)
	    )
    )	 
) ;

allMatroids_3 = genMatroids 3;
allMatroids_4 = genMatroids 4;
allMatroids_5 = genMatroids 5;
allMatroids_6 = genMatroids 6;


---------------------
-- input: DAG
-- output: pairwise conditional independence statements

CIs = method()
CIs Digraph := G -> (
    -- find all conditional independent statements
    someCIs := globalMarkov G; 
    -- needs to add more according to semigraphoid axioms
    pairCIs := sort unique flatten flatten apply(someCIs , S ->(
    	    flatten table(S#0, S#1, (i,j)-> (
		    apply(subsets(set(S#0) + set (S#1) - set(i,j)), U -> 
			{sort{i,j},sort toList(set(S#2) + U)})
		    ))));
    return(pairCIs);
    );
    
---------------------------------    
-- input: pairwise CI statements
-- output: matroids compatible to the CI statements
-- each matroid is encoded by bases and CD statements
CIs2matroids = method(Options => {Print => true}) 
CIs2matroids (ZZ, List) := opt -> (n, pairCIs) -> (
    count := 0;
    compatibleMatroids := {};
    scan(allMatroids_n, m -> ( -- * means set intersection
    	    if #(set(m_1)*set(pairCIs)) == 0 then 
    	    (   if (opt.Print == true) then 
		(<< count << ". "; printMatroid(m);); --print compatible matroid	      
		count = count+1;
		compatibleMatroids = compatibleMatroids |{m});	
	    )
    	);
    return compatibleMatroids;
    );

---------------------------------
-- input: matroid info, consisting of bases and a list of conditional dependences
-- This function prints the rank, the number bases, and the list of a bases.
printMatroid = m -> (
    << "Rank= " << #(m#0#0) << ". #Bases= " << #(m#0) << ". Bases: "<< apply(m#0, a->concatenate(a/toString)) <<endl; 
    );

printMatroidList = M -> scan(#M,i-> (<< i << ". "; printMatroid(M_i);));

------------------------------
-- input: CI statements, all compatible matroids, list of indices of matroids
-- output: true of the matroids are sufficient for constructing the DAG asociahedron

sufficient = (n, pairCIs, compatibleMatroids, selection) -> (
    condDepFromMatroids := {};
    if(#compatibleMatroids > 0) then (
    	condDepFromMatroids = fold((a,b)->a+b,
     	    apply(compatibleMatroids_(toList selection), m-> set(m_1))
     	    );
    );
    -- check that we get all CI/CD statements    
    return(#condDepFromMatroids + #pairCIs == binomial(n,2)*2^(n-2));
)

----------------------------------
-- "missingCDs" is like the function "sufficient" but returns a list of conditional dependences that are not covered by the matroids. (This probably runs more slowly than "sufficient".)

allCIs = n -> (
    L := {};
    scan(1..n-1, i -> 
	scan(i+1..n, j-> (
		scan(subsets(set(1..n)-set{i,j}), K -> 
		    (L = L | {{{i,j},sort toList K}};)
		    )
		)
   	)
    );
    return L;
)

missingCDs = (n, pairCIs, compatibleMatroids, selection) -> (
    condDepFromMatroids := {};
    if(#compatibleMatroids > 0) then (
    	condDepFromMatroids = fold((a,b)->a+b,
     	    apply(compatibleMatroids_(toList selection), m-> set(m_1))
     	    );
    );
    -- return conditional dependences that are not covered by selected matroids
    return (toList( set(allCIs n) - (set pairCIs) - condDepFromMatroids));
)



------------------
-- method function to check whether the DAG associahedron (or a collection of CI statements)
-- is a Minkowski sum of matroid polytopes
-- input G is either a Digraph or a List of pairwise CI statements

checkMSMP = method(Options => {Print => true}) 

checkMSMP (ZZ, List) := opt -> (n,pairCIs) -> ( 
    -- n is the size of ground set (# nodes in the di-graph)
    -- pairCIs is a list of elements of the form {{1,2},{3,4,5}}
    -- meaning "1 is independent of 2 given {3,4,5}"
    compatibleMatroids := CIs2matroids(n,pairCIs,opt);
    return (sufficient(n, pairCIs, compatibleMatroids, 0..#compatibleMatroids-1));
);

checkMSMP Digraph := opt-> G -> (
    return(checkMSMP(#(vertexSet G),CIs G, opt));
    );

-----------------------------------------------------



end;

-----------------------------------------------------
-----------------------------------------------------
-- START HERE ---------------------------------------
-----------------------------------------------------

restart
load "matroids.m2";

 G = digraph {{1,{3}}, {2,{3}},{3,{4}}};
 G = digraph {{1,{4}},{2,{4}},{3,{4}}};
 G = digraph {{1,{4}}, {2,{3}},{3,{4}}}
G = digraph {{1,{2,3}}, {2,{4}},{3,{4}}}
G = digraph {{1,{3,4}},{2,{3}},{3,{4}}}

G = digraph{{1,{3,4}}, {2,{3}},{3,{5}},{4,{5}}} -- 77 vertices and 22 facets in DAG associahedron
G = digraph{{1,{3}},{2,{3,4}},{3,{5}},{4,{5}} }
G = digraph{{1,{2,3}},{2,{5}},{3,{5}},{4,{5}} }
G = digraph{{1,{2}}, {2,{5}}, {3,{4}}, {4,{5}}}
G = digraph{{1,{3,4}}, {2,{5}},{3,{4}},{4,{5}}}

G = digraph{{1,{3,4}},{2,{3,6}},{3,{5}},{4,{5}},{6,{5}}} -- 422 vertices and 46 facets in DAG associahedron
G = digraph{{1,{2,3}},{2,{3,4}},{3,{4,5}},{4,{6}},{5,{6}}} -- 427 vertices and 45 facets
G = digraph{{1,{3}},{2,{3,4}},{4,{6}},{5,{6}}}
G = digraph{{1,{4,6}},{2,{4,5}},{3,{5,6}}}
G = digraph{{1,{4}},{2,{4,6}},{3,{5,6}},{4,{5,6}}}
G = digraph{{1,{4}},{2,{4,5}},{3,{5}},{4,{5}},{5,{6}}}
G = digraph{{1,{4,5}},{2,{6}},{3,{4,6,7}},{4,{5,6}},{6,{7}}}

checkMSMP G

scan(#compatibleMatroids, i -> (
	m = compatibleMatroids_i; 
	print(i,#(m#0#0), #m#0);
	))

-------------------------------------------------
-- GAUSSOIDS -----------------------------------
-- representable Gaussioids from the paper of Drton & Xiao
-- The rest come from undirected for DAG models, and
-- we've already checked all DAGS with 4 nodes
-------------------------------------------------

restart
load "matroids.m2";

gaussoid_13 = { {{1,2},{}}, {{1,4},{2,3}}, {{2,3},{1,4}} };

gaussoid_15 = { {{1,4},{}}, {{1,4},{2,3}}, {{2,3},{}}, {{2,3},{1,4}} };

gaussoid_28 = { {{1,3},{2}}, {{1,4},{}}, {{2,3},{1,4}}, {{2,4},{3}} };

gaussoid_35 = { {{1,2},{3}}, {{3,4},{1}} };

gaussoid_40 = { {{1,4},{}}, {{1,4},{3}}, {{2,3},{}}, {{2,3},{4}}, {{3,4},{}}, {{3,4},{1}}, {{3,4},{2}} };

gaussoid_43 = { {{2,3},{}}, {{2,3},{4}}, {{2,4},{1,3}}, {{3,4},{}}, {{3,4},{2}}};

gaussoid_44 = { {{1,4},{}}, {{2,3},{}} };

gaussoid_48 = { {{1,2},{4}}, {{2,3},{}} };

gaussoid_49 = { {{1,4},{3}}, {{2,3},{}} };

all({13,15,28,35,40,43,44,48,49}, i -> checkMSMP(4,gaussoid_i))

gaussoid_9 = { {{2,3},{1}}, {{2,3},{1,4}}, {{2,4},{3}}, {{3,4},{1}}, {{3,4},{1,2}} };

gaussoid_14 = { {{1,4},{}}, {{1,4},{2,3}}, {{2,3},{1,4}} };

gaussoid_15 = { {{1,4},{}}, {{1,4},{2,3}}, {{2,3},{}}, {{2,3},{1,4}} };

gaussoid_17 = { {{1,2},{3}}, {{2,3},{1,4}} };

gaussoid_19 = { {{1,2},{3}}, {{1,4},{2}}, {{2,3},{1,4}} }; -- not MSMP

gaussoid_20 = { {{1,4},{2}}, {{1,4},{3}}, {{2,3},{1,4}} }; -- not MSMP

gaussoid_21 = { {{1,3},{2}}, {{2,3},{1,4}}, {{2,4},{3}} };

gaussoid_23 = { {{1,2},{}}, {{1,4},{3}}, {{2,3},{1,4}} };

gaussoid_24 = { {{1,2},{}}, {{2,3},{1,4}}, {{2,4},{3}} };

gaussoid_25 = { {{1,2},{}}, {{2,3},{1,4}}, {{3,4},{2}} };

gaussoid_27 = { {{1,2},{3}}, {{1,4},{}}, {{2,3},{1,4}} };

gaussoid_28 = { {{1,3},{2}}, {{1,4},{}}, {{2,3},{1,4}}, {{2,4},{3}} };

gaussoid_29 = { {{2,3},{}}, {{2,3},{1,4}} };

gaussoid_30 = { {{1,4},{2}}, {{2,3},{}}, {{2,3},{1,4}} };

gaussoid_32 = { {{1,2},{3}}, {{1,2},{4}} };

gaussoid_33 = { {{1,2},{3}}, {{1,3},{4}} };

gaussoid_34 = { {{1,2},{3}}, {{1,4},{2}}, {{2,3},{4}} }; -- not MSMP

gaussoid_36 = { {{1,2},{3}}, {{1,2},{4}}, {{3,4},{1}} };

gaussoid_37 = { {{1,2},{3}}, {{1,2},{4}}, {{3,4},{1}}, {{3,4},{2}} };

gaussoid_42 = { {{2,3},{}}, {{2,3},{4}}, {{2,4},{1}}, {{3,4},{}}, {{3,4},{2}} };

gaussoid_45 = { {{1,2},{3,4}}, {{1,4},{}}, {{2,3},{}} };

gaussoid_46 = { {{1,4},{}}, {{1,4},{2,3}}, {{2,3},{}} };

gaussoid_50 = { {{1,2},{4}}, {{1,4},{3}}, {{2,3},{}} }; -- not MSMP

gaussoid_51 = { {{1,4},{2}}, {{1,4},{3}}, {{2,3},{}} }; -- not MSMP

gaussoid_52 = { {{1,3},{4}}, {{2,3},{}}, {{2,4},{1}} };

toCheck = {9,14,15,17,19,20,21,23,24,25,27,28,29,30,32,33,34,36,37,42,45,46,50,51,52};
MSMP = select(toCheck, i -> checkMSMP(4,gaussoid_i))
notMSMP = sort toList( set(toCheck) - set(MSMP))  -- {19, 20, 34, 50, 51}
scan(notMSMP, i -> (<< endl << i << ". CI statements:" << endl << gaussoid_i << endl<< "Matroids:" << endl << checkMSMP(4,gaussoid_i) << endl;)) 

-------------------------------------------------------------
-- NON-SUBMODULAR EXAMPLES from Three Counterexamples on Semigraphoids
-- Raymond Hemmecke, Jason Morton, Anne Shiu,
-- Bernd Sturmfels, and Oliver Wienand
-- http://arxiv.org/pdf/math/0610451v1.pdf
-------------------------------------------------------------

ex1 = { {{2,3},{1,4}}, {{1,4},{2,3}}, {{1,2},{}}, {{3,4},{}}};
checkMSMP(4,ex1)

ex2 = { {{1,5},{}}, {{2,3},{1,5}}, {{2,3},{1,4,5}}, {{2,4},{1,5}}, {{2,4},{1,3,5}}, {{3,4},{1,5}}, {{3,4},{1,2,5}}, {{1,2},{}}, {{2,5},{}}, {{1,3},{}}, {{3,5},{}}, {{1,4},{}}, {{4,5},{}}, {{2,3},{}}, {{2,4},{}}, {{2,3},{4}}, {{2,4},{3}}, {{1,5},{234}}, {{3,4},{}}, {{3,4},{2}}, {{1,2},{3}}, {{1,3},{2}}, {{2,3},{1}}, {{4,5},{123}}, {{2,5},{3}}, {{2,3},{5}}, {{3,5},{2}}, {{1,4},{2,3,5}}, {{3,4},{5}}, {{3,5},{4}}, {{4,5},{3}}, {{1,2},{345}}, {{1,2},{4}}, {{1,4},{2}}, {{2,4},{1}}, {{3,5},{124}}, {{2,4},{5}}, {{2,5},{4}}, {{4,5},{2}}, {{1,3},{245}}, {{1,4},{3}}, {{3,4},{1}}, {{1,3},{4}}, {{2,5},{134}} };
checkMSMP(5,ex2)

-------------------------------------------------------------
-- Check whether all DAGs up to 6 nodes are MSMP
-------------------------------------------------------------
restart
load "matroids.m2"
-- Comment out the "print" line in checkMSMP function

n= 4; 
allDAGS = toList(( 
	(set apply(subsets(toList(2..n)), S -> {1,S}))**
	(set apply(subsets(toList(3..n)), S -> {2,S}))**
	(set apply(subsets(toList(4..n)), S -> {3,S})) )/deepSplice@@toList);
#allDAGS


n= 5; 
allDAGS = toList(( 
	(set apply(subsets(toList(2..n)), S -> {1,S}))**
	(set apply(subsets(toList(3..n)), S -> {2,S}))**
	(set apply(subsets(toList(4..n)), S -> {3,S}))**
	(set apply(subsets(toList(5..n)), S -> {4,S}))
	)/deepSplice@@toList);
#allDAGS

n= 6;
allDAGS = toList(( 
	(set apply(subsets(toList(2..n)), S -> {1,S}))**
	(set apply(subsets(toList(3..n)), S -> {2,S}))**
	(set apply(subsets(toList(4..n)), S -> {3,S}))**
	(set apply(subsets(toList(5..n)), S -> {4,S}))**
	(set apply(subsets(toList(6..n)), S -> {5,S}))
	)/deepSplice@@toList);
#allDAGS

count = 0;
all(allDAGS, G -> (
	msmp = checkMSMP (digraph G, Print => false);
	if not msmp then print G; 
	count = count+1; 
	if count%100 == 0 then print count;
	msmp
	)
    )

-----------------------
---- n=7 --------------
restart
load "matroids.m2"

n=7;
-- generate matroids on {1,2,...,7} coming from connected matroids on 6 elements
allMatroids_7 = flatten apply(2..6, i -> 
    (
	S := subsets(1..n,i);
	flatten apply(S, s-> 
	    apply(M_i, m-> 
		(   
		    newBases := apply(m_0, l-> l/(j->s_(j-1))); 
		    loops := set(1..n)-s;
		    newCDs := flatten table(m_1, subsets loops, (l,L)->      	
			{
			    	    l_0/(j->s_(j-1)), 
			    	    sort(toList(L)|l_1/(j->s_(j-1)))
			    } );
		    {newBases, newCDs}
		    )
	    	)
	    )
    )	 
) ;

------------------------------------------------
distinctMatroids_n = (lines get ("connected-dependences-"|toString(n)))/value;
------------------------------------------------

G = digraph {{1,{2}},{2,{4}},{3,{4,7}},{6,{7}},{5,{6}}}
G = digraph {{1,{4,5}}, {2,{5,7}},{3,{6,7}},{5,{6,5}}} --ex1
G = digraph {{1,{2,5,7}},{2,{4}},{3,{5}},{5,{6}},{6,{7}}}--ex2
G = digraph {{1,{3,7}},{2,{4}},{3,{4,5}},{4,{5}},{5,{6}}}--ex3
G = digraph {{1,{5}},{2,{4,5,7}},{3,{4,6}},{4,{5}},{5,{7}}}--ex4
G = digraph {{1,{3,5}},{2,{3}},{3,{7}},{4,{5,6,7}},{5,{7}}}--ex5
G = digraph {{1,{4,5}},{2,{6}},{3,{4,6,7}},{4,{5,6}},{6,{7}}}--ex6

all(1000, a-> (
	-- random DAG on n elements
	G = digraph apply(toList(1..n-1), i -> 	    {i,select(toList(i+1..n), j-> random(0,1) == 1)} );
	pairCIs = CIs G;
	foundMatroids = CIs2matroids (n, pairCIs, Print=>false); -- this uses only matroids coming from connected matroids on subsets of size <= 6
	suff=sufficient(n,pairCIs,foundMatroids, (0..#foundMatroids-1));
       	-- if suff then << a << ". " << G << endl << "6 points suffice." << endl << endl;
	if not suff then (
	    << a << ". " << G << endl << "Found " << #foundMatroids <<" compatible matroids." << endl;
	    << "Missing "<< missingCDs(n,pairCIs,foundMatroids,0..#foundMatroids-1) << endl;
	    -- check connected matroids on n=7 elements
	    -- generate them but don't store in memory
	    count = #foundMatroids;
	    testedCount = 0;
	    scan(distinctMatroids_n, m -> (
		    scan(permutations toList(1..n), p -> (
			    -- testedCount = testedCount+1; 
			    -- if(testedCount%1000==0) then << testedCount << endl;    
			    bases = apply(m#0, b -> sort apply(b, i-> p#(i-1)));      
			    CDs = apply(m#1, L -> apply(L, l -> sort apply(l, i-> p#(i-1))));	  
			    -- <<"testing " << {bases,CDs} << endl;
			    if #(set(CDs)*set(pairCIs)) == 0 then 
			    (  
				bases = sort bases;
				CDs = sort CDs;
				if not member({bases,CDs}, set foundMatroids) then (
				    printMatroid({bases, CDs});
				    count = count+1;
				    foundMatroids = foundMatroids |{ {bases,CDs} };      	         );  
		 		--suff=sufficient(n,pairCIs,foundMatroids, (0..#foundMatroids-1));
--		      		if suff then break;
		      		);	
		  	    ));
--		    if suff then (
--			<< "That's enough. " << #foundMatroids << " total matroids suffice." << endl << endl; break;)
		    ));
	    );	
	suff=sufficient(n,pairCIs,foundMatroids, (0..#foundMatroids-1));
	if not suff then ( << G << endl << " not MSMP!! ------------ " << endl; ) else ( << a << ". MSMP" << endl;);
	suff
	)
    )
-----------------------------------------------------
-- make "interactive" selection of matroids 
-----------------------------------------------------
restart
load "matroids.m2"

G = digraph {{1,{2,3}},{2,{4}},{3,{4}}}

G = digraph {{1,{3,4}},{2,{3}},{3,{4}}} -- this example shows that the "nice" induced subgraphs do not suffice
C = CIs G;
M = CIs2matroids (#(vertexSet G),C, Print=>true);
missingCDs(#(vertexSet G), C, M, 0..7)

sufficient(#(vertexSet G), C, M, set(0..9)-set{8})
missingCDs(#(vertexSet G), C, M, (0..7)|(9..9))
missingCDs(#(vertexSet G), C, M, (0..6)|(9..10))

G = digraph {{1,{3}},{2,{3}},{3,{5}},{4,{5}}}
G = digraph {{1,{4}},{2,{4}},{3,{4}},{4,{5}}}
G = digraph {{1,{3,4}},{2,{3,4}},{3,{5}},{4,{5}}}
G = digraph{{1,{4}},{2,{4,6}},{3,{5,6}},{4,{5,6}}}

C = CIs G;
M = CIs2matroids (#(vertexSet G),C);
sufficient(#(vertexSet G), C, M, 0..66)

printMatroidList M

checkMSMP G

