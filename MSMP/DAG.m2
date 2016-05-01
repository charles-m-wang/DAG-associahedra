--restart
loadPackage "GraphicalModels";
loadPackage "Posets";


--------------------
-- FUNCTIONS
--------------------
-- input: a poset given as a list of linear extensions
-- output: non-trivial relations in the poset 
--    	   (intersection of input permutations)

rels = P -> (
    sort toList fold((c,d)->c*d,  -- "*" = set intersection
 	apply(P, p -> fold((a,b)->a+b, -- "+" = set union 
		 flatten apply(0..#p-2, i-> 
		     set apply(i+1..#p-1, j-> {p_j,p_i}))
		 )
	     )
	 )
    );

------------------------
-- input: a poset given as a list of linear extensions
-- output: the facets containing the vertex associated 
--    	   to the input poset in the generalized permutohedron
--         (the rays in the corresponding cone in the S_n fan)

facetsFromPoset = P -> ( 
    R = rels P;
    PP = poset R;
    H = graph coveringRelations PP;
    candidateFacets = unique sort flatten apply(P, p -> toList apply(1..#p-1, i->sort drop(p,i))); 
    select(candidateFacets, F -> (
	    H1 = inducedSubgraph(H, F);
	    H2 = inducedSubgraph(H, vertices H - set F);
	    isConnected(H1) and isConnected(H2)
	    )
	)
    );

----------------------------------------
-- END OF FUNCTIONS
---------------------------------------

 n=4; -- number of vertices of the DAG
-- G = digraph {{1,{3}}, {2,{3}},{3,{4}}};
-- G = digraph {{1,{4}},{2,{4}},{3,{4}}};
G = digraph {{1,{4}}, {2,{3}},{3,{4}}}
--G = digraph {{1,{2,3}}, {2,{4}},{3,{4}}}

--n = 5;
--G = digraph{{1,{3,4}}, {2,{3}},{3,{5}},{4,{5}}} -- 77 vertices and 22 facets in DAG associahedron

--n = 6;
--G = digraph{{1,{3,4}},{2,{3,6}},{3,{5}},{4,{5}},{6,{5}}} -- 422 vertices and 46 facets in DAG associahedron
--G = digraph{{1,{2,3}},{2,{3,4}},{3,{4,5}},{4,{6}},{5,{6}}} -- 427 vertices and 45 facets

n=7;
G = digraph{{1,{2,3,5}},{2,{6,7}},{3,{4}},{4,{5,7}},{5,{6}},{6,{7}}} --102 facets and 3196 vertices
G = digraph{{1,{4,7}},{2,{3,5}},{3,{5,6}},{4,{5,6}},{5,{7}}} --83 facets and 2361 vertices

----------------
-- find all conditional independent statements

CIs = globalMarkov G; -- needs to add more according to semigraphoid axioms
pairCIs = sort unique flatten flatten apply(CIs , S ->(
    flatten table(S#0, S#1, (i,j)-> (
	apply(subsets(set(S#0) + set (S#1) - set(i,j)), U -> {sort{i,j},sort toList(set(S#2) + U)})
	)
	)
    )
    );
 netList pack(2,sort pairCIs)

------------------
-- computing vertices of the DAG associahedron


V = permutations toList(1..n);
E = flatten apply(pairCIs, S ->
    (
        A := permutations S_1;
	C := permutations(toList(1..n)-(set(S_1)+set(S_0)));
        flatten table(A,C,(a,c)->{a|{S#0#0, S#0#1}|c, a|{S#0#1, S#0#0}|c})
	)
    );
VV = V/(i->v_i);
EE = E/(e->e/(i->v_i));
cc = connectedComponents graph(VV,EE);
print("The number of vertices = ", #cc); -- #vertices of the DAG associahedron
-- cc/(i->#i)
permClasses = cc/(C->(C/(i->i#1))); -- vertices of DAG associahedron
-- netList pack(1, sort permClasses) -- the poset corresponding to each vertex is given as the set of its linear extensions

-------------------
-- facets

F = unique sort flatten apply(permClasses, facetsFromPoset);
-- netList pack(3,F) -- facets of DAG associahedron
print("The number of facets = ", #F); -- #facets of the DAG associahedron

--
--nonFacets =  sort toList(set(subsets(toList(1..n))) - F);
--netList pack(3,nonFacets)
