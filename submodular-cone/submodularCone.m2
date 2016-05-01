
-------------------------------
toPolymakeFormat = M -> (
     if M === null then ""
     else(
     	  S := "";
     	  if numRows M > 0 then
	     S = S|replace("\\|", "", toString net M);
     	  S
     	  )
     )
----------------------------------------------------

n = 5;

P = subsets n;
R = ZZ[apply(P, S -> x_S)];
ineqs = {};
scan(n,i -> apply(i+1..n-1, j-> 
	scan(subsets(toList(0..n-1)-set{i,j}), K ->(
	    ineqs = ineqs |{apply(gens R, X -> coefficient(X, x_(sort(K|{i})) + x_(sort(K|{j})) - x_(sort(K|{i,j})) - x_K ))};
	    ))
	)
)
M = matrix ineqs
rank M
N = transpose gens kernel M
A = M||N||-N
A = (matrix apply(numRows A, i-> {0})) | A
"submodularReduced"|n|".poly" << "INEQUALITIES" << endl << toPolymakeFormat (A) << endl << close;

