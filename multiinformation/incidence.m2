load "5vertex.vertices"  -- should get 77 vertices and 22 facets
load "5vertex.ineqs"

load "6vertex.vertices"  -- should get 422 vertices and 46 facets
load "6vertex.ineqs"

load "6vertexB.vertices"  -- should get 427 vertices and 45 facets
load "6vertexB.ineqs"

load "7vertex.vertices"  -- should get 3196 vertices and 102 facets
load "7vertex.ineqs"

load "7vertexB.vertices"  -- should get 2361 vertices and 83 facets
load "7vertexB.ineqs"
load "7vertexB.vertices"

#V
#I

verts = matrix V;
ineqs = matrix drop(I,-2); -- drop the last two inequalities, which are equations
slack = entries(verts * transpose ineqs);

maxSlack = max apply(slack, v-> max(v/abs)) -- max absolute value of entry
BIG = 2^35/maxSlack

minSlack = min apply(BIG*slack, v-> min v) -- should be <= 0
epsilon = max(1,-1*minSlack) -- usually 1

incidence = apply(BIG*slack, r-> apply(r, i-> -- rounding happens here
	(if i > epsilon then 1 else if i < -epsilon then -1 else 0)));

FthruV = unique incidence; -- facets thru vertices
<< endl << "number of vertices = " << #FthruV << endl;

VinF = unique entries matrix transpose FthruV; -- vertices in facets

-- look for maximal faces (i.e. remove redundant inequalities)
FF = select(VinF, f -> (
	maximal = true;
	for g in VinF do 
	  if min(f-g) >= 0 and g =!= f 
	  then (maximal = false; break); 
	maximal
	)
    ) ;   

<< endl<< "number of facets = " << #FF << endl; 
<< endl<< "number of vertices = "<<#unique entries transpose matrix FF<<endl;




