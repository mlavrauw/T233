#############################################################################
##
##  RankDistribution.gi              T233 package
##                                                            Nour Alnajjarine
##                                                            Michel Lavrauw
##                                                        
##
##  Copyright 2019	Sabanci University
##                  
##  Implementation stuff for rank distributions of contraction spaces
##
#############################################################################

InstallGlobalFunction(MatrixOfPoint,function(x,m,n)
# turns a point of a projective space into an (mxn)-matrix containing the coordinates
	local list,i;
	list:=Coordinates(x);
	return List([1..m],i->list{(i-1)*n+[1..n]});
end);

InstallGlobalFunction(RankOfPoint,function(x,m,n)
# returns the rank of MatrixOfPoint(x,m,n)
    local list,mat,i;
	mat:=MatrixOfPoint(x,m,n);
	return RankMat(mat);
end);

InstallGlobalFunction(RankDistribution,function(subspace,m,n)
	local list;
	list:=List(Points(subspace),x->RankOfPoint(x,m,n));
	return Collected(list);
end);


