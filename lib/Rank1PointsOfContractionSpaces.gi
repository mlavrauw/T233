
#############################################################################
##
##  Rank1PointsOfContractionSpaces.gi         T233 package
##                                                            Nour Alnajjarine
##                                                            Michel Lavrauw

##
##  Copyright 2019	Sabanci University
##                
##
##  Implementation stuff for rank one points of contraction subspaces
##
#############################################################################

InstallGlobalFunction(Rank1PtsOftheContractionSubspace,function(Ai,m,n)
	local x, list;
	list:=[];
	for x in Points(Ai) do
	if RankOfPoint(x,m,n)=1 then list:=Union(list,[x]);
	fi;
od;
return(list);
end);
