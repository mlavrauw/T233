
#############################################################################
##
##  Rank2,3PointsOfContractionSpaces.gd         T233 package
##                                                            Nour Alnajjarine
##                                                            Michel Lavrauw

##
##  Copyright 2019	Sabanci University
##                
##
##  Implementation stuff for rank 2 and 3 points of contraction subspaces
##
#############################################################################

InstallGlobalFunction(Rank2PtsOftheContractionSubspace,function(Ai,m,n)
	local x, list;
	list:=[];
	for x in Points(Ai) do
	if RankOfPoint(x,m,n)=2 then list:=Union(list,[x]);
	fi;
od;
return(list);
end);

InstallGlobalFunction(Rank3PtsOftheContractionSubspace,function(Ai,m,n)
	local x, list;
	list:=[];
	for x in Points(Ai) do
	if RankOfPoint(x,m,n)=3 then list:=Union(list,[x]);
	fi;
od;
return(list);
end);