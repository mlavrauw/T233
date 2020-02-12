
#############################################################################
##
##  SubspaceOfContractions.gi              T233 package
##                                                            Nour Alnajjarine
##                                                            Michel Lavrauw
##                                                        
##
##  Copyright 2019	Sabanci University
##                  
##  Implementation stuff for contraction subspaces
##
#############################################################################

InstallGlobalFunction(CubicalArrayFromPointInTensorProductSpace,function(point,n1,n2,n3)
# this obvious depends on how we choose the coordinates
local vec,array,i,cube,lijst;
	if not Length(Coordinates(point))=n1*n2*n3 then
		Error("Incompatible arguments","\n");
	fi;
	vec:=Coordinates(point);
	array:=List([1..n1],i->vec{[1+(n2*n3)*(i-1)..(n2*n3)*i]});
	cube:=[];
	for lijst in array do
		Add(cube,List([1..n2],i->lijst{[1+n3*(i-1)..n3*i]}));
	od;
	return cube;
end);

InstallGlobalFunction(ContractionOfPointInTensorProductSpace,function(vec,position,point,n1,n2,n3)
local fld,q,cube,mat,mats,i,j,k;
	fld:=BaseField(AmbientSpace(point));
	q:=Size(fld);
	cube:=CubicalArrayFromPointInTensorProductSpace(point,n1,n2,n3);
	if position=1 then
		mat:=Sum([1..n1],i->vec[i]*cube[i]);
		return(VectorSpaceToElement(PG(n2*n3-1,q),Flat(mat)));
	else
		if position=2 then
			mats:=[];
			for j in [1..n2] do
				Add(mats,List([1..n1],i->cube[i][j]));
			od;
			mat:=Sum([1..n2],j->vec[j]*mats[j]);
			return(VectorSpaceToElement(PG(n1*n3-1,q),Flat(mat)));
		else # in this case position is 3
			mats:=[];
			for k in [1..n3] do
				Add(mats,List([1..n1],i->TransposedMat(cube[i])[k]));
			od;
			mat:=Sum([1..n3],k->vec[k]*mats[k]);
			return(VectorSpaceToElement(PG(n1*n2-1,q),Flat(mat)));
		fi;
	fi;
end);

# Next we make a subspace in a two-fold tensor product space from a tensor in a 3-fold tensor space
InstallGlobalFunction(SubspaceOfContractions,function(position,point,n1,n2,n3)
local nlist,n,q,V,bas,list,vec;
	nlist:=[n1,n2,n3];
	n:=nlist[position];
	q:=Size(BaseField(AmbientSpace(point)));
	V:=GF(q)^n;
	bas:=Basis(V);
	list:=List(bas,
        vec->ContractionOfPointInTensorProductSpace(vec,position,point,n1,n2,n3));
	return(Span(list));
end);
