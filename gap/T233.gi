#
# T233: Algorithms for tensors.
#
# Implementations
#



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

InstallGlobalFunction(Rank1PtsOftheContractionSubspace,function(Ai,m,n)
	local x, list;
	list:=[];
	for x in Points(Ai) do
		if RankOfPoint(x,m,n)=1 then list:=Union(list,[x]);
		fi;
	od;
	return(list);
end);


InstallGlobalFunction(RepO10odd,function(q)
	local eta,v1,v2,ext;
	eta:=PrimitiveElement(GF(q));
	v1:=[1,0,0,0,eta,0,0,0,0]*One(GF(q));
	v2:=[0,1,0,1,0,0,0,0,0]*One(GF(q));
	ext:=[MatrixOfPoint(VectorSpaceToElement(PG(8,q),v1),3,3),MatrixOfPoint(VectorSpaceToElement(PG(8,q),v2),3,3)];
	return ext;
end);

InstallGlobalFunction(AlternativeRepresentationOfFiniteFieldElements,function(coeff,q)
	local w,list,c;
	w:=PrimitiveElement(GF(q));
	list:=[];
	for c in coeff do
		if c=Zero(GF(q)) then
			Add(list,Zero(GF(q)));
		else
			Add(list,LogFFE(c,w));
		fi;
	od;
	return List(list,i->w^i);
end);

InstallGlobalFunction(RepO10even,function(q)
	local f, coeff, c, u, v, v1, v2, ext;
  f:=MinimalPolynomial(GF(q),Z(q^2));
  coeff:=CoefficientsOfUnivariatePolynomial(f);
  c:=AlternativeRepresentationOfFiniteFieldElements(coeff,q);
  v:=-1/c[1];
  u:=c[2];
  v1:=[v,0,0,0,1,0,0,0,0]*One(GF(q));
  v2:=[0,1,0,1,u,0,0,0,0]*One(GF(q));
  ext:=[MatrixOfPoint(VectorSpaceToElement(PG(8,q),v1),3,3),MatrixOfPoint(VectorSpaceToElement(PG(8,q),v2),3,3)];
	return ext;
end);


InstallGlobalFunction(RepO15odd,function(q)
	local eta,v1,v2,ext;
	eta:=PrimitiveElement(GF(q));
	v1:=[1,0,0,0,eta,0,0,0,1]*One(GF(q));
	v2:=[0,1,0,1,0,0,0,0,0]*One(GF(q));
	ext:=[MatrixOfPoint(VectorSpaceToElement(PG(8,q),v1),3,3),MatrixOfPoint(VectorSpaceToElement(PG(8,q),v2),3,3)];
	return ext;
end);


InstallGlobalFunction(RepO15even,function(q)
        local f, coeff, c, u, v, v1, v2, ext;
  f:=MinimalPolynomial(GF(q),Z(q^2));
  coeff:=CoefficientsOfUnivariatePolynomial(f);
  c:=AlternativeRepresentationOfFiniteFieldElements(coeff,q);
  v:=-1/c[1];
  u:=c[2];
  v1:=[v,0,0,0,1,0,0,0,1]*One(GF(q));
  v2:=[0,1,0,1,u,0,0,0,0]*One(GF(q));
  ext:=[MatrixOfPoint(VectorSpaceToElement(PG(8,q),v1),3,3),MatrixOfPoint(VectorSpaceToElement(PG(8,q),v2),3,3)];
	return ext;
end);


InstallGlobalFunction(OrbitOfTensor,function(A)
	local fld, q, A1, A2, A3, R1, R2, R3, orbit, rank, representative, p, x1, x2, M2, v1, v2, B, sv, sm,
	t, z, y1, y2, z1, z2, p1, Ip1,	p2, Ip2, p3, Ip3, p4, Ip4, V, U, x, y, M1, Mw, w, K, F;
	fld:=BaseField(AmbientSpace(A));
	q:=Size(fld);
	A1:=SubspaceOfContractions(1,A,2,3,3);
	R1:=RankDistribution(A1,3,3);
	if q=2 then
		if R1=[[1,1]] then
			orbit:=1;
			representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]],[[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]]];
			return([orbit, representative]);
		else if R1=[[2,1]] then
			orbit:=2;
			representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]],[[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]]];
			return([orbit, representative]);
		else if R1=[[3,1]] then
			orbit:=3;
			representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),Z(q)^0]],[[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]]];
			return([orbit, representative]);
		else if R1=[[1,q+1]] then
			orbit:=4;
			representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]],[[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]]];
			return([orbit, representative]);
		else if R1=[[1,2],[2,q-1]] then
			orbit:=5;
			representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]],[[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]]];
			return([orbit, representative]);
		else if R1=[[1,1],[2,1],[3,q-1]] then
			orbit:=8;
			representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]],[[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),Z(q)^0]]];
			return([orbit, representative]);
		else if R1=[[1,1],[3,q]] then
			orbit:=9;
			representative:=[[[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)],[Z(q)^0,0*Z(q),0*Z(q)]],[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),Z(q)^0]]];
			return([orbit, representative]);
		else if R1=[[2,2],[3,q-1]] then
			orbit:=13;
			representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]],[[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),Z(q)^0]]];
			return([orbit, representative]);

		else if R1=[[3,q+1]] then
			orbit:=17;
			K:=GF(q);
			F:=GF(q^3);
			w:=PrimitiveElement(F);
			Mw:=CompanionMat(MinimalPolynomial(K,w));
			M1:=One(K^[3,3]);
			x:=VectorSpaceToElement(PG(8,q),Flat(M1));
			y:=VectorSpaceToElement(PG(8,q),Flat(Mw));
			representative:=[M1,Mw];
			return([orbit, representative]);
		else
			if R1=[[1,1],[2,q]]  then
				A2:=SubspaceOfContractions(2,A,2,3,3);
				A3:=SubspaceOfContractions(3,A,2,3,3);
				R2:=RankDistribution(A2,2,3);
				R3:=RankDistribution(A3,2,3);
				if R2=[[1,1],[2,q]] then
					if R3=[[1,1],[2,q]] then
						orbit:=6;
						representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]],[[0*Z(q),Z(q)^0,0*Z(q)],[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]]];
					else orbit:=7;
						representative:=[[[0*Z(q),0*Z(q),Z(q)^0],[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]],[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]]];
					fi;
				else
					orbit:=7;
					representative:=[[[0*Z(q),0*Z(q),Z(q)^0],[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]],[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]]];
				fi;
				return([orbit, representative]);

			else
				if R1=[[2,q+1]] then
					A2:=SubspaceOfContractions(2,A,2,3,3);
					A3:=SubspaceOfContractions(3,A,2,3,3);
					R2:=RankDistribution(A2,2,3);
					R3:=RankDistribution(A3,2,3);
					if R2=[[2,q+1]] then
						if R3=[[2,q+1]] then orbit:=10;
							representative:=RepO10even(q);
						else orbit:=11;
							representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]],[[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),Z(q)^0],[0*Z(q),0*Z(q),0*Z(q)]]];
						fi;
					else

						if R3=[[2,q+1]] then
							orbit:=11;
							representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]],[[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),Z(q)^0],[0*Z(q),0*Z(q),0*Z(q)]]];
						else if ProjectiveDimension(Span(Rank1PtsOftheContractionSubspace(A2,2,3)))=2 then
							orbit:=14;
							representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]],[[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),Z(q)^0]]];
						else	 orbit:=12;
							representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]],[[0*Z(q),0*Z(q),Z(q)^0],[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)]]];
						fi;
					fi;
				fi;
				return([orbit, representative]);
			else if R1=[[2,1],[3,q]] then
				x1:=First(Points(A1),p->RankOfPoint(p,3,3)=3);
				#x1:=Rank3PtsOftheContractionSubspace(A1,3,3)[1];
				x2:=First(Points(A1),p->RankOfPoint(p,3,3)=2);
				#x2:=Rank2PtsOftheContractionSubspace(A1,3,3)[1];
				M2:=MatrixOfPoint(x2,3,3);
				B:=Basis(Subspace(GF(q)^3,M2));
				v1:=List([1..3],i->Coefficients(B,M2[i])[1]);
				v2:=List([1..3],i->Coefficients(B,M2[i])[2]);
				t:=TransposedMat([v1])*[BasisVectors(B)[1]]; #y=t
				z:=TransposedMat([v2])*[BasisVectors(B)[2]];
				y1:=VectorSpaceToElement(PG(2,q),v1);
				y2:=VectorSpaceToElement(PG(2,q),BasisVectors(B)[1]);
				z1:=VectorSpaceToElement(PG(2,q),v2);
				z2:=VectorSpaceToElement(PG(2,q),BasisVectors(B)[2]);
				sv:=SegreVariety([PG(2,q),PG(2,q)]);
				sm:=SegreMap(sv);
				p1:=Cartesian(Points(y1),Points(y2));
				Ip1:=ImagesSet(sm,p1)[1];
				p2:=Cartesian(Points(y1),Points(z2));
				Ip2:=ImagesSet(sm,p2)[1];
				p3:=Cartesian(Points(z1),Points(y2));
				Ip3:=ImagesSet(sm,p3)[1];
				p4:=Cartesian(Points(z1),Points(z2));
				Ip4:=ImagesSet(sm,p4)[1];
				V:=Span([Ip1,Ip2,Ip3,Ip4]);
				U:=Span(V,x1);
				#if First(Points(sv),x->x in U and not x in V)=fail then
				if not ForAny(Points(U),x->not x in V and RankOfPoint(x,3,3)=1) then
					orbit:=16;
					representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),Z(q)^0]],[[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),Z(q)^0],[0*Z(q),0*Z(q),0*Z(q)]]];
				else
					orbit:=15;
					representative:=RepO15even(q);
				fi;
				return([orbit, representative]);


			fi;
		fi;
	fi;
fi;
fi;
fi;
fi;
fi;
fi;
fi;
fi;
fi;



else

	if R1=[[1,1]] then
		orbit:=1;
		representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]],[[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]]];
		return([orbit, representative]);
	else if R1=[[2,1]] then
		orbit:=2;
		representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]],[[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]]];
		return([orbit, representative]);
	else if R1=[[3,1]] then
		orbit:=3;
		representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),Z(q)^0]],[[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]]];
		return([orbit, representative]);
	else if R1=[[1,q+1]] then
		orbit:=4;
		representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]],[[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]]];
		return([orbit, representative]);
	else if R1=[[1,2],[2,q-1]] then
		orbit:=5;
		representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]],[[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]]];
		return([orbit, representative]);
	else if R1=[[1,1],[2,1],[3,q-1]] then
		orbit:=8;
		representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]],[[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),Z(q)^0]]];
		return([orbit, representative]);
	else if R1=[[1,1],[3,q]] then
		orbit:=9;
		representative:=[[[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)],[Z(q)^0,0*Z(q),0*Z(q)]],[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),Z(q)^0]]];
		return([orbit, representative]);
	else if R1=[[2,2],[3,q-1]] then
		orbit:=13;
		representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]],[[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),Z(q)^0]]];
		return([orbit, representative]);

	else if R1=[[3,q+1]] then
		orbit:=17;
		K:=GF(q);
		F:=GF(q^3);
		w:=PrimitiveElement(F);
		Mw:=CompanionMat(MinimalPolynomial(K,w));
		M1:=One(K^[3,3]);
		x:=VectorSpaceToElement(PG(8,q),Flat(M1));
		y:=VectorSpaceToElement(PG(8,q),Flat(Mw));
		representative:=[M1,Mw];
		return([orbit, representative]);
	else
		if R1=[[1,1],[2,q]]  then
			A2:=SubspaceOfContractions(2,A,2,3,3);
			A3:=SubspaceOfContractions(3,A,2,3,3);
			R2:=RankDistribution(A2,2,3);
			R3:=RankDistribution(A3,2,3);
			if R2=[[1,1],[2,q]] then
				if R3=[[1,1],[2,q]] then
					orbit:=6;
					representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]],[[0*Z(q),Z(q)^0,0*Z(q)],[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]]];
				else orbit:=7;
					representative:=[[[0*Z(q),0*Z(q),Z(q)^0],[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]],[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]]];
				fi;
			else
				orbit:=7;
				representative:=[[[0*Z(q),0*Z(q),Z(q)^0],[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]],[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]]];
			fi;
			return([orbit, representative]);
		else if R1=[[2,3],[3,q-2]] then
			orbit:=14;
			representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]],[[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),Z(q)^0]]];
			return([orbit, representative]);
		else
			if R1=[[2,q+1]] then
				A2:=SubspaceOfContractions(2,A,2,3,3);
				A3:=SubspaceOfContractions(3,A,2,3,3);
				R2:=RankDistribution(A2,2,3);
				R3:=RankDistribution(A3,2,3);
				if R2=[[2,q+1]] then
					if R3=[[2,q+1]] then orbit:=10;
						if 2 in PrimeDivisors(q) then
							representative:=RepO10even(q);
						else
							representative:=RepO10odd(q);
						fi;
					else orbit:=11;
						representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]],[[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),Z(q)^0],[0*Z(q),0*Z(q),0*Z(q)]]];
					fi;
				else
					if R3=[[2,q+1]] then
						orbit:=11;
						representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]],[[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),Z(q)^0],[0*Z(q),0*Z(q),0*Z(q)]]];
					else orbit:=12;
						representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),0*Z(q)]],[[0*Z(q),0*Z(q),Z(q)^0],[0*Z(q),0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)]]];
					fi;
				fi;
				return([orbit, representative]);
			else if R1=[[2,1],[3,q]] then
				x1:=First(Points(A1),p->RankOfPoint(p,3,3)=3);
				#x1:=Rank3PtsOftheContractionSubspace(A1,3,3)[1];
				x2:=First(Points(A1),p->RankOfPoint(p,3,3)=2);
				#x2:=Rank2PtsOftheContractionSubspace(A1,3,3)[1];
				M2:=MatrixOfPoint(x2,3,3);
				B:=Basis(Subspace(GF(q)^3,M2));
				v1:=List([1..3],i->Coefficients(B,M2[i])[1]);
				v2:=List([1..3],i->Coefficients(B,M2[i])[2]);
				t:=TransposedMat([v1])*[BasisVectors(B)[1]]; #y=t
				z:=TransposedMat([v2])*[BasisVectors(B)[2]];
				y1:=VectorSpaceToElement(PG(2,q),v1);
				y2:=VectorSpaceToElement(PG(2,q),BasisVectors(B)[1]);
				z1:=VectorSpaceToElement(PG(2,q),v2);
				z2:=VectorSpaceToElement(PG(2,q),BasisVectors(B)[2]);
				sv:=SegreVariety([PG(2,q),PG(2,q)]);
				sm:=SegreMap(sv);
				p1:=Cartesian(Points(y1),Points(y2));
				Ip1:=ImagesSet(sm,p1)[1];
				p2:=Cartesian(Points(y1),Points(z2));
				Ip2:=ImagesSet(sm,p2)[1];
				p3:=Cartesian(Points(z1),Points(y2));
				Ip3:=ImagesSet(sm,p3)[1];
				p4:=Cartesian(Points(z1),Points(z2));
				Ip4:=ImagesSet(sm,p4)[1];
				V:=Span([Ip1,Ip2,Ip3,Ip4]);
				U:=Span(V,x1);
				 if not ForAny(Points(U),x->not x in V and RankOfPoint(x,3,3)=1) then
				#if First(Points(sv),x->x in U and not x in V)=fail then
					orbit:=16;
					representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),Z(q)^0]],[[0*Z(q),Z(q)^0,0*Z(q)],[0*Z(q),0*Z(q),Z(q)^0],[0*Z(q),0*Z(q),0*Z(q)]]];
				else
					orbit:=15;
					if 2 in PrimeDivisors(q) then
						representative:=RepO15even(q);
					else
						representative:=RepO15odd(q);
					fi;

				fi;
				return([orbit, representative]);


			fi;
		fi;
	fi;
fi;
fi;
fi;
fi;
fi;
fi;
fi;
fi;
fi;
fi;

fi;
end);


InstallGlobalFunction(RankOfTensor,function(A)
	local fld, q, O;
	fld:=BaseField(AmbientSpace(A));
	q:=Size(fld);
	O:=OrbitOfTensor(A)[1];
	if O=1 then
		return(1);
	else if O=2 then
		return(2);
	else if O=3 then
		return(3);
	else if O=4 then
		return(2);
	else if O=5 then
		return(2);
	else if O=6 then
		return(3);
	else if O=7 then
		return(3);
	else if O=8 then
		return(3);
	else if O=9 then
		return(4);
	else if O=10 then
		return(3);
	else if O=11 then
		return(3);
	else if O=12 then
		return(4);
	else if O=13 then
		return(4);
	else if O=14 then
		return(3);
	else if O=15 then
		return(4);
	else if O=16 then
		return(3);
	else if O=17 then
		if q=2 then
			return(5);
		else
			return(4);
		fi;
	fi;
fi;
fi;
fi;
fi;
fi;
fi;
fi;
fi;
fi;
fi;
fi;
fi;
fi;
fi;
fi;
fi;
end);
