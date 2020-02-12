

#############################################################################
##
##  OrbitOfTensor.gi         T233 package
##                                                            Nour Alnajjarine
##                                                            Michel Lavrauw

##
##  Copyright 2019	Sabanci University
##                
##
##  Implementation stuff for orbits of tensors in PG(GF(q)^2⊗GF(q)^3⊗GF(q)^3)
##
#############################################################################
# Let V:= F^2⊗F^3⊗F^3 where F denotes the finite field of order q, and consider the group G := GL(2, q) × GL(3, q) × GL(3, q) as a subgroup
# of GL(V) stabilising the set of fundamental tensors (i.e. tensors of rank one) in V.
# Note that, the set of fundamental tensors in PG(V) is projectively equivalent to the Segre variety:
# SegreVariety([PG(1,q),PG(2,q),PG(2,q)].
# The OrbitOfTensor function takes an arbitrary tensor in PG(V) and returns its orbit under the natural action of G on V and a representative of this orbit.



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