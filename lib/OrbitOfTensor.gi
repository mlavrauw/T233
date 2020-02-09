

#############################################################################
##
##  OrbitOfTensor.gd         T233 package
##                                                            Nour Alnajjarine
##                                                            Michel Lavrauw

##
##  Copyright 2019	Sabanci University
##                
##
##  Implementation stuff for orbits of tensors in PG(GF(q)^2⊗GF(q)^3⊗GF(q)^3)
##
#############################################################################


InstallGlobalFunction(OrbitOfTensor,function(A)
local fld, q, A1, A2, A3, R1, R2, R3, orbit, rank, representative,x1, x2, M2, f1, f2, sv, sm, p1, Ip1,
 p2, Ip2, p3, Ip3, p4, Ip4, V, U, x, y, M1, Mw, w, K, F;
fld:=BaseField(AmbientSpace(A));
q:=Size(fld);
A1:=SubspaceOfContractions(1,A,2,3,3);
A2:=SubspaceOfContractions(2,A,2,3,3);
A3:=SubspaceOfContractions(3,A,2,3,3);
R1:=RankDistribution(A1,3,3);
R2:=RankDistribution(A2,2,3);
R3:=RankDistribution(A3,2,3);
if R1=[[1,1]] then
  orbit:=1; rank:=1;
representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],
                          [0*Z(q),0*Z(q),0*Z(q)],
                          [0*Z(q),0*Z(q),0*Z(q)]],
                        [[0*Z(q),0*Z(q),0*Z(q)], 
                         [0*Z(q),0*Z(q),0*Z(q)],
                         [0*Z(q),0*Z(q),0*Z(q)]]];
return([orbit,rank, CubicalArrayFromPointInTensorProductSpace(A,2,3,3),
          representative]);
else
if R1=[[2,1]] then
 orbit:=2; rank:=2;
representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],
                          [0*Z(q),Z(q)^0,0*Z(q)],
                          [0*Z(q),0*Z(q),0*Z(q)]],
                         [[0*Z(q),0*Z(q),0*Z(q)],
                          [0*Z(q),0*Z(q),0*Z(q)],
                          [0*Z(q),0*Z(q),0*Z(q)]]];
return([orbit,rank, CubicalArrayFromPointInTensorProductSpace(A,2,3,3),
          representative]);
else
if R1=[[3,1]] then
 orbit:=3; rank:=3;
representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],
                          [0*Z(q),Z(q)^0,0*Z(q)],
                          [0*Z(q),0*Z(q),Z(q)^0]],
                         [[0*Z(q),0*Z(q),0*Z(q)],
                          [0*Z(q),0*Z(q),0*Z(q)],
                          [0*Z(q),0*Z(q),0*Z(q)]]];
return([orbit,rank, CubicalArrayFromPointInTensorProductSpace(A,2,3,3),
          representative]);
else
if R1=[[1,q+1]] then
 orbit:=4; rank:=2;
representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],
                          [0*Z(q),0*Z(q),0*Z(q)],
                          [0*Z(q),0*Z(q),0*Z(q)]],
                         [[0*Z(q),Z(q)^0,0*Z(q)],
                          [0*Z(q),0*Z(q),0*Z(q)],
                          [0*Z(q),0*Z(q),0*Z(q)]]];
return([orbit,rank, CubicalArrayFromPointInTensorProductSpace(A,2,3,3),
          representative]);
else
if R1=[[1,2],[2,q-1]] then
 orbit:=5; rank:=2;
representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],
                          [0*Z(q),0*Z(q),0*Z(q)],
                          [0*Z(q),0*Z(q),0*Z(q)]],
                         [[0*Z(q),0*Z(q),0*Z(q)],
                          [0*Z(q),Z(q)^0,0*Z(q)],
                          [0*Z(q),0*Z(q),0*Z(q)]]];
return([orbit,rank, CubicalArrayFromPointInTensorProductSpace(A,2,3,3),
          representative]);
else
if R1=[[1,1],[2,1],[3,q-1]] then
 orbit:=8; rank:=3;
representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],
                          [0*Z(q),0*Z(q),0*Z(q)],
                          [0*Z(q),0*Z(q),0*Z(q)]],
                         [[0*Z(q),0*Z(q),0*Z(q)],
                          [0*Z(q),Z(q)^0,0*Z(q)],
                          [0*Z(q),0*Z(q),Z(q)^0]]];
return([orbit,rank, CubicalArrayFromPointInTensorProductSpace(A,2,3,3),
          representative]);
else 
if R1=[[1,1],[3,q]] then
 orbit:=9; rank:=4;
representative:=[[[0*Z(q),0*Z(q),0*Z(q)],
                          [0*Z(q),0*Z(q),0*Z(q)],
                          [Z(q)^0,0*Z(q),0*Z(q)]],
                         [[Z(q)^0,0*Z(q),0*Z(q)],
                          [0*Z(q),Z(q)^0,0*Z(q)],
                          [0*Z(q),0*Z(q),Z(q)^0]]];
return([orbit,rank, CubicalArrayFromPointInTensorProductSpace(A,2,3,3),
          representative]);
else
if R1=[[2,2],[3,q-1]] then
 orbit:=13; rank:=4;
representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],
                          [0*Z(q),Z(q)^0,0*Z(q)],
                          [0*Z(q),0*Z(q),0*Z(q)]],
                         [[0*Z(q),Z(q)^0,0*Z(q)],
                          [0*Z(q),0*Z(q),0*Z(q)],
                          [0*Z(q),0*Z(q),Z(q)^0]]];
return([orbit,rank, CubicalArrayFromPointInTensorProductSpace(A,2,3,3),
          representative]);
else
if R1=[[2,3],[3,q-2]] then
 orbit:=14; rank:=3;
representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],
                          [0*Z(q),Z(q)^0,0*Z(q)],
                          [0*Z(q),0*Z(q),0*Z(q)]],
                         [[0*Z(q),0*Z(q),0*Z(q)],
                          [0*Z(q),Z(q)^0,0*Z(q)],
                          [0*Z(q),0*Z(q),Z(q)^0]]];
return([orbit,rank, CubicalArrayFromPointInTensorProductSpace(A,2,3,3),
          representative]);
else 
if R1=[[3,q+1]] then
  orbit:=17;
  K:=GF(q);
  F:=GF(q^3);
  w:=PrimitiveElement(F);
  Mw:=CompanionMat(MinimalPolynomial(K,w));
  M1:=One(K^[3,3]);
  x:=VectorSpaceToElement(PG(8,q),Flat(M1));
  y:=VectorSpaceToElement(PG(8,q),Flat(Mw));
  representative:=[M1,Mw];
  if q=2 then rank:=5;
  else 
  rank:=4; 
  fi;
return([orbit,rank, CubicalArrayFromPointInTensorProductSpace(A,2,3,3),
          representative]);
else 
if R1=[[1,1],[2,q]]  then
 if R3=[[1,1],[2,q]] then
 orbit:=6; rank:=3;
 representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],
                           [0*Z(q),0*Z(q),0*Z(q)],
                           [0*Z(q),0*Z(q),0*Z(q)]],
                          [[0*Z(q),Z(q)^0,0*Z(q)],
                           [Z(q)^0,0*Z(q),0*Z(q)],
                           [0*Z(q),0*Z(q),0*Z(q)]]];
 else orbit:=7; rank:=3;
 representative:=[[[0*Z(q),0*Z(q),Z(q)^0],
                           [0*Z(q),0*Z(q),0*Z(q)],
                           [0*Z(q),0*Z(q),0*Z(q)]],
                          [[Z(q)^0,0*Z(q),0*Z(q)],
                           [0*Z(q),Z(q)^0,0*Z(q)],
                           [0*Z(q),0*Z(q),0*Z(q)]]];
 fi;
 return([orbit,rank, CubicalArrayFromPointInTensorProductSpace(A,2,3,3),
          representative]);
else
if R1=[[2,q+1]] then
 if R2=[[1,q+1],[2,q^2]] then
 orbit:=12; rank:=4;
 representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],
                           [0*Z(q),Z(q)^0,0*Z(q)],
                           [0*Z(q),0*Z(q),0*Z(q)]],
                          [[0*Z(q),0*Z(q),Z(q)^0],
                           [0*Z(q),0*Z(q),0*Z(q)],
                           [0*Z(q),Z(q)^0,0*Z(q)]]];
 else
  if R3=[[2,q+1]] then
  orbit:=10; rank:=3;
  representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],
                            [0*Z(q),Z(q)^0,0*Z(q)],
                            [0*Z(q),0*Z(q),0*Z(q)]],
                           [[0*Z(q),Z(q)^0,0*Z(q)],
                            [0*Z(q),0*Z(q),0*Z(q)],
                            [0*Z(q),0*Z(q),0*Z(q)]]];

  else
  orbit:=11; rank:=3;
  representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],
                            [0*Z(q),Z(q)^0,0*Z(q)],
                            [0*Z(q),0*Z(q),0*Z(q)]],
                           [[0*Z(q),Z(q)^0,0*Z(q)],
                            [0*Z(q),0*Z(q),Z(q)^0],
                            [0*Z(q),0*Z(q),0*Z(q)]]];
  fi;
 fi;
return([orbit,rank, CubicalArrayFromPointInTensorProductSpace(A,2,3,3),
          representative]);
else
x1:=Rank3PtsOftheContractionSubspace(A1,3,3)[1];
x2:=Rank2PtsOftheContractionSubspace(A1,3,3)[1];
M2:=TriangulizedMat(MatrixOfPoint(x2,3,3));
f1:=VectorSpaceToElement(PG(2,q),M2[1]);
f2:=VectorSpaceToElement(PG(2,q),M2[2]);
sv:=SegreVariety([PG(2,q),PG(2,q)]);
sm:=SegreMap(sv);
p1:=Cartesian(Points(f1),Points(f1));
Ip1:=ImagesSet(sm,p1)[1];
p2:=Cartesian(Points(f1),Points(f2));
Ip2:=ImagesSet(sm,p2)[1];
p3:=Cartesian(Points(f2),Points(f1));
Ip3:=ImagesSet(sm,p3)[1];
p4:=Cartesian(Points(f2),Points(f2));
Ip4:=ImagesSet(sm,p4)[1];
V:=Span([Ip1,Ip2,Ip3,Ip4]);
U:=Span(V,x1);
 if First(Points(sv),x->x in U and not x in V)=fail then
 orbit:=16; rank:=4;
 representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],
                           [0*Z(q),Z(q)^0,0*Z(q)],
                           [0*Z(q),0*Z(q),Z(q)^0]],
                          [[0*Z(q),Z(q)^0,0*Z(q)],
                           [0*Z(q),0*Z(q),Z(q)^0],
                           [0*Z(q),0*Z(q),0*Z(q)]]];
 else
 orbit:=15; rank:=4;
 representative:=[[[Z(q)^0,0*Z(q),0*Z(q)],
                           [0*Z(q),Z(q)^0,0*Z(q)],
                           [0*Z(q),0*Z(q),Z(q)^0]],
                          [[0*Z(q),Z(q)^0,0*Z(q)],
                           [0*Z(q),0*Z(q),0*Z(q)],
                           [0*Z(q),0*Z(q),0*Z(q)]]];
return([orbit,rank, CubicalArrayFromPointInTensorProductSpace(A,2,3,3),
          representative]);
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
