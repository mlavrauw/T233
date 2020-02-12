
#############################################################################
##
## RepO10.gi              T233 package
##                                                            Nour Alnajjarine
##                                                            Michel Lavrauw
##                                                        
##
##  Copyright 2019	Sabanci University
##                  
##  Implementation stuff for representative of the orbit o10
##
#############################################################################


InstallGlobalFunction(RepO10odd,function(q)
	local eta,v1,v2,ext;
	eta:=PrimitiveElement(GF(q));
	v1:=[1,0,0,0,eta,0,0,0,0]*One(GF(q));
	v2:=[0,1,0,1,0,0,0,0,0]*One(GF(q));
	ext:=[MatrixOfPoint(VectorSpaceToElement(PG(8,q),v1),3,3),MatrixOfPoint(VectorSpaceToElement(PG(8,q),v2),3,3)];
	return ext;
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