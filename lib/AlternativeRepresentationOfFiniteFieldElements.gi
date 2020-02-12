
#############################################################################
##
## AlternativeRepresentationOfFiniteFieldElements.gi              T233 package
##                                                            Nour Alnajjarine
##                                                            Michel Lavrauw
##                                                        
##
##  Copyright 2019	Sabanci University
##                  
##  Implementation stuff for alternative representation of finite field elements
##
#############################################################################

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