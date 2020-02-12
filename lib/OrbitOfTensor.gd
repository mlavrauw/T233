

#############################################################################
##
##  OrbitOfTensor.gd         T233 package
##                                                            Nour Alnajjarine
##                                                            Michel Lavrauw

##
##  Copyright 2019	Sabanci University
##                
##
##  Declaration stuff for orbits of tensors in PG(GF(q)^2⊗GF(q)^3⊗GF(q)^3)
##
#############################################################################

# Let V:= F^2⊗F^3⊗F^3 where F denotes the finite field of order q, and consider the group G := GL(2, q) × GL(3, q) × GL(3, q) as a subgroup
# of GL(V) stabilising the set of fundamental tensors (i.e. tensors of rank one) in V.
# Note that, the set of fundamental tensors in PG(V) is projectively equivalent to the Segre variety:
# SegreVariety([PG(1,q),PG(2,q),PG(2,q)].
# The OrbitOfTensor function takes an arbitrary tensor in PG(V) and returns its orbit under the natural action of G on V and a representative of this orbit.

DeclareGlobalFunction( "OrbitOfTensor" );

