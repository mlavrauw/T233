
#############################################################################
##
##  RankOfTensor.gi              T233 package
##                                                            Nour Alnajjarine
##                                                            Michel Lavrauw
##                                                        
##
##  Copyright 2019	Sabanci University
##                  
##  Implementation stuff for rank of tensors in PG(F^2⊗F^3⊗F^3).
##
#############################################################################

# Let V:= F^2⊗F^3⊗F^3 where F denotes the finite field of order q, and consider the group G := GL(2, q) × GL(3, q) × GL(3, q) as a subgroup
# of GL(V) stabilising the set of fundamental tensors (i.e. tensors of rank one) in V.
# Note that, the set of fundamental tensors in PG(V) is projectively equivalent to the Segre variety:
# SegreVariety([PG(1,q),PG(2,q),PG(2,q)].
# The OrbitOfTensor function takes an arbitrary tensor in PG(V) and returns its orbit under the natural action of G on V and a representative of this orbit.
# The RankOfTensor function takes a tensor in PG(V) and returns its rank by using the classification of G-orbits of tensor in PG(V) mentioned in the 
# OrbitOfTensor function.
 
InstallGlobalFunction(RankOfTensor,function(A)
# We should turn this into a method which can also be applied to a tensor represented in a different way.
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
	                                    	return(4);
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
