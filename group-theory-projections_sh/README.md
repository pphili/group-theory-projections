# group-theory-projections
Numerical application of the group theoretic projection operators as defined in 'Dresselhaus, M. S., Dresselhaus, G., & Jorio, A. (2007). Group theory: application to the physics of condensed matter. Springer Science & Business Media, Chapter 4'. The projection operators are applied to eigenstates of angular momentum. An arbitrary state |J, m, l, s>, is defined by 4 quantum numbers: J is the total angular momentum, m is the angular momentum along the axis of quantization, l is the orbital angular momentum and s is the spin. For the purposes of this script, we fix s=1/2.

##Inputs
When running projection_sh.py, there needs to be another file in the same directory that contains the Euler angles of all the transformations of the group as well as if the transformation contains and inversion (-1) or not (1). Furthermore it must contain the representation for which the projection operators will be defined. See G8Td.py as an example of the file for the Gamma 8 representation of the Td double group. 

Once projection_sh.py is run, it will prompt the user for values of J, m and l as well as the basis state the user would like to project to and from. 
