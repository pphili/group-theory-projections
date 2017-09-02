# group-theory-projections
Numerical application of the group theoretic projection operators for the T_d group. We use the tetrahedral single (double) group projection operators to projected basis states of the \Gamma_4 (\Gamma_8) representation onto each other. 

##Mathematica Version

First run Td_realspace_rot.nb. This initializes all of the matrices needed for the projections.
The Check.nb file has a default test data set that is projected accurately onto the basis partners. If you would like to use the projectors on your own data, load them into the Mathematica file under the variable name `T`. `T` must be a table with each entry containing the position coordinates and the spinor associated to that point. 

```python
T = Table[{x,y,z}, {Fup(x,y,z), Fdown(x,y,z)}, {x, x_min, x_max}, {y, y_min, y_max}, {z, z_min, z_max}]
```

When you are ready to use the projection operators, run the cells under the heading "Projections". The `m1` and `m2` variables indicate the state you are projecting from and the state you are projecting to, respectively. Finally `Ftab` represents a table  
```python
Ftab = Table[{x,y,z}, {PFup(x,y,z), PFdown(x,y,z)}, {x, x_min, x_max}, {y, y_min, y_max}, {z, z_min, z_max}]
```
where `PFup(x,y,z)` [`PFdown(x,y,z)`] is the spin up (down) component of the projected function. 
For visualization puposes the variable `view[x,y,z]` is an interpolation of the projected data that can be easily plotted

##Python Version

Coming Soon
