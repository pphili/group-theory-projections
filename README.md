# group-theory-projections
Numerical application of the group theoretic projection operators for the T_d group. 

We use the tetrahedral single (double) group projection operators to projected basis states of the \Gamma_4 (\Gamma_8) representation onto each other. 

For the Mathematica version of this code, first run Td_realspace_rot.nb. This initializes all of the matrices needed for the projections.
The Check.nb file has a default test data set that is projected accurately onto the basis partners. If you would like to use the projectors on your own data, load them into the Mathematica file under the variable name "T" and run the cells under the heading "Discrete Test".

```python
T = Table[{x,y,z}, {Fup(x,y,z), Fdown(x,y,z)}]
```
