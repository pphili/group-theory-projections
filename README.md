# group-theory-projections
Numerical application of the group theoretic projection operators as defined in 'Dresselhaus, M. S., Dresselhaus, G., & Jorio, A. (2007). Group theory: application to the physics of condensed matter. Springer Science & Business Media, Chapter 4'. Both versions of the code (Mathematica and Python) work in the same way, and are applied to the example of the \Gamma_8 representation of the tetrahedral double group. However, the python version is easier to apply to any representation of any group, while the Mathematica version should only be used for visualization of states sampled on a small number of points (small datasets). 

The input data is the same for both versions of the code. The expected input is a list of points (x,y,z coordinates) followed by the spin-up and spin-down components of the state at that point. The sampling points must be closed under the symmetry operations of the group. In other words, 
```math #yourmathlabel
a + b = c
```



Each new point should be on a newline. The file exampledata_npts21.OUT is an example input. It contains a description of a valence-band state of GaAs, in the cubic unit cell. This wavefunction is sampled at 21 equally spaced points along each coordinate. 

## Mathematica Version
This version is recommended only for testing and visulizing the projections for small data files.

First run Td_realspace_rot.nb. This initializes all of the matrices needed for the projections.

### Single Group 
We consruct the projection operators for the \Gamma_4 representation of the tetrahedral group. The CheckSG.nb file has a default input test state that is one of the basis states of this representation. Once your function is loaded into the notebook under the variable 

```python
T = Table[{x,y,z}, F(x,y,z), {x, x_min, x_max}, {y, y_min, y_max}, {z, z_min, z_max}],
```
the rest of the notebook can be run. The function `AR[m, n, r, \[Theta], \[Phi]]`, will be the projected function, where `m` and `n` variables indicate the state you are projecting from and the state you are projecting to, respectively and `r`, `\[Theta]`, `\[Phi]` are the spherical coordinates.

### Double group
We consruct the projection operators for the \Gamma_8 representation of the tetrahedral double group.

The CheckDG.nb file has a default test data set that is projected accurately onto the basis partners. If you would like to use the projectors on your own data, load them into the Mathematica file under the variable name `T`. `T` must be a table with each entry containing the position coordinates and the spinor associated to that point. 

```python
T = Table[{x,y,z}, {Fup(x,y,z), Fdown(x,y,z)}, {x, x_min, x_max}, {y, y_min, y_max}, {z, z_min, z_max}].
```
Important: Imported discrete function must be defined on a cubic grid centered at the origin. Furthermore the grid points must be equally spaced and symmetric (same number of grid points along the positive and negative directions).

When you are ready to use the projection operators, run the cells under the heading "Projections". The `m1` and `m2` variables indicate the state you are projecting from and the state you are projecting to, respectively. Finally `Ftab` represents a table  
```python
Ftab = Table[{x,y,z}, {PFup(x,y,z), PFdown(x,y,z)}, {x, x_min, x_max}, {y, y_min, y_max}, {z, z_min, z_max}]
```
where `PFup(x,y,z)` [`PFdown(x,y,z)`] is the spin up (down) component of the projected function. 
For visualization puposes the variable `view[x,y,z]` is an interpolation of the projected data that can be easily plotted

## Python Version
This version of the script can be used to output the each basis state of a representation of a group given an intial arbitrary linear combination of all basis states. This is done by implementing the group-theoretic projection operators defined in Dresselhaus, M. S., Dresselhaus, G., & Jorio, A. (2007). Group theory: application to the physics of condensed matter. Springer Science & Business Media, Chapter 4. There are two files required to perform the projections, projection.py and a file containing the matrices necessary to construct the projection operators which we call reps.py.  

### reps.py
The reps.py script should contain all the matrices required for the projections. All matrices should be stored as numpy arrays. `R` should be a list of numpy arrays determining the transformation of the coordinate system under the operations of the group, `wD12` should be a list of numpy arrays describing the transformation of a spinor under the operations of the group and finally, `G` should be the representation for which we want the basis states. You can also include a list of `labels` for each basis state in the representation. The default for the `labels` will be numbers starting from 0. The ordering of the lists matters as each group element should have the same index in each list. 

As an example, G8Td.py is the reps.py file that contains the appropriate matrices to perform the projections onto the four basis states of the \Gamma_8 representation (G8) of the tetrahedral (Td) double group. In this case, the labels represent heavy-hole and light-hole states. 
```
labels = ['HHd', 'LHd', 'LHu', HHu']
```

### projection.py
The projection.py script contains the functions necessary to perform the projections. It imports reps.py. 

### Usage
Input file: The input file should be a list of points on a cubic grid, that are equally spaced in all three dimensions, centered at the origin, (0,0,0). Furthermore, the endpoints should be included. Each line of the input file should also contain the two complex numbers identifying the amplitudes of each component of the spinor at that point. This file should describe the state that will be projected onto the basis states of the desired representation. A typical line in the input file should look like:

```
x-coordinate y-coordinate z-coordinate (spin-up amplitude) (spin-down amplitude) 
```

#### Default run:

```
python projection.py 
```
You will be prompted to enter the filename describing the state you wish to apply the projection operators on. Here the file is called `state_to_project`

```
filename: 'state_to_project'
```

Finally, you will be prompted to enter the number of cores for which you would like to run this process on. In this example we pick 4

```
number of cores: 4
```

#### Default output:
The output will be a series of files each named `'proj'+labels[i]+'.OUT'` for the different elements of `labels`, each one representing a projected state written out in the same convention as the input state. 
