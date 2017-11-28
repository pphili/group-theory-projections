# group-theory-projections
Numerical application of the group theoretic projection operators as defined in 'Dresselhaus, M. S., Dresselhaus, G., & Jorio, A. (2007). Group theory: application to the physics of condensed matter. Springer Science & Business Media, Chapter 4'. Both versions of the code (Mathematica and Python) work in the same way and are applied to the example of the \Gamma_8 representation of the tetrahedral double group. The python version is easier to apply to any representation of any group, while the Mathematica version should only be used for visualization of states sampled on a small number of points (small datasets) because of its large memory requirements compared to python. 

The input data is the same for both versions of the code. The expected input is a list of points (x,y,z coordinates) followed by the spin-up and spin-down components of the state at that point. Each new point should be on a new line. The sampling points must be closed under the symmetry operations of the group. In other words, any symmetry operation of the group applied to any point in you input data yields another point in your input data. This requirement is necessary for the projection operators to yield accurate states. 

The file exampledata_npts21.OUT is an example input. It contains a description of a valence-band state of GaAs, in the cubic unit cell. This wavefunction is sampled at 21 equally spaced points along each coordinate, with points sampled along the boundary of the cubic unit cell. The reason we chose the cubic unit cell and sampled points on the entire boundary of the cubic unit cell is to conform to the requirement that the sampling points must be closed under the symmetries of the group (see above paragraph). We will use the projection operators to output the four basis states of the \Gamma_8 representation of the tetrahedral double group which correspond to the two light-hole states and the two heavy-hole states.  

## Mathematica Version
This version is recommended only for testing and visualizing the projections for small data files. The notbooks are constructed specifically for projecting onto the basis states of the \Gamma_8 representation of the tetrahedral double group, but can be generalized to other representations of the same or any other group. 

First run Td_realspace_rot.nb. This initializes all of the matrices needed for the projections. In this file, `L` is the order  (number of elements) of the group under consideration, the matrices `R[i]` represent the transformation of the coordinates under the application of the `i`th group element, `rs[i]` represent the corresponding transformation of the spinors, and `\[CapitalRho][i]` are the corresponding representations. These variables can be changed to store matrices corresponding to a different representation or a different group and thus implement other projection operators. 

### Demonstration
#### Single Group 
We consruct the projection operators for the \Gamma_4 representation of the tetrahedral group. The CheckSG.nb file has a default input test state that is one of the basis states of this representation. It is loaded into the variable `T`

```python
T = Table[{x,y,z}, F(x,y,z), {x, x_min, x_max}, {y, y_min, y_max}, {z, z_min, z_max}],
```
where `F(x,y,z)` is the value of the wavefunction at the point `(x,y,z)`. It is not a spinor becasue we are considering single groups here. 
The function `AR[m, n, r, \[Theta], \[Phi]]`, will be the projected function, where `m` and `n` variables indicate the state you are projecting from and the state you are projecting to, respectively and `r`, `\[Theta]`, `\[Phi]` are the spherical coordinates.

The plots at the end of the notebook attempt to demonstrate the agreement of the projected states with the expected projections.

#### Double group
We consruct the projection operators for the \Gamma_8 representation of the tetrahedral double group.

The CheckDG.nb file has a default test data set that is projected accurately onto the basis partners. It is loaded into the variable `T` which is different now compared to the single group case because there is a spinor associated to each point in space

```python
T = Table[{x,y,z}, {Fup(x,y,z), Fdown(x,y,z)}, {x, x_min, x_max}, {y, y_min, y_max}, {z, z_min, z_max}],
```
where `Fup(x,y,z)` and `Fdown(x,y,z)` are the spin-up and spin-down components of the states respectively. Again, the plots at the end of the notebook are meant to demonstrate the agreement between the projected states and the expected output. 

### Example - Double Group

When you are ready to use the projection operators on your data, input the path of your input file into the inic=tialization cell of the `pathToInput` variable. Again, remember to update the Td_realspace_rot.nb notebook to put the appropriate matrices for the projection operator you are interested in.  Run the cells under the heading "Projections". `Proj[ii,jj,r]` contains the projected states. The `ii` and `jj` variables indicate the state you are projecting from and the state you are projecting to, respectively, and `r={x,y,z}` is the coordinate. The functions defined by `Proj[ii,jj,r]` are interpolated using Mathematica's `Interpolation` function. 

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
