# group-theory-projections
Numerical application of the group theoretic projection operators as defined in 'Dresselhaus, M. S., Dresselhaus, G., & Jorio, A. (2007). Group theory: application to the physics of condensed matter. Springer Science & Business Media, Chapter 4'. Both versions of the code (Mathematica and Python) work in the same way and are applied to the example of the \Gamma_8 representation of the tetrahedral double group. The python version is easier to apply to any representation of any group, while the Mathematica version should only be used for visualization of states sampled on a small number of points (small datasets) because of its large memory requirements compared to python. 

The input data is the same for both versions of the code. The expected input is a list of points (x,y,z coordinates) followed by the spin-up and spin-down components of the state at that point. Each new point should be on a new line. The sampling points must be closed under the symmetry operations of the group. In other words, any symmetry operation of the group applied to any point in you input data yields another point in your input data. This requirement is necessary for the projection operators to yield accurate states. 

The file exampledata_npts21.OUT is an example input. It contains a description of a valence-band state of GaAs, in the cubic unit cell. This wavefunction is sampled at 21 equally spaced points along each coordinate, with points sampled along the boundary of the cubic unit cell. The reason we chose the cubic unit cell and sampled points on the entire boundary of the cubic unit cell is to conform to the requirement that the sampling points must be closed under the symmetries of the group (see above paragraph). We will use the projection operators to output the four basis states of the \Gamma_8 representation of the tetrahedral double group which correspond to the two light-hole states and the two heavy-hole states.  

For the projection operators applied to angular momentum eigenstates, see the `group-theory-projections_sh` directory

## Mathematica Version
This version is recommended only for testing and visualizing the projections for small data files. The notbooks are constructed specifically for projecting onto the basis states of the \Gamma_8 representation of the tetrahedral double group, but can be generalized to other representations of the same or any other group. 

First run Td_realspace_rot.nb. This initializes all of the matrices needed for the projections. In this file, `L` is the order  (number of elements) of the group under consideration, the matrices `R[i]` represent the transformation of the coordinates under the application of the `i`th group element, `rs[i]` represent the corresponding transformation of the spinors, and `\[CapitalRho][i]` are the corresponding representations. These variables can be changed to store matrices corresponding to a different representation or a different group and thus implement other projection operators. 

### Demonstration
#### Single Group 
We consruct the projection operators for the \Gamma_4 representation of the tetrahedral group. The CheckSG.nb file has a default input test state that is one of the basis states of this representation. It is loaded into the variable `T`

```python
T = Table[{x,y,z}, F(x,y,z), {x, x_min, x_max}, {y, y_min, y_max}, {z, z_min, z_max}]
```
where `F(x,y,z)` is the value of the wavefunction at the point `(x,y,z)`. It is not a spinor becasue we are considering single groups here. 
The function `AR[m, n, r, \[Theta], \[Phi]]`, will be the projected function, where `m` and `n` variables indicate the state you are projecting from and the state you are projecting to, respectively and `r`, `\[Theta]`, `\[Phi]` are the spherical coordinates.

The plots at the end of the notebook attempt to demonstrate the agreement of the projected states with the expected projections.

#### Double group
We consruct the projection operators for the \Gamma_8 representation of the tetrahedral double group.

The CheckDG.nb file has a default test data set that is projected accurately onto the basis partners. It is loaded into the variable `T` which is different now compared to the single group case because there is a spinor associated to each point in space

```python
T = Table[{x,y,z}, {Fup(x,y,z), Fdown(x,y,z)}, {x, x_min, x_max}, {y, y_min, y_max}, {z, z_min, z_max}]
```
where `Fup(x,y,z)` and `Fdown(x,y,z)` are the spin-up and spin-down components of the states respectively. Again, the plots at the end of the notebook are meant to demonstrate the agreement between the projected states and the expected output. 

### Example - Double Group

When you are ready to use the projection operators on your data, input the path of your input file into the inic=tialization cell of the `pathToInput` variable. Again, remember to update the Td_realspace_rot.nb notebook to put the appropriate matrices for the projection operator you are interested in.  Run the cells under the heading "Projections". `Proj[ii,jj,r]` contains the projected states. The `ii` and `jj` variables indicate the state you are projecting from and the state you are projecting to, respectively, and `r={x,y,z}` is the coordinate. The functions defined by `Proj[ii,jj,r]` are interpolated using Mathematica's `Interpolation` function. 

## Python Version
Similarly to the Mathematica version of the script, there are two files required to perform the projections, projection.py and a file containing the matrices necessary to construct the projection operators which we call reps.py.  

#### Requirements 
* python 2.7+
* transforms3d (https://pypi.python.org/pypi/transforms3d)
* multiprocessing (https://docs.python.org/2/library/multiprocessing.html)
* functools (https://docs.python.org/2/library/functools.html)

### reps.py
The reps.py script should contain all the matrices required for the projections. All matrices should be stored as numpy arrays. `R` should be a list of numpy arrays determining the transformation of the coordinate system under the operations of the group, `wD12` should be a list of numpy arrays describing the transformation of a spinor under the operations of the group and finally, `G` should be the representation for which we want the basis states. You can also include a list of `labels` for each basis state in the representation. The default for the `labels` will be numbers starting from 0. The ordering of the lists matters as each group element should have the same index in each list. 

As an example, G8Td.py is the reps.py file that contains the appropriate matrices to perform the projections onto the four basis states of the \Gamma_8 representation (G8) of the tetrahedral (Td) double group. In this case, the labels represent heavy-hole and light-hole states. 
```
labels = ['HHd', 'LHd', 'LHu', HHu']
```

### projection.py
The projection.py script contains the functions necessary to perform the projections. It imports reps.py. 

#### Default run:

```
python projection.py 
```
You will be prompted to enter the filename describing the state you wish to apply the projection operators on. The default state (no input) is the state defined in the file exampledata_npts21.OUT. Here the file is called `state_to_project`

```
filename: 'state_to_project'
```

Finally, you will be prompted to enter the number of cores for which you would like to run this process on. The default is 1. In this example we pick 4

```
number of cores: 4
```

#### Default output:
The output will be a series of files each named `'proj'+labels[i]+'.OUT'` for the different elements of `labels`, each one representing a projected state written out in the same convention as the input state. In the default case the output will be four files each representing a different basis state of the \Gamma_8 representation of the tetrahedral double group. More specifically, the four output states will be the heavy-hole and light-hole states at the \Gamma point of GaAs. 

### matelems.py
This module can be used to compute the matrix elements of position, momentum, spin and angular momentum between two states defined on a grid in real space. 

`psi(coords , wf)` 
creates a `psi` object with coordinates `coords` and spinor at each coordinate determined by `wf`. `coords` is a list of lists, while `wf` is a list of numpy arrays. We take here the coordinates to be uniformly spaced. `dV` contains the volume element determined from the coordinates spacing to be used for numerical integrations. 

#### Methods
`psi.sameDomain(psi2)`
compares the coords attribute of `psi` and `psi2` and returns `True` if they are the same and `False` if they differ.

`psi.position(psi2)`
returns a tuple corresponding to the matrix elements of (position) x, y and z, determined from a Riemann sum over all the coordinates if `psi.sameDomain(psi2) = True`.

`psi.p(psi2)`
returns a tuple corresponding to the matrix elements of (momentum) p_x, p_y and p_z, determined from a Riemann sum over all the coordinates if `psi.sameDomain(psi2) = True`.

`psi.spin(psi2)`
returns a tuple corresponding to the matrix elements of (spin) S_x, S_y and S_z, determined from a Riemann sum over all the coordinates if `psi.sameDomain(psi2) = True`.

`psi.L(psi2)`
returns a tuple corresponding to the matrix elements of (angular momentum) L_x, L_y and L_z, determined from a Riemann sum over all the coordinates if `psi.sameDomain(psi2) = True`.

`psiFromFile( fileName, r = 6 )`
returns a `psi` object from a file with name `filename`. `r` determines to which decimal the coordinates should be rounded. 


