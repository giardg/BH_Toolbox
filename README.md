# PowerEquivalentModel

The power-equivalent model is a numerical package used to create equivalent permeability curves for nonlinear and hysteretic materials in the harmonic domain from magnetic measurements.

## Usage
There are two steps to follow for creating an equivalent permeability curve:
1) Run FEM_transitoire_fortran\SlabProblem_IO.m with the electric and magnetic properties of a given material. This code solves the 1-D time-transient slab diffusion problem with the scalar Preisach model to compute the eddy current and hysteresis losses as a function of x. The problem is solved in Fortran 90 using the finite element method (FEM). The losses curve are then stored in the folder. The current values used in our software represent magnetic measurements done on AISI4340 steel at Polytechnique Montreal. Corresponding loss curves were stored in Resultats_transitoire_fortran_aisi4340. A user could describe other materials with this model.
2) Run FEM_PEM_1D\main.m The losses created from 1) are used as inputs in the power-equivalent model. The equivalent permeability curve is stored in FEM_PEM_1D\Resultats_mu_aisi4340.

## Citation
Please refer to CITATION.cff to cite this repository
