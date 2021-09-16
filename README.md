# PowerEquivalentModel

The power-equivalent model is a numerical package used to create equivalent permeability curves for nonlinear and hysteretic materials in the harmonic domain from magnetic measurements.

## Usage
PEM_IO.m is the main script to run the Power-equivalent model. The input parameters are:
- freq: the frequency of the magnetic field used in the model
- L: the length of the 1-D slab domain
- rho: the isotropic resistivity of the material
- H0_list: a list of all the H0 values to apply at x = 0 (Dirichelet condition)
- mattype: the type of magnetic material to model
- output: the output folder to store the results

The variable "mattype" is a number from 1 to 7, each describing one permeability model:
1) Linear permeability (B = mu0.mu_r)
2) Arctangent anhysteretic curve (B = mu0.H + (2.bsat/pi).atan((0.5 .pi.mu0.(murmax-1)/bsat).H))
3) Hysteretic curve with 3n coeff Preisach model (B+ = Sum ai.atan((H+ci)/bi), B- = Sum ai.atan((H-ci)/bi))
4) Limit case for high field (B = mu0.H+Bsat if H > Hsat, B = mu0.H-Bsat si H < -Hsat)
5) Hysteretic with the major curve being an ellipse and described by the EFG formulation with Preisach's model
6) Hysteretic with the major and minor curves being ellipses and described by Preisach's model
7) Hysteretic curve with 4 parameters Preisach model (Br, Bsat, Hc, s <-> Wh)

The magnetic properties used in any model come from experimental measurements on a given material. The user can choose a specific model with "mattype" and modify the corresponding magnetic properties in PEM_IO.m. It is also possible to list the properties in order to create equivalent permeability curves at different temperatures, for example.

## Algorithm
There are two main steps in the algorithm for the Power-equivalent model:
### FEM_transitoire_fortran
The first step of the code solves the 1-D time-transient slab diffusion problem with the scalar Preisach model to compute the eddy current and hysteresis losses as a function of x. The problem is solved in Fortran 90 using the finite element method (FEM). The losses curve are then stored in the folder. The current values used in our software represent magnetic measurements done on AISI4340 steel at Polytechnique Montreal. Corresponding loss curves are stored in Results\Resultats_transitoire_fortran. A user could describe other materials with this model.

### FEM_PEM_1D
The second step uses the losses created from step 1 as inputs to solve the power-equivalent model. The equivalent permeability curves are then stored in Results\Resultats_mu.

## Citation
Please refer to CITATION.cff to cite this repository
