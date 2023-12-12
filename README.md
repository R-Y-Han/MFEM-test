# README

Small test cases for MFEM

## 1. GMRES

The Newton iteration method with a GMRES solver. This solver solves the equation (linear or nonlinear)
$$
f(x) = 0
$$
and use an approximate to calculate the Jacobian of $f$.

The nonlinear form $f$ need to be defined by user as a derived-class of operator. The Initial value should be chosen carefully or it may not give the wanted value.

TODO:
  - [ ] Preconditioner
  - [ ] The chosen of $\delta x_{k0}$

## 2. HeatCondution

The test case solves the PDE $\frac{dE}{dt} = \nabla \cdot (K \nabla T)$, where $E = cv T + \mathcal{P} T^4$.

The basis use DG basis, the rhs is discreted with IP method. The iteration uses the GMRES solver