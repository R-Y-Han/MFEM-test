# README

Small test cases for MFEM

## 1. HeatCondution

Test case mainly for nonlinear GMRES solver. The iteration class NonlinearGMRES use a Newton method with GMRES solver.

The test case solves the PDE $\frac{dE}{dt} = \nabla \cdot (K \nabla T)$, where $E = cv T + \mathcal{P} T^4$.

The basis use DG basis, the rhs is discreted with IP method.