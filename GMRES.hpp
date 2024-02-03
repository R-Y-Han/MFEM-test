#ifndef _GMRES_
#define _GMRES_

#include "mfem.hpp"

using std::cout;
namespace mfem{

// Class for solving the nonlinear equation f(x) = 0 using Newton method with
// GMRES solver, in this method the Jacobian of f(x) is to be approximated
class NonlinearGMRES : public IterativeSolver
{
protected:
  int max_Krylov_iter;
  double linear_tor; // for termination of Krylov iteration

  mutable bool Krylov_converged;
  mutable bool Krylov_basis_updated;
  mutable int basis_level; // point the number of basis generated, start from 0
  mutable Vector fk, dx_k, dx0, r0;
  mutable Vector beta, Jm_r;
  mutable Vector z;
  mutable DenseMatrix Jr;

  mutable CGSolver B_solver; // solver for inverting matrix in LS when computing beta
  mutable DSmoother B_prec;

  void GetBeta(const Vector &f, const Vector &x, const int l) const;
  void GetKrylovBasis(const Vector &x, const Vector &r, const int l) const;
  void JacobiMult(const Vector &x, const Vector &y, Vector &Jy) const;

public:
  NonlinearGMRES()
  {
    max_Krylov_iter = 100;
    linear_tor = 1e-2;

    B_solver.iterative_mode = false;
    B_solver.SetRelTol(1e-8);
    B_solver.SetAbsTol(0.0);
    B_solver.SetMaxIter(30);
    B_solver.SetPrintLevel(-1);
    B_solver.SetPreconditioner(B_prec);
  }

  void SetMaxKrylovIter(int step) { max_Krylov_iter = step; }
  void SetLinearTor(int tor) { linear_tor = tor; }

  // Solve for f(y) = 0, x is the initial input
  virtual void Mult(const Vector &x, Vector &y) const;
};

void NonlinearGMRES::Mult(const Vector &x, Vector &y) const
{
  MFEM_ASSERT(oper != NULL, "the Operator is not set (use SetOperator).");
  MFEM_ASSERT(prec != NULL, "the Solver is not set (use SetSolver).");

  int k;
  double norm0, norm, norm_goal, norm_Krylov, norm_Krylov_goal;
  const bool have_x = (x.Size() == Height());
  dx_k.SetSize(Height());
  dx0.SetSize(Height());
  r0.SetSize(Height());
  y.SetSize(Height());
  Jr.SetSize(max_Krylov_iter+1, Height());

  // Initial state
  if (have_x) { y = x; }
  else { y = 0.0; }

  oper->Mult(y, fk);
  norm0 = norm = initial_norm = Norm(fk);
  if (print_options.first_and_last && !print_options.iterations)
  {
    mfem::out << "Newton iteration " << std::setw(2) << 0
              << " : ||f(x)|| = " << norm << "...\n";
  }

  norm_goal = std::max(rel_tol * norm, abs_tol);

  // prec->iterative_mode = false;

  // Newton iteration, x_{k+1} = x_k + dx_k
  for (k = 0; true; k++) {
    MFEM_ASSERT(IsFinite(norm), "norm = " << norm);
    if (print_options.iterations) {
      mfem::out << "Newton iteration " << std::setw(2) << k
                << " : ||f(x)|| = " << norm;
      if (k > 0) {
        mfem::out << ", ||f(x)||/||f(x)_0|| = " << norm/norm0;
      }
      mfem::out << '\n';
    }
    Monitor(k, norm, fk, y);

    if (norm <= norm_goal) {
      converged = true;
      break;
    }
    if (k >= max_iter) {
      converged = false;
      break;
    }

    // Krylov iteration for getting dx_k
    norm_Krylov_goal = std::max(linear_tor * norm, abs_tol);
    dx0 = 1e-7;
    JacobiMult(dx0, y, r0);
    r0 += fk;
    r0.Neg(); // r0 = -f(x_k) - J_k(dx_k0)
    Krylov_basis_updated = false;

    for (int l = 0; true; l++) {
      if (Norm(r0) <= norm_Krylov_goal) {
        dx_k = dx0;
        Krylov_converged = true;
        break;
      }
      GetKrylovBasis(y, r0, l);
      GetBeta(fk, y, l);

      // Get dx_kl = dx_0 + sum_{m=0}^{l-1} beta_m J_k^m(r0)
      dx_k = dx0;
      for (int m=0; m < l; m++) {
        Jr.GetRow(m, Jm_r);
        add(dx_k, beta(m), Jm_r, dx_k);
      }

      // Calculate residual ||f(x_k) + J_k (dx_kl)||_2
      JacobiMult(y, dx_k, z);
      z += fk;

      norm_Krylov = Norm(z);
      MFEM_ASSERT(IsFinite(norm_Krylov), "Krylov norm = " << norm_Krylov);
      if (print_options.iterations) {
        mfem::out << "\t|-" << l << ". Krylov iteration " << std::setw(2) << l
                  << " : ||f(x_k) + J_k(dx_kl)|| / ||f(x_k)|| = "
                  << norm_Krylov / norm << '\n';
      }
      Monitor(l, norm_Krylov, z, y);

      if (norm_Krylov <= norm_Krylov_goal) {
        Krylov_converged = true;
        break;
      }
      if (l >= max_Krylov_iter) {
        Krylov_converged = false;
        break;
      }
      // break if Krylov basis has 0 term
      Vector bl;
      Jr.GetRow(basis_level, bl);
      if (Norm(bl) <= 1e-6) {
        break;
      }
    }

    add(y, dx_k, y); // y = x_{k+1}

    oper->Mult(y, fk);
    norm = Norm(fk);
  }
}

void NonlinearGMRES::JacobiMult(const Vector &x, const Vector &y, Vector &Jy) const
{
  const int size = x.Size();
  const double eta = 1e-5;

  Jy.SetSize(size);

  double norm_y = sqrt(y * y);
  if (norm_y == 0.0) {
    Jy = 0.;
    return ;
  }
  double sum = std::max(x.Sum(), 1e-6);
  double eps = eta * sum / (size * norm_y);
  eps = std::max(1e-6, eps);

  // Jy = [f(x+ eps y) - f(x)] / eps
  Vector a(size), f(size);
  oper->Mult(x, f);
  add(x, eps, y, a);
  oper->Mult(a, Jy); // f(x + eps y)
  subtract(Jy, f, Jy);
  Jy /= eps;
}

void NonlinearGMRES::GetKrylovBasis(const Vector &x, const Vector &r, const int l) const
{
  // Jr(m) = J^m (r0)
  if (!Krylov_basis_updated) {
    Jr = 0.0;
    Jr.SetRow(0, r);
    basis_level = 0;
  }
  Vector a, b;
  for (int m = basis_level+1; m <= l+1; m++) {
    Jr.GetRow(m-1, a); // a = J^{m-1} (r0)
    JacobiMult(x, a, b);
    Jr.SetRow(m, b);
    basis_level = m;
  }

  Krylov_basis_updated = true;
}

void NonlinearGMRES::GetBeta(const Vector &f, const Vector &x, const  int l) const
{
  // min ||f(x) + J(dx_k)||_2^2
  if (!l) { return ; }
  beta.SetSize(l);
  SparseMatrix M(l,l);
  Vector B(l);
  Vector b1, b2, c;

  JacobiMult(x, dx0, c);
  c += fk; // dx0, fk have been computed

  for (int i = 0; i < l; i++) {
    Jr.GetRow(i+1, b1);
    for (int j = 0; j < l; j++) {
      Jr.GetRow(j+1, b2);
      M.Set(i,j, b1 * b2);
      // M(i,j) = b1 * b2;
    }

    B(i) = -1. * (b1 * c);
  }

  M.Finalize();

  if (l <= 10) {
    DenseMatrix Mden, Minv;
    M.ToDenseMatrix(Mden);
    DenseMatrixInverse inv(&Mden);
    inv.Factor();
    inv.GetInverseMatrix(Minv);
    Minv.Mult(B, beta);
  } else {
    B_solver.SetOperator(M);
    B_solver.Mult(B, beta);
  }
}

} // mfem

#endif