#include "mfem.hpp"
#include "GMRES.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;


class NonlinearEq : public Operator
{
private:
  int d;
public:
  NonlinearEq(int dim) : Operator(dim), d(dim) { }
  virtual void Mult(const Vector &x, Vector &y) const override
  {
     y.SetSize(d);
     if (d == 1) {
       y(0) = exp(x(0)) - 2;
     }
     else if (d == 2) {
       y(0) = x(0) * x(0) + x(1) * x(1) - 1.;
       y(1) = x(0) * x(0) - x(1);
     }
     else {
       y = 0.0;
     }
  }
};

class ConductionOperator : public TimeDependentOperator
{
protected:
  FiniteElementSpace &fespace;
  Array<int> ess_tdof_list; // this list remains empty for pure Neumann b.c.

  BilinearForm *M;
  BilinearForm *D;

  SparseMatrix Mmat, Dmat;
  class DiscreteForm : public Operator
  {
    private:
      SparseMatrix *Mmat, *Dmat;
      Vector u;
      double dt;

      mutable Vector z; // auxiliary vector
    public:
      DiscreteForm(SparseMatrix &M, SparseMatrix &D,
                   Vector &u_, const double dt_)
        : Operator(M.Height())
      {
        Mmat = &M;
        Dmat = &D;
        u = u_;
        dt = dt_;
      }
      virtual void Mult(const Vector &k, Vector &y) const override
      {
        y.SetSize(k.Size());
        z.SetSize(k.Size());
        // y = (M + dt D) k - D u
        SparseMatrix *F;
        F = Add(1.0, *Mmat, dt, *Dmat);
        Mmat->Mult(k, y);
        F->Mult(k, y);
        Dmat->Mult(u, z);
        add(y, z, y);
      }
  };
  DiscreteForm *T;
  double current_dt;

  NonlinearGMRES T_solver; // Implicit solver for T = M + dt D
  //  DSmoother T_prec;  // Preconditioner for the implicit solver

   double kappa;
  //  FunctionCoefficient *K;

public:
   ConductionOperator(FiniteElementSpace &f, double kappa,
                      const Vector &u);

   /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
       This is the only requirement for high-order SDIRK implicit integration.*/
   virtual void ImplicitSolve(const double dt, const Vector &u, Vector &k);

   /// Update the diffusion BilinearForm K using the given true-dof vector `u`.
   void SetParameters(const Vector &u);

   virtual ~ConductionOperator();
};

const int dim = 2;
const double cv = 1.0;
double kappa = 0.1;
double IP = 10.;
double P_ = 0.;

double InitialTemperature(const Vector &x);
double ComputeK(const double &T);

int main(int argc, char *argv[])
{
  // 1. Parse command-line options.
  const char *mesh_file = "./2D.mesh";
  int ref_levels = 0;
  int order = 2;
  int ode_solver_type = 1;
  double t_final = 0.1;
  double dt = 1.0e-3;
  bool visualization = true;
  bool visit = false;
  int vis_steps = 5;

  const char *outputname = "./GMRES.gf";

  int precision = 8;
  cout.precision(precision);

  OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh",
                 "Mesh file to use.");
  args.AddOption(&ref_levels, "-r", "--refine",
                 "Number of times to refine the mesh uniformly.");
  args.AddOption(&order, "-o", "--order",
                 "Order (degree) of the finite elements.");
  args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                 "ODE solver: 1 - Backward Euler, 2 - SDIRK2, 3 - SDIRK3,\n\t"
                 "\t   11 - Forward Euler, 12 - RK2, 13 - RK3 SSP, 14 - RK4.");
  args.AddOption(&t_final, "-tf", "--t-final",
                 "Final time; start time is 0.");
  args.AddOption(&dt, "-dt", "--time-step",
                 "Time step.");
  args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                 "--no-visualization",
                 "Enable or disable GLVis visualization.");
  args.AddOption(&visit, "-visit", "--visit-datafiles", "-no-visit",
                 "--no-visit-datafiles",
                 "Save data files for VisIt (visit.llnl.gov) visualization.");
  args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                 "Visualize every n-th timestep.");
  args.AddOption(&IP, "-IP", "--IP",
                 "Interior Penalty parameter.");
  args.Parse();
  if (!args.Good())
  {
     args.PrintUsage(cout);
     return 1;
  }
  args.PrintOptions(cout);

  // 2. Read the mesh from the given mesh file. We can handle triangular,
  //    quadrilateral, tetrahedral and hexahedral meshes with the same code.
  Mesh *mesh = new Mesh(mesh_file, 1, 1);
  int dim = mesh->Dimension();

  // 3. Define the ODE solver used for time integration. Several implicit
  //    singly diagonal implicit Runge-Kutta (SDIRK) methods, as well as
  //    explicit Runge-Kutta methods are available.
  ODESolver *ode_solver;
  switch (ode_solver_type)
  {
     // Implicit L-stable methods
     case 1:  ode_solver = new BackwardEulerSolver; break;
     case 2:  ode_solver = new SDIRK23Solver(2); break;
     case 3:  ode_solver = new SDIRK33Solver; break;
     // Explicit methods
     case 11: ode_solver = new ForwardEulerSolver; break;
     case 12: ode_solver = new RK2Solver(0.5); break; // midpoint method
     case 13: ode_solver = new RK3SSPSolver; break;
     case 14: ode_solver = new RK4Solver; break;
     case 15: ode_solver = new GeneralizedAlphaSolver(0.5); break;
     // Implicit A-stable methods (not L-stable)
     case 22: ode_solver = new ImplicitMidpointSolver; break;
     case 23: ode_solver = new SDIRK23Solver; break;
     case 24: ode_solver = new SDIRK34Solver; break;
     default:
        cout << "Unknown ODE solver type: " << ode_solver_type << '\n';
        delete mesh;
        return 3;
  }

  // 4. Refine the mesh to increase the resolution. In this example we do
  //    'ref_levels' of uniform refinement, where 'ref_levels' is a
  //    command-line parameter.
  for (int lev = 0; lev < ref_levels; lev++)
  {
     mesh->UniformRefinement();
  }

  cout << "\n1. 1D nonlinear equation case.\n";
  NonlinearGMRES g_solver;
  g_solver.iterative_mode = false;
  g_solver.SetRelTol(0.0);
  g_solver.SetAbsTol(1e-7);
  g_solver.SetMaxIter(100);
  g_solver.SetPrintLevel(0);
  g_solver.SetMaxKrylovIter(100);
  NonlinearEq f(1);
  g_solver.SetOperator(f);
  Vector x(1), y(1);
  x = 0.;
  g_solver.Mult(x,y);
  cout << " x = ";
  y.Print();
  
  cout << "\n2. Nonlinear equations case\n";
  NonlinearEq A(2);
  g_solver.SetOperator(A);
  x.SetSize(2);
  x = 1.;
  g_solver.Mult(x, y);
  cout << "x = ";
  y.Print();
  cout << "\n";



  // 5. Define the vector finite element space representing the current and the
  //    initial temperature, u_ref.
  // H1_FECollection fe_coll(order, dim);
  DG_FECollection dg_coll(order, dim);
  FiniteElementSpace fespace(mesh, &dg_coll);

  int fe_size = fespace.GetTrueVSize();
  cout << "Number of temperature unknowns: " << fe_size << endl;

  GridFunction u_gf(&fespace);

  // 6. Set the initial conditions for u. All boundaries are considered
  //    natural.
  FunctionCoefficient u_0(InitialTemperature);
  u_gf.ProjectCoefficient(u_0);
  Vector u;
  u_gf.GetTrueDofs(u);

  // 7. Initialize the conduction operator and the visualization.
   ConductionOperator oper(fespace, kappa, u);

  u_gf.SetFromTrueDofs(u);

   {
      ofstream omesh("./GMRES.mesh");
      omesh.precision(precision);
      mesh->Print(omesh);
      ofstream osol("./GMRES-init.gf");
      osol.precision(precision);
      u_gf.Save(osol);
   }

   VisItDataCollection visit_dc("GMRES", mesh);
   visit_dc.RegisterField("Energy", &u_gf);
   if (visit)
   {
      visit_dc.SetCycle(0);
      visit_dc.SetTime(0.0);
      visit_dc.Save();
   }

   socketstream sout;
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      sout.open(vishost, visport);
      if (!sout)
      {
         cout << "Unable to connect to GLVis server at "
              << vishost << ':' << visport << endl;
         visualization = false;
         cout << "GLVis visualization disabled.\n";
      }
      else
      {
         sout.precision(precision);
         sout << "solution\n" << *mesh << u_gf;
         sout << "pause\n";
         sout << flush;
         cout << "GLVis visualization paused."
              << " Press space (in the GLVis window) to resume it.\n";
      }
   }

   // 8. Perform time-integration (looping over the time iterations, ti, with a
   //    time-step dt).
   ode_solver->Init(oper);
   double t = 0.0;

   bool last_step = false;
   for (int ti = 1; !last_step; ti++)
   {
      if (t + dt >= t_final - dt/2)
      {
         last_step = true;
      }

      ode_solver->Step(u, t, dt);

      if (last_step || (ti % vis_steps) == 0)
      {
         cout << "step " << ti << ", t = " << t << endl;

         u_gf.SetFromTrueDofs(u);

         if (visualization)
         {
            sout << "solution\n" << *mesh << u_gf << flush;
         }

         if (visit)
         {
            visit_dc.SetCycle(ti);
            visit_dc.SetTime(t);
            visit_dc.Save();
         }
      }
      oper.SetParameters(u);
   }

   // 9. Save the final solution. This output can be viewed later using GLVis:
   //    "glvis -m ex16.mesh -g ex16-final.gf".
   {
      ofstream osol(outputname);
      osol.precision(precision);
      u_gf.Save(osol);
   }

   // 10. Free the used memory.
   delete ode_solver;
   delete mesh;

   return 0;
}

ConductionOperator::ConductionOperator(FiniteElementSpace &f,
                                       double kap, const Vector &u)
   : TimeDependentOperator(f.GetTrueVSize(), 0.0), fespace(f), M(NULL), D(NULL),
     T(NULL), current_dt(0.0)
{
  const double rel_tol = 1e-8;

  GridFunction coeff_gf(&fespace);
  for (int i = 0; i < fespace.GetNDofs(); i++) {
   coeff_gf(i) = cv;
  }
  GridFunctionCoefficient Coeff_T(&coeff_gf);

  M = new BilinearForm(&fespace);
  M->AddDomainIntegrator(new MassIntegrator(Coeff_T));
  M->Assemble();
  M->FormSystemMatrix(ess_tdof_list, Mmat);

  kappa = kap;

  T_solver.iterative_mode = false;
  T_solver.SetRelTol(rel_tol);
  T_solver.SetAbsTol(0.0);
  T_solver.SetMaxIter(50);
  T_solver.SetPrintLevel(0);
  T_solver.SetMaxKrylovIter(100);

  SetParameters(u);
}

void ConductionOperator::ImplicitSolve(const double dt,
                                       const Vector &u, Vector &du_dt)
{
   // Solve the equation:
   //    du_dt = M^{-1}*[D(u + dt*du_dt)]
   // for du_dt, where K is linearized by using u from the previous timestep
   if (!T)
   {
      Vector u_(u);
      T = new DiscreteForm(Mmat, Dmat, u_, dt); // Dmat is infact -D_T
      current_dt = dt;
      T_solver.SetOperator(*T);
   }
   MFEM_VERIFY(dt == current_dt, ""); // SDIRK methods use the same dt
   T_solver.Mult(u, du_dt);
  //  for (int i = 0; i < fespace.GetFE(0)->GetDof(); i++) {
  //   cout << du_dt(i) * dt << "\t";
  //  }
  //  cout<<"\n";
}

void ConductionOperator::SetParameters(const Vector &u)
{
   GridFunction K_gf(&fespace);
   K_gf.SetFromTrueDofs(u);
   for (int i = 0; i < K_gf.Size(); i++)
   {
      // K_gf(i) = kappa + alpha*K_gf(i);
      K_gf(i) = ComputeK(K_gf(i));
   }

   delete D;
   D = new BilinearForm(&fespace);

   GridFunctionCoefficient u_coeff(&K_gf);

   D->AddDomainIntegrator(new DiffusionIntegrator(u_coeff));
   D->AddInteriorFaceIntegrator(new DGDiffusionIntegrator(u_coeff, -1, IP));
   D->Assemble();
   D->Finalize();
   D->FormSystemMatrix(ess_tdof_list, Dmat);
   delete T;
   T = NULL; // re-compute T on the next ImplicitSolve
}

ConductionOperator::~ConductionOperator()
{
   delete T;
   delete M;
   delete D;
}


double InitialTemperature(const Vector &x)
{
   if (fabs(x(0)) < 0.1 && fabs(x(1)) < 0.1)
   {
      return 10.0;
   }
   else
   {
      return 0.1;
   }
}

double InitialEnergy(const Vector &x)
{
   double T = InitialTemperature(x);
   return cv * T + P_ * pow(T, 4);
}

double ComputeK(const double &T)
{
   return kappa;
}