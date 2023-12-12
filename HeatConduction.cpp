#include "mfem.hpp"
#include "GMRES.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

// class for solving dE/dt = div (K \nabla T), E = cv T + P_ T^4 (where \rho = 1)
class ImplicitConduction : public TimeDependentOperator
{
protected:
   FiniteElementSpace &fespace;
   Array<int> ess_tdof_list; // this list remains empty for pure Neumann b.c.

   BilinearForm *M; // mass matrix
   NonlinearForm *D; // \int div(K \nabla T) b
   SparseMatrix Mmat;

   // The nonlinear function for GMRES
   // Use backward Euler, M(E^n+1 - E^n) / dt = D(E^n+1)
   class EquationForm : public Operator
   {
   private:
      SparseMatrix *mass;
      NonlinearForm *diff;
      double dt;
      mutable Vector E;
      mutable Vector z;
   public:
      EquationForm(SparseMatrix &ma, NonlinearForm &di, const Vector &u_, double dt_)
         : Operator(u_.Size()), mass(&ma), diff(&di)
           { dt = dt_; E = u_; }
      // f(k) = M(k - E) - dt D(k) = 0
      virtual void Mult(const Vector &k, Vector &y) const override
      {
         z.SetSize(E.Size());
         diff->Mult(k, y);
         y *= -dt;
         subtract(k, E, z);
         mass->AddMult(z, y);
      }
   };
   EquationForm *F; // form the system F(E) = M(dE/dt) - D(E + dt dE/dt) = 0
   double current_dt;

   GMRESSolver T_solver; // Implicit solver for F(E) = 0
   DSmoother T_prec;  // Preconditioner for the implicit solver

   double kappa, P;

   mutable Vector z; // auxiliary vector

public:
   ImplicitConduction(FiniteElementSpace &f, const double kap, const double _P,
                      const Vector &u);

   virtual void Mult(const Vector &u, Vector &du_dt) const;
   /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
       This is the only requirement for high-order SDIRK implicit integration.*/
   virtual void ImplicitSolve(const double dt, const Vector &u, Vector &k);

   /// Update the diffusion BilinearForm K using the given true-dof vector `u`.
   void SetParameters(const Vector &u);

   virtual ~ImplicitConduction();
};

double ComputeT(const double &E);
void ComputeGradT(const double &E, const Vector &dE, Vector &dT);

class DomainIntegrator : public NonlinearFormIntegrator
{
private:
   const int dim;
   Coefficient *Q;
   void GetState(const Vector &U, const FiniteElement &el, DenseMatrix &state);

public:
   DomainIntegrator(const int dim_, Coefficient &coeff) : dim(dim_), Q(&coeff) { }
   virtual void AssembleElementVector(const FiniteElement &el,
                                      ElementTransformation &Tr,
                                      const Vector &elfun, Vector &elvect);
};

class FaceIntegrator : public NonlinearFormIntegrator
{
private:
   const double dim;
   Coefficient *Q;
   Vector shape1;
   Vector shape2;
   Vector funval1;
   Vector funval2;
   Vector nor;
   Vector flux;
   double fluxN;

public:
   FaceIntegrator(const int dim_, Coefficient &Q_)
      : dim(dim_), nor(dim), Q(&Q_),
        funval1(1), funval2(1), flux(dim) { }

   virtual void AssembleFaceVector(const FiniteElement &el1,
                                   const FiniteElement &el2,
                                   FaceElementTransformations &Tr,
                                   const Vector &elfun, Vector &elvect);
};

double InitialTemperature(const Vector &x);
double InitialEnergy(const Vector &x);

class NonlinearEq : public Operator
{
public:
   NonlinearEq(int dim) : Operator(dim) { }
   // y = f(x) = x^2 - 4
   virtual void Mult(const Vector &x, Vector &y) const override
   {
      y.SetSize(1);
      y(0) = sin(x(0)) - 1;
   }
};

const int dim = 2;
const double cv = 1.0;
double kappa = 0.1;
double IP = 10.;
double P_ = 0.;

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

   const char *outputname = "./HeatConduction.gf";

   int precision = 8;
   cout.precision(precision);

   OptionsParser args(argc, argv);
   args.AddOption(&IP, "-IP", "--IPparam",
                  "IP method panelty param.");
   args.AddOption(&P_, "-P", "--P",
                  "radiation effect for E_r = \\nu P T^4.");
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
   args.AddOption(&kappa, "-k", "--kappa",
                  "Kappa coefficient offset.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&visit, "-visit", "--visit-datafiles", "-no-visit",
                  "--no-visit-datafiles",
                  "Save data files for VisIt (visit.llnl.gov) visualization.");
   args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                  "Visualize every n-th timestep.");
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

NonlinearEq *f = new NonlinearEq(1);
NonlinearGMRES g_solver;
g_solver.SetRelTol(0.0);
g_solver.SetAbsTol(1e-7);
g_solver.SetMaxIter(30);
g_solver.SetPrintLevel(0);
g_solver.SetMaxKrylovIter(2);
g_solver.SetOperator(*f);
Vector x(1), y(1);
x = 1.;
g_solver.Mult(x,y);
y.Print();
exit(0);


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
   // ConductionOperator oper(*fespace, kappa, u);
   ImplicitConduction oper(fespace, kappa, P_, u);

   u_gf.SetFromTrueDofs(u);

   {
      ofstream omesh("./HeatConduction.mesh");
      omesh.precision(precision);
      mesh->Print(omesh);
      ofstream osol("./HeatConduction-init.gf");
      osol.precision(precision);
      u_gf.Save(osol);
   }

   VisItDataCollection visit_dc("HeatConduction", mesh);
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

ImplicitConduction::ImplicitConduction(FiniteElementSpace &f, const double kap,
                                       const double _P, const Vector &u)
   : TimeDependentOperator(f.GetTrueVSize(), 0.0), fespace(f), M(NULL), D(NULL),
     F(NULL), current_dt(0.0), z(height)
{
   const double rel_tol = 1e-8;

   M = new BilinearForm(&fespace);
   M->AddDomainIntegrator(new MassIntegrator());
   M->Assemble();
   M->FormSystemMatrix(ess_tdof_list, Mmat);

   kappa = kap;
   P = _P;

   T_solver.iterative_mode = false;
   T_solver.SetRelTol(rel_tol);
   T_solver.SetAbsTol(0.0);
   T_solver.SetMaxIter(100);
   T_solver.SetPrintLevel(0);
   // T_solver.SetPreconditioner(T_prec);

   SetParameters(u);
}

void ImplicitConduction::Mult(const Vector &u, Vector &du_dt) const
{
   MFEM_ABORT("Mult not implemented.");
}

void ImplicitConduction::ImplicitSolve(const double dt,
                                       const Vector &u, Vector &du_dt)
{
   // Solve F(du/dt) = 0, F(du/dt) = M du/dt - D(u + dt * du/dt)
   SetParameters(u);
   if (!F)
   {
      F = new EquationForm(Mmat, *D, u, dt);
      current_dt = dt;
      T_solver.SetOperator(*F);
   }
   MFEM_VERIFY(dt == current_dt, "");
   Vector b(fespace.GetTrueVSize());
   b = 1.;
   T_solver.Mult(b, du_dt);
}

void ImplicitConduction::SetParameters(const Vector &u)
{
   GridFunction K_gf(&fespace);
   K_gf.SetFromTrueDofs(u);
   for (int i = 0; i < K_gf.Size(); i++)
   {
      // K_gf(i) = kappa + alpha*K_gf(i);
      K_gf(i) = ComputeK(K_gf(i));
   }

   delete D;
   D = new NonlinearForm(&fespace);

   GridFunctionCoefficient u_coeff(&K_gf);

   D->AddDomainIntegrator(new DomainIntegrator(dim, u_coeff));
   D->AddInteriorFaceIntegrator(new FaceIntegrator(dim, u_coeff));
   delete F;
   F = NULL; // re-compute T on the next ImplicitSolve
}


ImplicitConduction::~ImplicitConduction()
{
   delete M;
   delete D;
   delete F;
}

void DomainIntegrator::AssembleElementVector(const FiniteElement &el,
                                             ElementTransformation &Tr,
                                             const Vector &elfun, Vector &elvect)
{
  const int elem = Tr.ElementNo;
  const int dof = el.GetDof();

  elvect.SetSize(dof);
  elvect = 0.0;
  DenseMatrix elvect_mat(elvect.GetData(), dof, 1);

  Vector shape(dof);
  DenseMatrix dshape(dof,dim), gshape(dof,dim);
  DenseTensor flux;
  DenseMatrix Jadj(dim); // J^*
  double w;

  int intorder = 2 * el.GetOrder();
  const IntegrationRule *ir = &IntRules.Get(Tr.GetGeometryType(), intorder);

  // Note elfun here is just E
  for (int i = 0; i < ir->GetNPoints(); i++) {
    const IntegrationPoint &ip = ir->IntPoint(i);

    el.CalcShape(ip, shape);
    double E = elfun * shape;
    el.CalcDShape(ip, dshape); // derivative on reference cell
    Vector dE(dim);
    dshape.MultTranspose(elfun, dE);

    Tr.SetIntPoint(&ip);
    CalcAdjugate(Tr.Jacobian(), Jadj);
    w = Tr.Weight();
    w = ip.weight / w;
   //  if (Q) {
   //    w *= Q->Eval(Tr, ip);
   //  }
   double T = ComputeT(E);
   w *= ComputeK(T);

    gshape.SetSize(dshape.Height(), dshape.Width());
    Mult(dshape, Jadj, gshape); // \nabla_x b |J|

    Vector dT;
    ComputeGradT(E, dE, dT);
    // \nabla_x T = \nabla_\xi T Jadj / |J|
    Jadj.MultTranspose(dT, dT);

    dT *= w;
    gshape.AddMult(dT, elvect);
  } /* for (int i = 0; i < IntRule->GetNPoints(); i++) */
}

void FaceIntegrator::AssembleFaceVector(const FiniteElement &el1,
                                        const FiniteElement &el2,
                                        FaceElementTransformations &Tr,
                                        const Vector &elfun, Vector &elvect)
{
   // Compute the term <F.n(u),[w]> on the interior faces.
   const int dof1 = el1.GetDof();
   const int dof2 = el2.GetDof();

   shape1.SetSize(dof1);
   shape2.SetSize(dof2);

   elvect.SetSize(dof1 + dof2);
   elvect = 0.0;

   DenseMatrix elfun1_mat(elfun.GetData(), dof1, 1);
   DenseMatrix elfun2_mat(elfun.GetData() + dof1, dof2, 1);
   Vector Evec1(elfun.GetData(), dof1);
   Vector Evec2(elfun.GetData(),dof2);

   Vector elvect1_mat(elvect.GetData(), dof1);
   Vector elvect2_mat(elvect.GetData() + dof1, dof2);

   DenseMatrix dshape1(dof1,dim), dshape2(dof1,dim);
   DenseMatrix Jadj1(dim), Jadj2(dim);

   // Integration order calculation from DGTraceIntegrator
   int intorder; // 这个intorder是积分的阶数
   if (Tr.Elem2No >= 0)
      intorder = (min(Tr.Elem1->OrderW(), Tr.Elem2->OrderW()) +
                  2*max(el1.GetOrder(), el2.GetOrder()));
   else
   {
      intorder = Tr.Elem1->OrderW() + 2*el1.GetOrder();
   }
   if (el1.Space() == FunctionSpace::Pk)
   {
      intorder++;
   }
   const IntegrationRule *ir = &IntRules.Get(Tr.GetGeometryType(), intorder);

   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);

      Tr.SetAllIntPoints(&ip); // set face and element int. points

      // Calculate basis functions on both elements at the face
      el1.CalcShape(Tr.GetElement1IntPoint(), shape1);
      el2.CalcShape(Tr.GetElement2IntPoint(), shape2);

      elfun1_mat.MultTranspose(shape1, funval1);
      elfun2_mat.MultTranspose(shape2, funval2);

      // Get the normal vector and the flux on the face
      CalcOrtho(Tr.Jacobian(), nor);
      
      // Get flux = {\nabla T} + mu [T]
      // here funval = E
      double T1 = ComputeT(funval1(0));
      double T2 = ComputeT(funval2(0));
      double mu = IP * Tr.Weight() / Tr.Elem1->Weight();
      el1.CalcDShape(Tr.GetElement1IntPoint(), dshape1);
      el2.CalcDShape(Tr.GetElement2IntPoint(), dshape2);
      CalcAdjugate(Tr.Elem1->Jacobian(), Jadj1);
      CalcAdjugate(Tr.Elem2->Jacobian(), Jadj2);
      Vector dE1(dim), dE2(dim), dT1(dim), dT2(dim);
      Jadj1.MultTranspose(Evec1, dE1);
      Jadj2.MultTranspose(Evec2, dE2);
      ComputeGradT(funval1(0), dE1, dT1);
      ComputeGradT(funval2(0), dE2, dT2);
      // \nabla_x T = \nabla_xi T Jadj / |J|
      Jadj1.MultTranspose(dT1, dT1);
      dT1 /= Tr.Elem1->Weight();
      Jadj2.MultTranspose(dT2, dT2);
      dT2 /= Tr.Elem2->Weight();
      add(0.5, dT1, dT2, flux);
      add(flux, mu*(T1-T2), nor, flux); // [T] = T1 nor1 + T2 nor2

      fluxN = flux * nor;


      double w = ip.weight;
      // if (Q) { w *= Q->Eval(Tr, ip); }
      double K = ComputeK(T1);
      w *= K;

      fluxN *= w;

      for (int s = 0; s < dof1; s++)
      {
         elvect1_mat(s) += fluxN * shape1(s); // 注意N的方向，elvect=\int f_k*n b_s
      }
      for (int s = 0; s < dof2; s++)
      {
         elvect2_mat(s) -= fluxN * shape2(s);
      }
   }
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

double ComputeT(const double &E)
{
   if (P_ < 1e-6) {
   double b1, b2, b3;
      b1 = (E - 0.0) / cv; // note velocity = 0 here
      b1 = std::max(b1, 1e-9);
      b2 = -pow(b1, 4) / (1.0 * cv); // rho = 1
      b3 = -(2. * pow(b1,3) * b2) / (1. * cv);
      return b1 + b2 * P_ + b3 * P_ * P_;
   } else {
      double c1, c2, temp1, temp2, temp3, s;
      c1 = (1. * cv) / P_;
      c2 = -1. * (E - 0.) / P_;

      temp3 = pow(c1,4) / 256.0 - pow(c2,3) / 27.0;
      temp1 = c1 * c1 / 16.0 + sqrt(temp3);
      temp2 = -c1 * c1 / 16.0 + sqrt(temp3);
      s = pow(temp1,(1./3.)) - pow(temp2,(1./3.));

      double T = - 2. * s + 2. * c1 / sqrt(2. * s);
      T = 0.5 * (-sqrt(2. * s) + sqrt(T));

      if (T < 1e-6) {
         T = (E - 0.) / cv;
      }
      return T;
   }
   return 0;
}

void ComputeGradT(const double &E, const Vector &dE, Vector &dT)
{
   double T = ComputeT(E);
   dT = dE;
   dT /= (cv + 4 * P_ *pow(T,3));
}

double ComputeK(const double &T)
{
   double K;
   K = kappa;
   return K;
}