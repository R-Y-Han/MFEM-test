#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

real_t NORM2(real_t *x)
{
  return sqrt(x[0]*x[0]+x[1]*x[1]);
}
int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = "./icf.mesh";
   int ref_levels = 0;
   int order = 2;
   bool visualization = false;
   const char *outputname = "./opt.mesh";

   bool normalize = false;

   int precision = 8;
   cout.precision(precision);

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh uniformly.");
   args.AddOption(&order, "-o", "--order",
                  "Order (degree) of the finite elements.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&normalize, "-nor", "--normalize", "-no-nor",
                  "--no-normalization",
                  "Enable or disable normalization.");
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

   // 4. Refine the mesh to increase the resolution. In this example we do
   //    'ref_levels' of uniform refinement, where 'ref_levels' is a
   //    command-line parameter.
   for (int lev = 0; lev < ref_levels; lev++)
   {
      mesh->UniformRefinement();
   }
   const int NE = mesh->GetNE();

   real_t vert[9][2];
   const int vid[9] = {14, 110, 109, 108, 6, 102, 101, 100, 2};
   for (int i = 0; i < 9; i++)
   {
      mesh->GetNode(vid[i], vert[i]);
      for (int j = 0; j < dim; j++)
      {
         vert[i][j] *= 1.2;
      }
      mesh->AddVertex(vert[i]);
   }

   const int newv[9] = {NE, NE + 1, NE + 2, NE + 3, NE + 4, NE + 5, NE + 6, NE + 7};

   const int el[8][4] = {
      {14, newv[0], newv[1], 110}, {110, newv[1], newv[2], 109},
      {109, newv[2], newv[3], 108}, {108, newv[3], newv[4], 6},
      {6, newv[4], newv[5], 102}, {102, newv[5], newv[6], 101},
      {101, newv[6], newv[7], 100}, {100, newv[7], newv[8], 2}
   };

   for (int i = 0; i < 8; i++)
   {
      mesh->AddQuad(el[i],1);
   }
   mesh->FinalizeQuadMesh(1,1,true);


   H1_FECollection fec(order, dim);
   FiniteElementSpace fes(mesh, &fec, dim, 0);
   mesh->SetNodalFESpace(&fes);

   GridFunction x(&fes);
   mesh->SetNodalGridFunction(&x);


   mesh->Save(outputname);


   // 10. Free the used memory.
   delete mesh;

   return 0;
}