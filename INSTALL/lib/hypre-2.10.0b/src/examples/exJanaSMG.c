/*
   Example Jana

   Interface:      Structured interface (Struct)

   Compile with:   make exJana

   Sample run:     mpirun -np 1 exJana -n 64 -solver 10 -K 3 -B 0 -C 1 -U0 2 -F 4

   To see options: exJana -help

   Description:    This example differs from the previous structured example
                   (Example 3) in that a more sophisticated stencil and
                   boundary conditions are implemented. The method illustrated
                   here to implement the boundary conditions is much more general
                   than that in the previous example.  Also symmetric storage is
                   utilized when applicable.

                   This code solves the Helmhotz equation:
                    -----------------
                    -grad(K * grad (u) )  + u*K = F*K;
                    u(0) = 0
 
                    K - phase field function
                    eps - smoothenitnhg width
                    boundary condition u = U0.
 
                    optionK = 0 - test on squaer
                    optionK = 1 - test on sphere
                    optionK = 2 - test on sphere (inverse)
 
 
 
                    The domain is split into N x N
                   processor grid.  Thus, the given number of processors should
                   be a perfect square. Each processor has a n x n grid, with
                   nodes connected by a 5-point stencil. Note that the struct
                   interface assumes a cell-centered grid, and, therefore, the
                   nodes are not shared.

                   To incorporate the boundary conditions, we do the following:
                   Let x_i and x_b be the interior and boundary parts of the
                   solution vector x. If we split the matrix A as
                             A = [A_ii A_ib; A_bi A_bb],
                   then we solve
                             [A_ii 0; 0 I] [x_i ; x_b] = [b_i - A_ib u_0; u_0].
                   Note that this differs from the previous example in that we
                   are actually solving for the boundary conditions (so they
                   may not be exact as in ex3, where we only solved for the
                   interior).  This approach is useful for more general types
                   of b.c.
*/

#include <math.h>
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE_struct_ls.h"

#ifdef M_PI
  #define PI M_PI
#else
  #define PI 3.14159265358979
#endif

#include "vis.c"

/* Macro to evaluate a function F in the grid point (i,j) */
#define Eval(F,i,j) (F( (ilower[0]+(i))*h, (ilower[1]+(j))*h ))
#define bcEval(F,i,j) (F( (bc_ilower[0]+(i))*h, (bc_ilower[1]+(j))*h ))

int optionK, optionB, optionPSI, optionU0, optionF;

/* Phase field coefficient */
double K(double x, double y)
{
    // allocate variables
    double radius = 0.25;
    double eps = 0.007;
    double tau = 1.e-10;
    double r, psi, rx, ry, dist;
    
    
    // get sign distance function r at (x,y) for given test case=optionK

    switch (optionK)
    {
        case 0:  // square
            // shif to [0,0] centred square
            rx = x - 0.5;
            ry = y - 0.5;
            
            if (fabs(rx) > fabs(ry))
                r = fabs(rx) - radius;
            else
                r = fabs(ry) - radius;
            
            break;
            
        case 1:  // sphere
            dist = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5));
            r = dist - radius;   // sign dist function, <0 in domain, >0 out
            break;
            
        case 2:  // sphere
            dist = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5));
            r = dist - radius;   // sign dist function, <0 in domain, >0 out
            break;

        default:
            break;
    }
    
    // compute phase field function at (x,y)
    psi = 0.5 * (1.0 - tanh(3.*r/eps) ) + tau;
    
    return psi;
}


/* Boundary condition */
double U0(double x, double y)
{
    return 0.0;
}

/* Right-hand side */
double F(double x, double y)
{
    double tau = 1.e-10;
    double psi = K(x,y) - tau;
    
    // allocate variables
    double dist, c, alpha, rhs, ifrac;
    
    switch (optionK) {
        case 0:  // square
            c = 4.;
            rhs = - (2. * c*c * PI*PI + 1.) * cos(c*PI*x) * cos(c*PI*y);
            break;
            
        case 1:  // sphere
            dist = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5));
            c = 4.;
            
            ifrac = 1./(c*c*PI*PI);
            alpha = c*PI*dist;
            rhs = sin(alpha) * (2.0 - ifrac ) + alpha*cos(alpha)*(1.0 + ifrac);
            break;
            
        case 2:  // sphere
            dist = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5));
            c = 4.;
            
            ifrac = 1./(c*c*PI*PI);
            alpha = c*PI*dist;
            rhs = sin(alpha) * (-2.0 + ifrac ) - alpha*cos(alpha)*(1.0 + ifrac) - 1.;
            break;

        default:
            break;
    }
    
    return -1.0 * psi * rhs ;
}

int main (int argc, char *argv[])
{
   int i, j, k;

   int myid, num_procs;

   int n, N, pi, pj;
   double h, h2;
   int ilower[2], iupper[2];

   int solver_id;
   int n_pre, n_post;
   int rap, relax, skip, sym;
   int time_index;

   int num_iterations;
   double final_res_norm;

   int vis;

   HYPRE_StructGrid     grid;
   HYPRE_StructStencil  stencil;
   HYPRE_StructMatrix   A;
   HYPRE_StructVector   b;
   HYPRE_StructVector   x;
   HYPRE_StructSolver   solver;

   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

   /* Set default parameters */
   n         = 256;
   optionK   = 0;
    optionB   = 0;
    optionPSI = 0;
   optionU0  = 0;
   optionF   = 0;
   solver_id = 0;
   n_pre     = 1;
   n_post    = 1;
   rap       = 0;
   relax     = 1;
   skip      = 0;
   sym       = 0;

   vis       = 0;

   /* Parse command line */
   {
      int arg_index = 0;
      int print_usage = 0;

      while (arg_index < argc)
      {
         if ( strcmp(argv[arg_index], "-n") == 0 )
         {
            arg_index++;
            n = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-K") == 0 )
         {
            arg_index++;
            optionK = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-PSI") == 0 )
         {
             arg_index++;
             optionPSI = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-U0") == 0 )
         {
            arg_index++;
            optionU0 = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-F") == 0 )
         {
            arg_index++;
            optionF = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-solver") == 0 )
         {
            arg_index++;
            solver_id = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-v") == 0 )
         {
            arg_index++;
            n_pre = atoi(argv[arg_index++]);
            n_post = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-rap") == 0 )
         {
            arg_index++;
            rap = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-relax") == 0 )
         {
            arg_index++;
            relax = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-skip") == 0 )
         {
            arg_index++;
            skip = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-sym") == 0 )
         {
            arg_index++;
            sym = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-vis") == 0 )
         {
            arg_index++;
            vis = 1;
         }
         else if ( strcmp(argv[arg_index], "-help") == 0 )
         {
            print_usage = 1;
            break;
         }
         else
         {
            arg_index++;
         }
      }

      if ((print_usage) && (myid == 0))
      {
         printf("\n");
         printf("Usage: %s [<options>]\n", argv[0]);
         printf("\n");
         printf("  -n  <n>             : problem size per processor (default: 8)\n");
         printf("  -K  <K>             : choice for the diffusion coefficient (default: 1)\n");
         printf("  -U0 <U0>            : choice for the boundary condition (default: 0)\n");
         printf("  -F  <F>             : choice for the right-hand side (default: 1) \n");
         printf("  -solver <ID>        : solver ID\n");
         printf("                        0  - SMG \n");
         printf("                        1  - PFMG\n");
         printf("                        10 - CG with SMG precond (default)\n");
         printf("                        11 - CG with PFMG precond\n");
         printf("                        17 - CG with 2-step Jacobi\n");
         printf("                        18 - CG with diagonal scaling\n");
         printf("                        19 - CG\n");
         printf("                        30 - GMRES with SMG precond\n");
         printf("                        31 - GMRES with PFMG precond\n");
         printf("                        37 - GMRES with 2-step Jacobi\n");
         printf("                        38 - GMRES with diagonal scaling\n");
         printf("                        39 - GMRES\n");
         printf("  -v <n_pre> <n_post> : number of pre and post relaxations\n");
         printf("  -rap <r>            : coarse grid operator type\n");
         printf("                        0 - Galerkin (default)\n");
         printf("                        1 - non-Galerkin ParFlow operators\n");
         printf("                        2 - Galerkin, general operators\n");
         printf("  -relax <r>          : relaxation type\n");
         printf("                        0 - Jacobi\n");
         printf("                        1 - Weighted Jacobi (default)\n");
         printf("                        2 - R/B Gauss-Seidel\n");
         printf("                        3 - R/B Gauss-Seidel (nonsymmetric)\n");
         printf("  -skip <s>           : skip levels in PFMG (0 or 1)\n");
         printf("  -sym <s>            : symmetric storage (1) or not (0)\n");
         printf("  -vis                : save the solution for GLVis visualization\n");
         printf("\n");
      }

      if (print_usage)
      {
         MPI_Finalize();
         return (0);
      }
   }

   /* Convection produces non-symmetric matrices */
   if (optionB && sym)
      optionB = 0;

   /* Figure out the processor grid (N x N).  The local
      problem size is indicated by n (n x n). pi and pj
      indicate position in the processor grid. */
   N  = sqrt(num_procs);
   h  = 1.0 / (N*n-1);
   h2 = h*h;
   pj = myid / N;
   pi = myid - pj*N;
    

   /* Define the nodes owned by the current processor (each processor's
      piece of the global grid) */
   ilower[0] = pi*n;
   ilower[1] = pj*n;
   iupper[0] = ilower[0] + n-1;
   iupper[1] = ilower[1] + n-1;

   /* 1. Set up a grid */
   {
      /* Create an empty 2D grid object */
      HYPRE_StructGridCreate(MPI_COMM_WORLD, 2, &grid);

      /* Add a new box to the grid */
      HYPRE_StructGridSetExtents(grid, ilower, iupper);

      /* This is a collective call finalizing the grid assembly.
         The grid is now ``ready to be used'' */
      HYPRE_StructGridAssemble(grid);
   }

   /* 2. Define the discretization stencil */

      /* Define the geometry of the stencil */
      int offsets[5][2] = {{0,0}, {-1,0}, {1,0}, {0,-1}, {0,1}};

      /* Create an empty 2D, 5-pt stencil object */
      HYPRE_StructStencilCreate(2, 5, &stencil);

      /* Assign stencil entries */
      for (i = 0; i < 5; i++)
         HYPRE_StructStencilSetElement(stencil, i, offsets[i]);

   /* 3. Set up Struct Vectors for b and x */
   {
      double *values;

      /* Create an empty vector object */
      HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &b);
      HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &x);

      /* Indicate that the vector coefficients are ready to be set */
      HYPRE_StructVectorInitialize(b);
      HYPRE_StructVectorInitialize(x);

      values = calloc((n*n), sizeof(double));

      /* Set the values of b in left-to-right, bottom-to-top order */
      for (k = 0, j = 0; j < n; j++)
         for (i = 0; i < n; i++, k++)
            values[k] = h2 * Eval(F,i,j);
      HYPRE_StructVectorSetBoxValues(b, ilower, iupper, values);

      /* Set x = 0 */
      for (i = 0; i < (n*n); i ++)
         values[i] = 0.0;
      HYPRE_StructVectorSetBoxValues(x, ilower, iupper, values);

      free(values);

      /* Assembling is postponed since the vectors will be further modified */
   }

   /* 4. Set up a Struct Matrix */
   {
      /* Create an empty matrix object */
      HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &A);

      /* Indicate that the matrix coefficients are ready to be set */
      HYPRE_StructMatrixInitialize(A);

      /* Set the stencil values in the interior. Here we set the values
         at every node. We will modify the boundary nodes later. */

         int stencil_indices[5] = {0, 1, 2, 3, 4}; /* labels correspond
                                                      to the offsets */
         double *values;

         values = calloc(5*(n*n), sizeof(double));
       
       
       // variables for arithmetic mean
       double avgE,avgN,avgS,avgW;
       double varE, varN, varS, varW, varC;

       
       if (optionPSI == 0)
       {
           for (k = 0, j = 0; j < n; j++)
               for (i = 0; i < n; i++, k+=5)
               {
                   // arithmetic mean of variables
                   varW = Eval(K,i-0.5,j);
                   varE = Eval(K,i+0.5,j);
                   varS = Eval(K,i,j-0.5);
                   varN = Eval(K,i,j+0.5);
                   varC = Eval(K,i,j);
                   
                   
                   values[k+1] = - varW;
                   values[k+2] = - varE;
                   values[k+3] = - varS;
                   values[k+4] = - varN;
                   
                   values[k] = h2 * varC
                   + varW + varE + varS + varN;
               }
       }
       
       // arithmetic
       if(optionPSI == 1)
       {
         /* The order is left-to-right, bottom-to-top */
         for (k = 0, j = 0; j < n; j++)
            for (i = 0; i < n; i++, k+=5)
            {
                // arithmetic mean of variables
                varW = Eval(K,i-1,j);
                varE = Eval(K,i+1,j);
                varS = Eval(K,i,j-1);
                varN = Eval(K,i,j+1);
                varC = Eval(K,i,j);
                
                avgW = 0.5 * (varC + varW);
                avgE = 0.5 * (varC + varE);
                avgS = 0.5 * (varC + varS);
                avgN = 0.5 * (varC + varN);
                
                values[k+1] = - avgW;
                values[k+2] = - avgE;
                values[k+3] = - avgS;
                values[k+4] = - avgN;
                
                values[k] = h2 * varC
                + avgW + avgE + avgS + avgN;
            }
       }
       
       
       // harmonic
       if(optionPSI == 2)
       {
           /* The order is left-to-right, bottom-to-top */
           for (k = 0, j = 0; j < n; j++)
               for (i = 0; i < n; i++, k+=5)
               {
                   // arithmetic mean of variables
                   varW = Eval(K,i-1,j);
                   varE = Eval(K,i+1,j);
                   varS = Eval(K,i,j-1);
                   varN = Eval(K,i,j+1);
                   varC = Eval(K,i,j);
                   
                   avgW = 2. * varC * varW / (varC + varW);
                   avgE = 2. * varC * varE / (varC + varE);
                   avgS = 2. * varC * varS / (varC + varS);
                   avgN = 2. * varC * varN / (varC + varN);
                   
                   values[k+1] = - avgW;
                   values[k+2] = - avgE;
                   values[k+3] = - avgS;
                   values[k+4] = - avgN;
                   
                   values[k] = h2 * varC
                   + avgW + avgE + avgS + avgN;
               }
       }
   

         HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 5,
                                        stencil_indices, values);

         free(values);


   }

   /* 5. Set the boundary conditions, while eliminating the coefficients
         reaching ouside of the domain boundary. We must modify the matrix
         stencil and the corresponding rhs entries. */
   {
      int bc_ilower[2];
      int bc_iupper[2];

      int stencil_indices[5] = {0, 1, 2, 3, 4};
      double *values, *bvalues;

      int nentries;
      if (sym == 0)
         nentries = 5;
      else
         nentries = 3;

      values  = calloc(nentries*n, sizeof(double));
      bvalues = calloc(n, sizeof(double));

      /* The stencil at the boundary nodes is 1-0-0-0-0. Because
         we have I x_b = u_0; */
      for (i = 0; i < nentries*n; i += nentries)
      {
         values[i] = 1.0;
         for (j = 1; j < nentries; j++)
            values[i+j] = 0.0;
      }

      /* Processors at y = 0 */
      if (pj == 0)
      {
         bc_ilower[0] = pi*n;
         bc_ilower[1] = pj*n;

         bc_iupper[0] = bc_ilower[0] + n-1;
         bc_iupper[1] = bc_ilower[1];

         /* Modify the matrix */
         HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, nentries,
                                        stencil_indices, values);

         /* Put the boundary conditions in b */
         for (i = 0; i < n; i++)
            bvalues[i] = bcEval(U0,i,0);

         HYPRE_StructVectorSetBoxValues(b, bc_ilower, bc_iupper, bvalues);
      }

      /* Processors at y = 1 */
      if (pj == N-1)
      {
         bc_ilower[0] = pi*n;
         bc_ilower[1] = pj*n + n-1;

         bc_iupper[0] = bc_ilower[0] + n-1;
         bc_iupper[1] = bc_ilower[1];

         /* Modify the matrix */
         HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, nentries,
                                        stencil_indices, values);

         /* Put the boundary conditions in b */
         for (i = 0; i < n; i++)
            bvalues[i] = bcEval(U0,i,0);

         HYPRE_StructVectorSetBoxValues(b, bc_ilower, bc_iupper, bvalues);
      }

      /* Processors at x = 0 */
      if (pi == 0)
      {
         bc_ilower[0] = pi*n;
         bc_ilower[1] = pj*n;

         bc_iupper[0] = bc_ilower[0];
         bc_iupper[1] = bc_ilower[1] + n-1;

         /* Modify the matrix */
         HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, nentries,
                                        stencil_indices, values);

         /* Put the boundary conditions in b */
         for (j = 0; j < n; j++)
            bvalues[j] = bcEval(U0,0,j);

         HYPRE_StructVectorSetBoxValues(b, bc_ilower, bc_iupper, bvalues);
      }

      /* Processors at x = 1 */
      if (pi == N-1)
      {
         bc_ilower[0] = pi*n + n-1;
         bc_ilower[1] = pj*n;

         bc_iupper[0] = bc_ilower[0];
         bc_iupper[1] = bc_ilower[1] + n-1;

         /* Modify the matrix */
         HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, nentries,
                                        stencil_indices, values);

         /* Put the boundary conditions in b */
         for (j = 0; j < n; j++)
            bvalues[j] = bcEval(U0,0,j);

         HYPRE_StructVectorSetBoxValues(b, bc_ilower, bc_iupper, bvalues);
      }

      /* Recall that the system we are solving is:
         [A_ii 0; 0 I] [x_i ; x_b] = [b_i - A_ib u_0; u_0].
         This requires removing the connections between the interior
         and boundary nodes that we have set up when we set the
         5pt stencil at each node. We adjust for removing
         these connections by appropriately modifying the rhs.
         For the symm ordering scheme, just do the top and right
         boundary */

      /* Processors at y = 0, neighbors of boundary nodes */
      if (pj == 0)
      {
         bc_ilower[0] = pi*n;
         bc_ilower[1] = pj*n + 1;

         bc_iupper[0] = bc_ilower[0] + n-1;
         bc_iupper[1] = bc_ilower[1];

         stencil_indices[0] = 3;

         /* Modify the matrix */
         for (i = 0; i < n; i++)
            bvalues[i] = 0.0;

         if (sym == 0)
            HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, 1,
                                           stencil_indices, bvalues);

         /* Eliminate the boundary conditions in b */
         for (i = 0; i < n; i++)
            bvalues[i] = bcEval(U0,i,-1) * (bcEval(K,i,-0.5));

         if (pi == 0)
            bvalues[0] = 0.0;

         if (pi == N-1)
            bvalues[n-1] = 0.0;

         /* Note the use of AddToBoxValues (because we have already set values
            at these nodes) */
         HYPRE_StructVectorAddToBoxValues(b, bc_ilower, bc_iupper, bvalues);
      }

      /* Processors at x = 0, neighbors of boundary nodes */
      if (pi == 0)
      {
         bc_ilower[0] = pi*n + 1;
         bc_ilower[1] = pj*n;

         bc_iupper[0] = bc_ilower[0];
         bc_iupper[1] = bc_ilower[1] + n-1;

         stencil_indices[0] = 1;

         /* Modify the matrix */
         for (j = 0; j < n; j++)
            bvalues[j] = 0.0;

         if (sym == 0)
            HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, 1,
                                           stencil_indices, bvalues);

         /* Eliminate the boundary conditions in b */
         for (j = 0; j < n; j++)
            bvalues[j] = bcEval(U0,-1,j) * (bcEval(K,-0.5,j));

         if (pj == 0)
            bvalues[0] = 0.0;

         if (pj == N-1)
            bvalues[n-1] = 0.0;

         HYPRE_StructVectorAddToBoxValues(b, bc_ilower, bc_iupper, bvalues);
      }

      /* Processors at y = 1, neighbors of boundary nodes */
      if (pj == N-1)
      {
         bc_ilower[0] = pi*n;
         bc_ilower[1] = pj*n + (n-1) -1;

         bc_iupper[0] = bc_ilower[0] + n-1;
         bc_iupper[1] = bc_ilower[1];

         if (sym == 0)
            stencil_indices[0] = 4;
         else
            stencil_indices[0] = 2;

         /* Modify the matrix */
         for (i = 0; i < n; i++)
            bvalues[i] = 0.0;

         HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, 1,
                                        stencil_indices, bvalues);

         /* Eliminate the boundary conditions in b */
         for (i = 0; i < n; i++)
            bvalues[i] = bcEval(U0,i,1) * (bcEval(K,i,0.5));

         if (pi == 0)
            bvalues[0] = 0.0;

         if (pi == N-1)
            bvalues[n-1] = 0.0;

         HYPRE_StructVectorAddToBoxValues(b, bc_ilower, bc_iupper, bvalues);
      }

      /* Processors at x = 1, neighbors of boundary nodes */
      if (pi == N-1)
      {
         bc_ilower[0] = pi*n + (n-1) - 1;
         bc_ilower[1] = pj*n;

         bc_iupper[0] = bc_ilower[0];
         bc_iupper[1] = bc_ilower[1] + n-1;

         if (sym == 0)
            stencil_indices[0] = 2;
         else
            stencil_indices[0] = 1;

         /* Modify the matrix */
         for (j = 0; j < n; j++)
            bvalues[j] = 0.0;

         HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, 1,
                                        stencil_indices, bvalues);

         /* Eliminate the boundary conditions in b */
         for (j = 0; j < n; j++)
            bvalues[j] = bcEval(U0,1,j) * (bcEval(K,0.5,j) );

         if (pj == 0)
            bvalues[0] = 0.0;

         if (pj == N-1)
            bvalues[n-1] = 0.0;

         HYPRE_StructVectorAddToBoxValues(b, bc_ilower, bc_iupper, bvalues);
      }

      free(values);
      free(bvalues);
   }

   /* Finalize the vector and matrix assembly */
   HYPRE_StructMatrixAssemble(A);
   HYPRE_StructVectorAssemble(b);
   HYPRE_StructVectorAssemble(x);

   /* 6. Set up and use a solver */
    
      /* Start timing */
      time_index = hypre_InitializeTiming("SMG Setup");
      hypre_BeginTiming(time_index);

      /* Options and setup */
      HYPRE_StructSMGCreate(MPI_COMM_WORLD, &solver);
      HYPRE_StructSMGSetMemoryUse(solver, 0);
      HYPRE_StructSMGSetMaxIter(solver, 100);
      HYPRE_StructSMGSetTol(solver, 1.0e-08);
      HYPRE_StructSMGSetRelChange(solver, 0);
      //HYPRE_StructSMGSetNonZeroGuess(solver);
      HYPRE_StructSMGSetNumPreRelax(solver, n_pre);
      HYPRE_StructSMGSetNumPostRelax(solver, n_post);
       
     /* Logging must be on to get iterations and residual norm info below */
      HYPRE_StructSMGSetPrintLevel(solver, 1);
      HYPRE_StructSMGSetLogging(solver, 1);
       
       /* Setup and solve */
      HYPRE_StructSMGSetup(solver, A, b, x);


      /* Finalize current timing */
      hypre_EndTiming(time_index);
      hypre_PrintTiming("Setup phase times", MPI_COMM_WORLD);
      hypre_FinalizeTiming(time_index);
      hypre_ClearTiming();

      /* Start timing again */
      time_index = hypre_InitializeTiming("SMG Solve");
      hypre_BeginTiming(time_index);

      /* Solve */
      HYPRE_StructSMGSolve(solver, A, b, x);

      /* Finalize current timing */
      hypre_EndTiming(time_index);
      hypre_PrintTiming("Solve phase times", MPI_COMM_WORLD);
      hypre_FinalizeTiming(time_index);
      hypre_ClearTiming();

      /* Get info and release memory */
      HYPRE_StructSMGGetNumIterations(solver, &num_iterations);
      HYPRE_StructSMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
    HYPRE_StructSMGDestroy(solver);
    
    
    /* ====== Save output to be read-in in Matlab ======= */
    if(myid == 0)
    {
        FILE *file;
        char filename[255];
        
        int nvalues = n*n;
        double *values = calloc(nvalues, sizeof(double));
        
        /* ==== solution =====*/
        HYPRE_StructVectorGetBoxValues(x, ilower, iupper, values);
        
        
        if(optionK==0)
            sprintf(filename, "%s", "exJanaSquare_u.dat");
        else
            sprintf(filename, "%s", "exJanaSphere_u.dat");

        if ((file = fopen(filename, "w")) == NULL)
        {
            printf("Error: can't open output file %s\n", filename);
            MPI_Finalize();
            exit(1);
        }
        
        int i,j,k;
        k = 0;
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
                fprintf(file, " %.6e", values[k++]);
            
            fprintf(file, "\n");
        }
        
        fflush(file);
        fclose(file);
        
        
        /* ======== save rhs = f*psi  ========*/
        HYPRE_StructVectorGetBoxValues(b, ilower, iupper, values);
        
        
        if(optionK==0)
            sprintf(filename, "%s", "exJanaSquare_rhs.dat");
        else
            sprintf(filename, "%s", "exJanaSphere_rhs.dat");
        
        if ((file = fopen(filename, "w")) == NULL)
        {
            printf("Error: can't open output file %s\n", filename);
            MPI_Finalize();
            exit(1);
        }
        
        /* save solution with global unknown numbers */
        k = 0;
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
                fprintf(file, " %.6e", values[k++]);
            
            fprintf(file, "\n");
        }
        
        fflush(file);
        fclose(file);
        //        free(values);
        
        /* ========= recompute and save psi =========*/
        
        if(optionK==0)
            sprintf(filename, "%s", "exJanaSquare_psi.dat");
        else
            sprintf(filename, "%s", "exJanaSphere_psi.dat");
        
        if ((file = fopen(filename, "w")) == NULL)
        {
            printf("Error: can't open output file %s\n", filename);
            MPI_Finalize();
            exit(1);
        }
        
        /* save solution with global unknown numbers */
        
        for (k = 0, j = 0; j < n; j++)
            for (i = 0; i < n; i++, k++)
                values[k] = Eval(K,i,j);
        
        k=0;
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
                fprintf(file, " %.6e", values[k++]);
            
            fprintf(file, "\n");
        }
        
        fflush(file);
        fclose(file);
        free(values);
    }
    

   if (myid == 0)
   {
      printf("\n");
      printf("Iterations = %d\n", num_iterations);
      printf("Final Relative Residual Norm = %e\n", final_res_norm);
      printf("\n");
   }


   /* Free memory */
   HYPRE_StructGridDestroy(grid);
   HYPRE_StructStencilDestroy(stencil);
   HYPRE_StructMatrixDestroy(A);
   HYPRE_StructVectorDestroy(b);
   HYPRE_StructVectorDestroy(x);


   /* Finalize MPI */
   MPI_Finalize();

   return (0);
}
