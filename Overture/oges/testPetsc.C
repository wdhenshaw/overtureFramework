//===============================================================================
//  TEST ROUTINE FOR PETSC
//
// This file was create to understand how version 13.8.2 deals with re-use or not
// of a preconditioner -- the interface changed at version 3.5 
//  
// Usage: 
//
// Version 1: Dec 2022 (WDH)
//==============================================================================

#include "mpi.h"
#include "Overture.h"
#include "Oges.h"

#include <petscksp.h>

// #include "CompositeGridOperators.h"
// #include "SquareMapping.h"
// #include "AnnulusMapping.h"
// #include "OGPolyFunction.h"
// #include "OGTrigFunction.h"
// #include "SparseRep.h"
// #include "PlotStuff.h"


// #define ForBoundary(side,axis)   for( axis=0; axis<mg.numberOfDimensions(); axis++ ) \
//                                  for( side=0; side<=1; side++ )
// bool measureCPU=TRUE;
// real
// CPU()
// // In this version of getCPU we can turn off the timing
// {
//   if( measureCPU )
//     return getCPU();
//   else
//     return 0;
// }


static char help[] = "Solves two linear systems in parallel with KSP.  The code\n\
illustrates repeated solution of linear systems with the same preconditioner\n\
method but different matrices (having the same nonzero structure).  The code\n\
also uses multiple profiling stages.  Input arguments are\n\
  -m <size> : problem size\n\
  -mat_nonsym : use nonsymmetric matrix (default is symmetric)\n\n";

int 
main(int argc, char *args[])
{
  printF("testPetsc start...\n");

  Overture::start(argc,args);  // initialize Overture

  // const int maxNumberOfGridsToTest=3;
  // int numberOfGridsToTest=maxNumberOfGridsToTest;
  // aString gridName[maxNumberOfGridsToTest] =   { "square5", "cic", "sib" };
  // aString commandFileName="";
  // bool plotOption=TRUE;
  
  // This macro will initialize the PETSc solver if OVERTURE_USE_PETSC is defined.
  // INIT_PETSC_SOLVER();

  Vec         x, b, u; /* approx solution, RHS, exact solution */
  Mat         A;       /* linear system matrix */
  KSP         ksp;     /* linear solver context */
  PetscReal   norm;    /* norm of solution error */
  PetscInt    i, j, Ii, J, Istart, Iend, m = 8, n = 7, its;
  PetscBool   flg;
  PetscScalar v;

  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &args, (char *)0, help));
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-m", &m, NULL));
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-n", &n, NULL));
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         Compute the matrix and right-hand-side vector that define
         the linear system, Ax = b.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create parallel matrix, specifying only its global dimensions.
     When using MatCreate(), the matrix format can be specified at
     runtime. Also, the parallel partitioning of the matrix is
     determined by PETSc at runtime.

     Performance tuning note:  For problems of substantial size,
     preallocation of matrix memory is crucial for attaining good
     performance. See the matrix chapter of the users manual for details.
  */
  PetscCall(MatCreate(PETSC_COMM_WORLD, &A));
  PetscCall(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, m * n, m * n));
  PetscCall(MatSetFromOptions(A));
  PetscCall(MatMPIAIJSetPreallocation(A, 5, NULL, 5, NULL));
  PetscCall(MatSeqAIJSetPreallocation(A, 5, NULL));
  PetscCall(MatSeqSBAIJSetPreallocation(A, 1, 5, NULL));
  PetscCall(MatMPISBAIJSetPreallocation(A, 1, 5, NULL, 5, NULL));
  PetscCall(MatMPISELLSetPreallocation(A, 5, NULL, 5, NULL));
  PetscCall(MatSeqSELLSetPreallocation(A, 5, NULL));

  /*
     Currently, all PETSc parallel matrix formats are partitioned by
     contiguous chunks of rows across the processors.  Determine which
     rows of the matrix are locally owned.
  */
  PetscCall(MatGetOwnershipRange(A, &Istart, &Iend));

  /*
     Set matrix elements for the 2-D, five-point stencil in parallel.
      - Each processor needs to insert only elements that it owns
        locally (but any non-local elements will be sent to the
        appropriate processor during matrix assembly).
      - Always specify global rows and columns of matrix entries.

     Note: this uses the less common natural ordering that orders first
     all the unknowns for x = h then for x = 2h etc; Hence you see J = Ii +- n
     instead of J = I +- m as you might expect. The more standard ordering
     would first do all variables for y = h, then y = 2h etc.

   */
  for (Ii = Istart; Ii < Iend; Ii++) {
    v = -1.0;
    i = Ii / n;
    j = Ii - i * n;
    if (i > 0) {
      J = Ii - n;
      PetscCall(MatSetValues(A, 1, &Ii, 1, &J, &v, ADD_VALUES));
    }
    if (i < m - 1) {
      J = Ii + n;
      PetscCall(MatSetValues(A, 1, &Ii, 1, &J, &v, ADD_VALUES));
    }
    if (j > 0) {
      J = Ii - 1;
      PetscCall(MatSetValues(A, 1, &Ii, 1, &J, &v, ADD_VALUES));
    }
    if (j < n - 1) {
      J = Ii + 1;
      PetscCall(MatSetValues(A, 1, &Ii, 1, &J, &v, ADD_VALUES));
    }
    v = 4.0;
    PetscCall(MatSetValues(A, 1, &Ii, 1, &Ii, &v, ADD_VALUES));
  }

  /*
     Assemble matrix, using the 2-step process:
       MatAssemblyBegin(), MatAssemblyEnd()
     Computations can be done while messages are in transition
     by placing code between these two statements.
  */
  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));

  /* A is symmetric. Set symmetric flag to enable ICC/Cholesky preconditioner */
  PetscCall(MatSetOption(A, MAT_SYMMETRIC, PETSC_TRUE));

  /*
     Create parallel vectors.
      - We form 1 vector from scratch and then duplicate as needed.
      - When using VecCreate(), VecSetSizes and VecSetFromOptions()
        in this example, we specify only the
        vector's global dimension; the parallel partitioning is determined
        at runtime.
      - When solving a linear system, the vectors and matrices MUST
        be partitioned accordingly.  PETSc automatically generates
        appropriately partitioned matrices and vectors when MatCreate()
        and VecCreate() are used with the same communicator.
      - The user can alternatively specify the local vector and matrix
        dimensions when more sophisticated partitioning is needed
        (replacing the PETSC_DECIDE argument in the VecSetSizes() statement
        below).
  */
  PetscCall(VecCreate(PETSC_COMM_WORLD, &u));
  PetscCall(VecSetSizes(u, PETSC_DECIDE, m * n));
  PetscCall(VecSetFromOptions(u));
  PetscCall(VecDuplicate(u, &b));
  PetscCall(VecDuplicate(b, &x));

  /*
     Set exact solution; then compute right-hand-side vector.
     By default we use an exact solution of a vector with all
     elements of 1.0;
  */
  PetscCall(VecSet(u, 1.0));
  PetscCall(MatMult(A, u, b));

  /*
     View the exact solution vector if desired
  */
  flg = PETSC_FALSE;
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-view_exact_sol", &flg, NULL));
  if (flg) PetscCall(VecView(u, PETSC_VIEWER_STDOUT_WORLD));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the linear solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));

  printF("Set reuse PC true\n");
  PetscCall(KSPSetReusePreconditioner(ksp,PETSC_TRUE));  

  /*
     Set operators. Here the matrix that defines the linear system
     also serves as the preconditioning matrix.
  */
  PetscCall(KSPSetOperators(ksp, A, A));

  /*
     Set linear solver defaults for this problem (optional).
     - By extracting the KSP and PC contexts from the KSP context,
       we can then directly call any KSP and PC routines to set
       various options.
     - The following two statements are optional; all of these
       parameters could alternatively be specified at runtime via
       KSPSetFromOptions().  All of these defaults can be
       overridden at runtime, as indicated below.
  */
  PetscCall(KSPSetTolerances(ksp, 1.e-2 / ((m + 1) * (n + 1)), 1.e-50, PETSC_DEFAULT, PETSC_DEFAULT));

  /*
    Set runtime options, e.g.,
        -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    These options will override those specified above as long as
    KSPSetFromOptions() is called _after_ any other customization
    routines.
  */
  PetscCall(KSPSetFromOptions(ksp));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Real cpu0= getCPU();
  PetscCall(KSPSolve(ksp, b, x));
  Real cpuSolve1a = getCPU()-cpu0;
  printF("Time to solve (1a) = %9.3e (s)\n",cpuSolve1a);


  cpu0= getCPU();
  PetscCall(KSPSolve(ksp, b, x));
  Real cpuSolve1b = getCPU()-cpu0;
  printF("Time to solve (1b) = %9.3e (s), setup=%9.3e(s)\n",cpuSolve1b,cpuSolve1a-cpuSolve1b);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Check the solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(VecAXPY(x, -1.0, u));
  PetscCall(VecNorm(x, NORM_2, &norm));
  PetscCall(KSPGetIterationNumber(ksp, &its));

  /*
     Print convergence information.  PetscPrintf() produces a single
     print statement from all processes that share a communicator.
     An alternative is PetscFPrintf(), which prints to a file.
  */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Norm of error %g iterations %" PetscInt_FMT "\n", (double)norm, its));



 // ******** CHANGE ENTRIES ***************
 // printF("Set reuse PC true\n");
 // PetscCall(KSPSetReusePreconditioner(ksp,PETSC_TRUE));

 for (Ii = Istart; Ii < Iend; Ii++) {
    v = -1.001;
    i = Ii / n;
    j = Ii - i * n;
    if (i > 0) {
      J = Ii - n;
      PetscCall(MatSetValues(A, 1, &Ii, 1, &J, &v, ADD_VALUES));
    }
    if (i < m - 1) {
      J = Ii + n;
      PetscCall(MatSetValues(A, 1, &Ii, 1, &J, &v, ADD_VALUES));
    }
    if (j > 0) {
      J = Ii - 1;
      PetscCall(MatSetValues(A, 1, &Ii, 1, &J, &v, ADD_VALUES));
    }
    if (j < n - 1) {
      J = Ii + 1;
      PetscCall(MatSetValues(A, 1, &Ii, 1, &J, &v, ADD_VALUES));
    }
    v = 4.004;
    PetscCall(MatSetValues(A, 1, &Ii, 1, &Ii, &v, ADD_VALUES));
  }


  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));

  PetscCall(VecSet(u, 1.0));
  PetscCall(MatMult(A, u, b));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  cpu0= getCPU();
  PetscCall(KSPSolve(ksp, b, x));
  Real cpuSolve2a = getCPU()-cpu0;
  printF("Time to solve (2a) = %9.3e (s)\n",cpuSolve2a);


  cpu0= getCPU();
  PetscCall(KSPSolve(ksp, b, x));
  Real cpuSolve2b = getCPU()-cpu0;
  printF("Time to solve (2b) = %9.3e (s), setup=%9.3e(s)\n",cpuSolve2b,cpuSolve2a-cpuSolve2b);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Check the solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(VecAXPY(x, -1.0, u));
  PetscCall(VecNorm(x, NORM_2, &norm));
  PetscCall(KSPGetIterationNumber(ksp, &its));

  /*
     Print convergence information.  PetscPrintf() produces a single
     print statement from all processes that share a communicator.
     An alternative is PetscFPrintf(), which prints to a file.
  */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Norm of error %g iterations %" PetscInt_FMT "\n", (double)norm, its));




  /*
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
  */
  PetscCall(KSPDestroy(&ksp));
  PetscCall(VecDestroy(&u));
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&b));
  PetscCall(MatDestroy(&A));

  /*
     Always call PetscFinalize() before exiting a program.  This routine
       - finalizes the PETSc libraries as well as MPI
       - provides summary and diagnostic information if certain runtime
         options are chosen (e.g., -log_view).
  */
  PetscCall(PetscFinalize());



//   KSP         ksp;         /* linear solver context */
//   Mat         C;           /* matrix */
//   Vec         x, u, b;     /* approx solution, RHS, exact solution */
//   PetscReal   norm, bnorm; /* norm of solution error */
//   PetscScalar v, none = -1.0;
//   PetscInt    Ii, J, ldim, low, high, iglobal, Istart, Iend;
//   PetscInt    i, j, m = 3, n = 2, its;
//   PetscMPIInt size, rank;
//   PetscBool   mat_nonsymmetric = PETSC_FALSE;
//   PetscBool   testnewC = PETSC_FALSE, testscaledMat = PETSC_FALSE;
// #if defined(PETSC_USE_LOG)
//   PetscLogStage stages[2];
// #endif
  
//   PetscFunctionBeginUser;
//   PetscCall(PetscInitialize(&argc, &args, (char *)0, help));
//   PetscCall(PetscOptionsGetInt(NULL, NULL, "-m", &m, NULL));
//   PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
//   PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
//   n = 2 * size;

//   printf("m=%d, n=%d, num-proc=size=%d\n",m,n,size);

//   /*
//      Set flag if we are doing a nonsymmetric problem; the default is symmetric.
//   */
//   PetscCall(PetscOptionsGetBool(NULL, NULL, "-mat_nonsym", &mat_nonsymmetric, NULL));

//   PetscCall(PetscOptionsGetBool(NULL, NULL, "-test_scaledMat", &testscaledMat, NULL));

//   /*
//      Register two stages for separate profiling of the two linear solves.
//      Use the runtime option -log_view for a printout of performance
//      statistics at the program's conlusion.
//   */
//   PetscCall(PetscLogStageRegister("Original Solve", &stages[0]));
//   PetscCall(PetscLogStageRegister("Second Solve", &stages[1]));

//   /* -------------- Stage 0: Solve Original System ---------------------- */
//   /*
//      Indicate to PETSc profiling that we're beginning the first stage
//   */
//   PetscCall(PetscLogStagePush(stages[0]));

//   /*
//      Create parallel matrix, specifying only its global dimensions.
//      When using MatCreate(), the matrix format can be specified at
//      runtime. Also, the parallel partitioning of the matrix is
//      determined by PETSc at runtime.
//   */
//   PetscCall(MatCreate(PETSC_COMM_WORLD, &C));
//   PetscCall(MatSetSizes(C, PETSC_DECIDE, PETSC_DECIDE, m * n, m * n));
//   PetscCall(MatSetFromOptions(C));
//   PetscCall(MatSetUp(C));

//   /*
//      Currently, all PETSc parallel matrix formats are partitioned by
//      contiguous chunks of rows across the processors.  Determine which
//      rows of the matrix are locally owned.
//   */
//   PetscCall(MatGetOwnershipRange(C, &Istart, &Iend));

//   /*
//      Set matrix entries matrix in parallel.
//       - Each processor needs to insert only elements that it owns
//         locally (but any non-local elements will be sent to the
//         appropriate processor during matrix assembly).
//       - Always specify global row and columns of matrix entries.
//   */
//   for (Ii = Istart; Ii < Iend; Ii++) {
//     v = -1.0;
//     i = Ii / n;
//     j = Ii - i * n;
//     if (i > 0) {
//       J = Ii - n;
//       PetscCall(MatSetValues(C, 1, &Ii, 1, &J, &v, ADD_VALUES));
//     }
//     if (i < m - 1) {
//       J = Ii + n;
//       PetscCall(MatSetValues(C, 1, &Ii, 1, &J, &v, ADD_VALUES));
//     }
//     if (j > 0) {
//       J = Ii - 1;
//       PetscCall(MatSetValues(C, 1, &Ii, 1, &J, &v, ADD_VALUES));
//     }
//     if (j < n - 1) {
//       J = Ii + 1;
//       PetscCall(MatSetValues(C, 1, &Ii, 1, &J, &v, ADD_VALUES));
//     }
//     v = 4.0;
//     PetscCall(MatSetValues(C, 1, &Ii, 1, &Ii, &v, ADD_VALUES));
//   }

//   /*
//      Make the matrix nonsymmetric if desired
//   */
//   if (mat_nonsymmetric) {
//     for (Ii = Istart; Ii < Iend; Ii++) {
//       v = -1.5;
//       i = Ii / n;
//       if (i > 1) {
//         J = Ii - n - 1;
//         PetscCall(MatSetValues(C, 1, &Ii, 1, &J, &v, ADD_VALUES));
//       }
//     }
//   } else {
//     PetscCall(MatSetOption(C, MAT_SYMMETRIC, PETSC_TRUE));
//     PetscCall(MatSetOption(C, MAT_SYMMETRY_ETERNAL, PETSC_TRUE));
//   }

//   /*
//      Assemble matrix, using the 2-step process:
//        MatAssemblyBegin(), MatAssemblyEnd()
//      Computations can be done while messages are in transition
//      by placing code between these two statements.
//   */
//   PetscCall(MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY));
//   PetscCall(MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY));

//   /*
//      Create parallel vectors.
//       - When using VecSetSizes(), we specify only the vector's global
//         dimension; the parallel partitioning is determined at runtime.
//       - Note: We form 1 vector from scratch and then duplicate as needed.
//   */
//   PetscCall(VecCreate(PETSC_COMM_WORLD, &u));
//   PetscCall(VecSetSizes(u, PETSC_DECIDE, m * n));
//   PetscCall(VecSetFromOptions(u));
//   PetscCall(VecDuplicate(u, &b));
//   PetscCall(VecDuplicate(b, &x));

//   /*
//      Currently, all parallel PETSc vectors are partitioned by
//      contiguous chunks across the processors.  Determine which
//      range of entries are locally owned.
//   */
//   PetscCall(VecGetOwnershipRange(x, &low, &high));

//   /*
//     Set elements within the exact solution vector in parallel.
//      - Each processor needs to insert only elements that it owns
//        locally (but any non-local entries will be sent to the
//        appropriate processor during vector assembly).
//      - Always specify global locations of vector entries.
//   */
//   PetscCall(VecGetLocalSize(x, &ldim));
//   for (i = 0; i < ldim; i++) {
//     iglobal = i + low;
//     v       = (PetscScalar)(i + 100 * rank);
//     PetscCall(VecSetValues(u, 1, &iglobal, &v, INSERT_VALUES));
//   }

//   /*
//      Assemble vector, using the 2-step process:
//        VecAssemblyBegin(), VecAssemblyEnd()
//      Computations can be done while messages are in transition,
//      by placing code between these two statements.
//   */
//   PetscCall(VecAssemblyBegin(u));
//   PetscCall(VecAssemblyEnd(u));

//   /*
//      Compute right-hand-side vector
//   */
//   PetscCall(MatMult(C, u, b));

//   /*
//     Create linear solver context
//   */
//   PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));

//   /*
//      Set operators. Here the matrix that defines the linear system
//      also serves as the preconditioning matrix.
//   */
//   PetscCall(KSPSetOperators(ksp, C, C));

//   /*
//      Set runtime options (e.g., -ksp_type <type> -pc_type <type>)
//   */
//   PetscCall(KSPSetFromOptions(ksp));

//   /*
//      Solve linear system.  Here we explicitly call KSPSetUp() for more
//      detailed performance monitoring of certain preconditioners, such
//      as ICC and ILU.  This call is optional, as KSPSetUp() will
//      automatically be called within KSPSolve() if it hasn't been
//      called already.
//   */
//   PetscCall(KSPSetUp(ksp));
//   PetscCall(KSPSolve(ksp, b, x));

//   /*
//      Check the residual
//   */
//   PetscCall(VecAXPY(x, none, u));
//   PetscCall(VecNorm(x, NORM_2, &norm));
//   PetscCall(VecNorm(b, NORM_2, &bnorm));
//   PetscCall(KSPGetIterationNumber(ksp, &its));
//   if (!testscaledMat || norm / bnorm > 1.e-5) PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Relative norm of the residual %g, Iterations %" PetscInt_FMT "\n", (double)norm / (double)bnorm, its));

//   /* -------------- Stage 1: Solve Second System ---------------------- */
//   /*
//      Solve another linear system with the same method.  We reuse the KSP
//      context, matrix and vector data structures, and hence save the
//      overhead of creating new ones.

//      Indicate to PETSc profiling that we're concluding the first
//      stage with PetscLogStagePop(), and beginning the second stage with
//      PetscLogStagePush().
//   */
//   printf("Solve again, adjust entries in the matrix...\n");

//   PetscCall(PetscLogStagePop());
//   PetscCall(PetscLogStagePush(stages[1]));

//   /*
//      Initialize all matrix entries to zero.  MatZeroEntries() retains the
//      nonzero structure of the matrix for sparse formats.
//   */
//   PetscCall(MatZeroEntries(C));

//   /*
//      Assemble matrix again.  Note that we retain the same matrix data
//      structure and the same nonzero pattern; we just change the values
//      of the matrix entries.
//   */
//   for (i = 0; i < m; i++) {
//     for (j = 2 * rank; j < 2 * rank + 2; j++) {
//       v  = -1.0;
//       Ii = j + n * i;
//       if (i > 0) {
//         J = Ii - n;
//         PetscCall(MatSetValues(C, 1, &Ii, 1, &J, &v, ADD_VALUES));
//       }
//       if (i < m - 1) {
//         J = Ii + n;
//         PetscCall(MatSetValues(C, 1, &Ii, 1, &J, &v, ADD_VALUES));
//       }
//       if (j > 0) {
//         J = Ii - 1;
//         PetscCall(MatSetValues(C, 1, &Ii, 1, &J, &v, ADD_VALUES));
//       }
//       if (j < n - 1) {
//         J = Ii + 1;
//         PetscCall(MatSetValues(C, 1, &Ii, 1, &J, &v, ADD_VALUES));
//       }
//       v = 6.0;
//       PetscCall(MatSetValues(C, 1, &Ii, 1, &Ii, &v, ADD_VALUES));
//     }
//   }
//   if (mat_nonsymmetric) {
//     for (Ii = Istart; Ii < Iend; Ii++) {
//       v = -1.5;
//       i = Ii / n;
//       if (i > 1) {
//         J = Ii - n - 1;
//         PetscCall(MatSetValues(C, 1, &Ii, 1, &J, &v, ADD_VALUES));
//       }
//     }
//   }
//   PetscCall(MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY));
//   PetscCall(MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY));

//   if (testscaledMat) {
//     PetscRandom rctx;

//     /* Scale a(0,0) and a(M-1,M-1) */
//     if (rank == 0) {
//       v  = 6.0 * 0.00001;
//       Ii = 0;
//       J  = 0;
//       PetscCall(MatSetValues(C, 1, &Ii, 1, &J, &v, INSERT_VALUES));
//     } else if (rank == size - 1) {
//       v  = 6.0 * 0.00001;
//       Ii = m * n - 1;
//       J  = m * n - 1;
//       PetscCall(MatSetValues(C, 1, &Ii, 1, &J, &v, INSERT_VALUES));
//     }
//     PetscCall(MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY));
//     PetscCall(MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY));

//     /* Compute a new right-hand-side vector */
//     PetscCall(VecDestroy(&u));
//     PetscCall(VecCreate(PETSC_COMM_WORLD, &u));
//     PetscCall(VecSetSizes(u, PETSC_DECIDE, m * n));
//     PetscCall(VecSetFromOptions(u));

//     PetscCall(PetscRandomCreate(PETSC_COMM_WORLD, &rctx));
//     PetscCall(PetscRandomSetFromOptions(rctx));
//     PetscCall(VecSetRandom(u, rctx));
//     PetscCall(PetscRandomDestroy(&rctx));
//     PetscCall(VecAssemblyBegin(u));
//     PetscCall(VecAssemblyEnd(u));
//   }

//   PetscCall(PetscOptionsGetBool(NULL, NULL, "-test_newMat", &testnewC, NULL));
//   if (testnewC) {
//     /*
//      User may use a new matrix C with same nonzero pattern, e.g.
//       ./ex5 -ksp_monitor -mat_type sbaij -pc_type cholesky -pc_factor_mat_solver_type mumps -test_newMat
//     */
//     Mat Ctmp;
//     PetscCall(MatDuplicate(C, MAT_COPY_VALUES, &Ctmp));
//     PetscCall(MatDestroy(&C));
//     PetscCall(MatDuplicate(Ctmp, MAT_COPY_VALUES, &C));
//     PetscCall(MatDestroy(&Ctmp));
//   }

//   PetscCall(MatMult(C, u, b));

//   /*
//      Set operators. Here the matrix that defines the linear system
//      also serves as the preconditioning matrix.
//   */
//   PetscCall(KSPSetOperators(ksp, C, C));

//   /*
//      Solve linear system
//   */
//   PetscCall(KSPSetUp(ksp));
//   PetscCall(KSPSolve(ksp, b, x));

//   /*
//      Check the residual
//   */
//   PetscCall(VecAXPY(x, none, u));
//   PetscCall(VecNorm(x, NORM_2, &norm));
//   PetscCall(KSPGetIterationNumber(ksp, &its));
//   if (!testscaledMat || norm / bnorm > PETSC_SMALL) PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Relative norm of the residual %g, Iterations %" PetscInt_FMT "\n", (double)norm / (double)bnorm, its));

//   /*
//      Free work space.  All PETSc objects should be destroyed when they
//      are no longer needed.
//   */
//   PetscCall(KSPDestroy(&ksp));
//   PetscCall(VecDestroy(&u));
//   PetscCall(VecDestroy(&x));
//   PetscCall(VecDestroy(&b));
//   PetscCall(MatDestroy(&C));

//   /*
//      Indicate to PETSc profiling that we're concluding the second stage
//   */
//   PetscCall(PetscLogStagePop());

//   PetscCall(PetscFinalize());


  Overture::finish();   
  printF("testPetsc done.\n");

  return(0);
}


