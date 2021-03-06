/**********************************************************************
 * Assignment 3 part 3 - Parallelization
 * By Alessandro Piccolo & Anton Sjöberg
 *
 **********************************************************************/

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void updateBoundary(double *A, int jmax, int imax, int n, int m);
void matrix_to_txt(double *matrix, int jmax, int imax, int m);
void print_matrix(double *matrix, int n, int m);

int main(int argc, char *argv[]) {
  int jmax = atoi(argv[1]);                 /* Inner number of rows           */
  int imax = atoi(argv[2]);                 /* Inner number of columns        */
  double omega = atoi(argv[3]);             /* Omega constant value 1 1.7 1.9 */
  int threads = atoi(argv[4]);              /* Give the number of threads     */ 

  int meshsize_i = imax/imax;               /* Size of block?                 */
  int meshsize_j = jmax/jmax;               /* Size of block?                 */
  int meshsize_i2 = meshsize_i*meshsize_i;  /* Meshsize of i                  */
  int meshsize_j2 = meshsize_j*meshsize_j;  /* Meshsize of j                  */
  int n = jmax + 2;                         /* Number of rows of the grid     */
  int m = imax + 2;                         /* Number of columns of the grid  */
  int total_inner_grid_elem = imax*jmax;    /* Tot inner grid elem            */
  int total_grid_elem = n*m;                /* Tot grid elem                  */

  double max_error = 0.008;                 /* Max error of L2 L2-norm        */
  float L2_norm_residual = max_error+1;     /* L2-norm of the residual        */ 
  int iterations = 1;                       /* Numbers of iterations          */
  int max_iterations = 1000;                /* Max iterations                 */
  double avrg, sum, time;                   /* Sum & avrg of elements & time  */
  int i,j;

  double meshConst = omega/(2/meshsize_i2 + 2/meshsize_j2);
  double omega_1 = 1-omega;

  double *A, *residual, *function_output;
  A = (double *)calloc(total_grid_elem, sizeof(double));
  residual = (double *)malloc(total_inner_grid_elem*sizeof(double));
  function_output = (double *)malloc(imax*sizeof(double));
  
  printf("...: Poisson parallel calculation :...\n");

  /* Fill inner matrix with ones */ 
  for(j = 1; j <= jmax; j++) { /* Row */
    for (i = 1; i <= imax; i++) { /* Column */
      A[j*m + i] = drand48();
    }
  }

  /* 1. Values of the inner points touching the bounday to the boundary layer */
  updateBoundary(A, jmax, imax, n, m);
  
  /* Start time */
  time = omp_get_wtime();

  /* Our right hand side function, gives sin outputs */
#pragma omp parallel num_threads(threads) 
{
  #pragma omp for
  for (i = 0; i < imax; i++) {
    function_output[i] = sin(3.14*2*(i+0.5));
  }
  
  while (L2_norm_residual > max_error && iterations < max_iterations) {
    /* 2. Perform an iteration of SOR on the inner points red points */ 
    #pragma omp for private(j,i)
    for(j = 1; j <= jmax; j++) { /* row */
      if (j%2 != 0) /* odd rows */
        for (i = 1; i <= imax; i+=2) { /* column */
          A[j*m+i] = omega_1*A[j*m+i] + meshConst * (
          (A[(j+1)*m+i] + A[(j-1)*m+i])/meshsize_i2 + 
          (A[j*m+(i+1)]+A[j*m+(i-1)])/meshsize_j2 
          - function_output[i-1] );
      }
      else { /* even rows */
        for (i = 2; i <= imax; i+=2) { /* column */
          A[j*m+i] = omega_1*A[j*m+i] + meshConst * (
          (A[(j+1)*m+i] + A[(j-1)*m+i])/meshsize_i2 + 
          (A[j*m+(i+1)]+A[j*m+(i-1)])/meshsize_j2 
          - function_output[i-1] );
        }
      }
    }

    /* 2. Perform an iteration of SOR on the inner points black points */ 
    #pragma omp for private(j,i)
    for(j = 1; j <= jmax; j++) { /* row */ 
      if (j%2 == 0) /* even rows */
        for (i = 1; i <= imax; i+=2) { /* column */
          A[j*m+i] = omega_1*A[j*m+i] + meshConst * (
          (A[(j+1)*m+i] + A[(j-1)*m+i])/meshsize_i2 + 
          (A[j*m+(i+1)]+A[j*m+(i-1)])/meshsize_j2 
          - function_output[i-1] );
      }
      else { /* odd rows */
        for (i = 2; i <= imax; i+=2) { /* column */
          A[j*m+i] = omega_1*A[j*m+i] + meshConst * (
          (A[(j+1)*m+i] + A[(j-1)*m+i])/meshsize_i2 + 
          (A[j*m+(i+1)]+A[j*m+(i-1)])/meshsize_j2 
          - function_output[i-1] );
        }
      }
    }

    /* 3. Update BC, same as in (1) */
    #pragma omp for nowait private(j) 
    for(j = 1; j <= imax; j++) {  /* Top and bottom part of the boundary */
      A[j] = A[m+j];
      A[(n-1)*m+j] = A[(n-2)*m+j];
    }
    #pragma omp for nowait private(i)
    for(i = 1; i<= jmax; i++) {   /* Left and right part of the boundary */
      A[i*m] = A[i*m+1];
      A[i*m+m-1] = A[i*m+m-2];
    }

    /* Normalize u every 10th iteration to be around zero */
    if (iterations % 10 == 0) {
      #pragma omp single
      sum = 0.0;
      #pragma omp for reduction(+:sum) private(j)
      for (j = 0; j < total_grid_elem; j++) {
        sum += A[j];
      }
      #pragma omp single 
      avrg = sum/ total_grid_elem;
      #pragma omp for private(j)
      for (j = 0; j < total_grid_elem; j++) {
        A[j] -= avrg;
      }
    }

    #pragma omp single nowait 
    L2_norm_residual = 0.0;

    /* 4. Compute the residual of the discrete equation for the inner grid */
    #pragma omp for private(j,i)
    for(j = 1; j <= jmax; j++) {
      //#pragma omp for private(i)
      for (i = 1; i <= imax; i++) {
        residual[(j-1)*imax+(i-1)] = function_output[i-1] - 
        (A[(j+1)*m+i]-2*A[j*m+i]+A[(j-1)*m+i])/meshsize_i2 -
        (A[j*m+(i+1)]-2*A[j*m+i]+A[j*m+(i-1)])/meshsize_j2;
      }
    }

    /* 4. L2-norm */
    #pragma omp for reduction(+:L2_norm_residual) private(j)
    for(j = 0; j < total_inner_grid_elem; j++) {
      L2_norm_residual += residual[j]*residual[j];
    }
    #pragma omp single 
    {
    L2_norm_residual = 1/sqrt(imax+jmax)*sqrt(L2_norm_residual);
    //printf("Iteration: %d, residual norm = %lf\n",iterations,L2_norm_residual);
    iterations += 1;
    }
  }
} /* End parallel */

  time = omp_get_wtime() - time;
  printf("\nElapsed time: %f \n",time);
  printf("Total number of iterations %d, final residual norm %lf\n",
  iterations-1,L2_norm_residual);

  /* Prints matrix to txt file */
  matrix_to_txt(A, jmax, imax, m);

  return 0;
}

/*
 * Update boundary values from inner grid
 */
void updateBoundary(double *A, int jmax, int imax, int n, int m) {
  int i, j;
  for(j = 1; j <= imax; j++) {  /* Top and bottom part of the boundary */
    A[j] = A[m+j];
    A[(n-1)*m+j] = A[(n-2)*m+j];
  }
  for(i = 1; i<= jmax; i++) {   /* Left and right part of the boundary */
    A[i*m] = A[i*m+1];
    A[i*m+m-1] = A[i*m+m-2];
  }
}

/*
 * Prints the n by m matrix to a txt file
 */
void matrix_to_txt(double *matrix, int jmax, int imax, int m) {
  FILE *f = fopen("output.txt", "w");
  if (f == NULL) {
      printf("Error opening file!\n");
      exit(1);
  }
  int i, j;
  for(j = 1; j <= jmax; j++) { /* Row */
    for (i = 1; i <= imax; i++) { /* Column */
      fprintf(f, "%f ", matrix[j*m + i]);
    }
    fprintf(f, "\n");
  }
  fclose(f);
}

/* 
 * Print given matrix that is n x m
 */
void print_matrix(double *matrix, int n, int m) {
  int i, j;
  for(i=0; i<n; i++) {
    for(j=0; j<m; j++) {
      printf("%f ", matrix[i*m+j]);
    }
    printf("\n");
  }
}
