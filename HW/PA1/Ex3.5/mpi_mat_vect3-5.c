/*
 * File:     prog3.5_mpi_mat_vect_col.c
 *
 * Purpose:  Implement matrix vector multiplication when the matrix
 *           has a block column distribution
 *
 * Compile:  mpicc -g -Wall -o prog3.5_mpi_mat_vect_col prog3.5_mpi_mat_vect_col.c
 * Run:      mpiexec -n <comm_sz> ./prog3.5_mpi_mat_vect_col
 *
 * Input:    order of matrix, matrix, vector
 * Output:   product of matrix and vector.  If DEBUG is defined, the 
 *           order, the input matrix and the input vector
 *
 * Notes:
 * 1.  The matrix should be square and its order should be evenly divisible
 *     by comm_sz
 * 2.  The program stores the local matrices as one-dimensional arrays
 *     in row-major order
 * 3.  The program uses a derived datatype for matrix input and output
 *
 * Author:   Jinyoung Choi
 *
 * IPP:      Programming Assignment 3.5
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>


void Check_for_error(int local_ok, char fname[], char message[], 
      MPI_Comm comm);
void Get_dims(int* m_p, int* local_m_p, int* n_p, int* local_n_p,
      int my_rank, int comm_sz, MPI_Comm comm, double* distribution_time);
void Allocate_arrays(double**, double**, double**, double**, double** local_A_pp, double** local_x_pp, 
      double** local_y_pp, int m, int n, int local_m, int local_n, 
      MPI_Comm comm, double* allocation_time);
void Build_derived_type(int m, int local_m, int n, int local_n,
      MPI_Datatype* block_col_mpi_t_p);
void Read_matrix(char prompt[], double* global_A, double local_A[], int m, int local_n, 
      int n, MPI_Datatype block_col_mpi_t, int my_rank, MPI_Comm comm, double* distribution_time);
void Print_matrix(char title[], double local_A[], int m, int local_n, 
      int n, MPI_Datatype block_col_mpi_t, int my_rank, MPI_Comm comm);
void Read_vector(char prompt[], double* global_x, double local_vec[], int n, int local_n, 
      int my_rank, MPI_Comm comm, double* distribution_time);
void Print_vector(char title[],double*, double local_vec[], int n,
      int local_n, int my_rank, MPI_Comm comm);
void Gather_result(double*, double local_vec[], int n,
      int local_n, int my_rank, MPI_Comm comm);
void Mat_vect_mult(double parallel_result[], double local_A[], double local_x[], 
      double local_y[], int local_m, int m, int n, int local_n, 
      int comm_sz, MPI_Comm comm, double* parallel_calc_time, double* distribution_time);
void Serial_mat_vect_mult(double*, double* global_A, double* global_x, int n);
double calcDifference(double* parallel_result, double* serial_result, int n);
/*-------------------------------------------------------------------*/

// const int ORDER = 20000;
int main(void) {
   double* global_A = NULL;
   double* global_x = NULL;
   double* parallel_result = NULL;
   double* serial_result = NULL;
   double* local_A;
   double* local_x;
   double* local_y;
   
   int m,n;
   int local_m, local_n;
   
   double allocation_time = 0.0;
   double distribution_time = 0.0;
   double parallel_calc_time = 0.0;

   int my_rank, comm_sz;
   MPI_Comm comm;
   MPI_Datatype block_col_mpi_t;
   
   MPI_Init(NULL, NULL);
   comm = MPI_COMM_WORLD;

   MPI_Comm_size(comm, &comm_sz);
   MPI_Comm_rank(comm, &my_rank);
   
   Get_dims(&m, &local_m, &n, &local_n, my_rank, comm_sz, comm, &distribution_time);
   Allocate_arrays(&parallel_result, &serial_result, &global_A, &global_x, &local_A, &local_x, &local_y, m, n, local_m, local_n, comm, &allocation_time);
   Build_derived_type(m, local_m, n, local_n, &block_col_mpi_t);
   Read_matrix("A", global_A, local_A, m, local_n, n, block_col_mpi_t, my_rank, comm, &distribution_time);

#  ifdef DEBUG
   Print_matrix("A", local_A, m, local_n, n, block_col_mpi_t, my_rank, comm);
#  endif
   
   Read_vector("x", global_x, local_x, n, local_n, my_rank, comm, &distribution_time);

   MPI_Barrier(comm);
   if (my_rank == 0) {
      printf("Calculating parallel...\n");
   }

   Mat_vect_mult(parallel_result, local_A, local_x, local_y, local_m, m, n, local_n, comm_sz, comm, &parallel_calc_time, &distribution_time);   
   if (my_rank == 0) {
      printf("parallel calculation complete, time cost = %f s\n", parallel_calc_time);
   }

   double parallel_distribution_time = 0.0;
   MPI_Reduce(&distribution_time, &parallel_distribution_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

   if (my_rank == 0) {
      printf("Calculating serial...\n");
      double serial_calc_start = MPI_Wtime();
      Serial_mat_vect_mult(serial_result, global_A, global_x, n);
      double serial_calc_end = MPI_Wtime();
      double serial_calc_time = serial_calc_end - serial_calc_start;
      printf("serial calculation complete, time cost = %f s\n", serial_calc_time);
      printf("\n\n");
      printf("==================Calculation Complete================\n");
      printf("========================REPORT========================\n");
      printf("Ex3.5\n");
      printf("Process number = %d, matrix order = %d\n", comm_sz, n);
      printf("L2 norm = %f\n", calcDifference(parallel_result, serial_result, n));
      printf("Parallel calc time = %f s, Serial calc time = %f s, SpeedUp Ratio = %f\n", parallel_calc_time, serial_calc_time, serial_calc_time / parallel_calc_time);
      printf("Parallel distribution time = %f s, parallel all time = %f, overall SpeedUp Ratio = %f\n", parallel_distribution_time, parallel_calc_time + parallel_distribution_time, serial_calc_time / (parallel_calc_time +parallel_distribution_time));
      printf("allocation time = %f\n", allocation_time);
      printf("======================REPORT END====================\n\n");
      free(global_A);
      free(global_x);
      free(parallel_result);
      free(serial_result);
   }

   free(local_A);
   free(local_x);
   free(local_y);
   MPI_Type_free(&block_col_mpi_t);
   MPI_Finalize();
   return 0;
}  /* main */

double calcDifference(double* parallel_result, double* serial_result, int n) {
   double result = 0.0;
   int i = 0;
#  ifdef DEBUG
   printf("parallel result is:\n");
   for (i = 0; i < n; i++)
      printf("%f ", parallel_result[i]);
   printf("\n");
   printf("serial result is:\n");
   for (i = 0; i < n; i++)
      printf("%f ", serial_result[i]);
   printf("\n");
   #endif
   for (i = 0; i < n; i++) {
      double diff = parallel_result[i] - serial_result[i];
      result += diff * diff;
   }
   result = sqrt(result);
   return result;
}

/*-------------------------------------------------------------------*/
void Check_for_error(
   int       local_ok   /* in */, 
   char      fname[]    /* in */,
   char      message[]  /* in */, 
   MPI_Comm  comm       /* in */) {
   
   int ok;
   
   MPI_Allreduce(&local_ok, &ok, 1, MPI_INT, MPI_MIN, comm);
   if (ok == 0) {
      int my_rank;
      MPI_Comm_rank(comm, &my_rank);
      if (my_rank == 0) {
         fprintf(stderr, "Proc %d > In %s, %s\n", my_rank, fname, 
               message);
         fflush(stderr);
      }
      MPI_Finalize();
      exit(-1);
   }
}  /* Check_for_error */


/*-------------------------------------------------------------------*/
void Get_dims(
      int*      m_p        /* out */, 
      int*      local_m_p  /* out */,
      int*      n_p        /* out */,
      int*      local_n_p  /* out */,
      int       my_rank    /* in  */,
      int       comm_sz    /* in  */,
      MPI_Comm  comm       /* in  */,
      double* distribution_time
      ) {
   
   int local_ok = 1;
   if (my_rank == 0) {
      printf("Enter the order of the matrix\n");
      scanf("%d", m_p);
   }
   MPI_Barrier(comm);
   double start_time = MPI_Wtime();
   MPI_Bcast(m_p, 1, MPI_INT, 0, comm);
   *distribution_time += MPI_Wtime() - start_time;

   *n_p = *m_p;
   // *n_p = ORDER;
   // *m_p = ORDER;
   if (*m_p <= 0 || *m_p % comm_sz != 0) local_ok = 0;
   Check_for_error(local_ok, "Get_dims",
               "m and n must be positive and evenly divisible by comm_sz", 
               comm);
   
   *local_m_p = *m_p/comm_sz;
   *local_n_p = *n_p/comm_sz;
}  /* Get_dims */


/*-------------------------------------------------------------------*/
void Allocate_arrays(
   double**  parallel_res_pp      ,
   double**  serial_res_pp        ,
   double**  global_A_pp          ,
   double**  global_x_pp          ,
   double**  local_A_pp  /* out */, 
   double**  local_x_pp  /* out */, 
   double**  local_y_pp  /* out */, 
   int       m           /* in  */,   
   int       n                    ,
   int       local_m     /* in  */, 
   int       local_n     /* in  */, 
   MPI_Comm  comm        /* in  */,
   double* allocation_time) {
   

   int local_ok = 1;
   MPI_Barrier(comm);
   double local_start_time = MPI_Wtime();
   *parallel_res_pp = malloc(n * sizeof(double));
   *serial_res_pp = malloc( n * sizeof(double));
   *global_A_pp = malloc(m * n * sizeof(double));
   *global_x_pp = malloc(n * sizeof(double));
   *local_A_pp = malloc(m*local_n*sizeof(double));
   *local_x_pp = malloc(local_n*sizeof(double));
   *local_y_pp = malloc(local_m*sizeof(double));
   double local_allocation_time = MPI_Wtime() - local_start_time;
   MPI_Reduce(&local_allocation_time, allocation_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
   if (*local_A_pp == NULL || local_x_pp == NULL ||
      local_y_pp == NULL) local_ok = 0;
   Check_for_error(local_ok, "Allocate_arrays",
               "Can't allocate local arrays", comm);
}  /* Allocate_arrays */


/*-------------------------------------------------------------------*/
void Build_derived_type(int m, int local_m, int n, int local_n,
      MPI_Datatype* block_col_mpi_t_p) {
   MPI_Datatype vect_mpi_t;

   /* m blocks each containing local_n elements */
   /* The start of each block is n doubles beyond the preceding block */
   MPI_Type_vector(m, local_n, n, MPI_DOUBLE, &vect_mpi_t);
   
   /* Resize the new type so that it has the extent of local_n doubles */
   MPI_Type_create_resized(vect_mpi_t, 0, local_n*sizeof(double),
         block_col_mpi_t_p);
   MPI_Type_commit(block_col_mpi_t_p);
}  /* Build_derived_type */


/*-------------------------------------------------------------------*/
void Read_matrix(
             char          prompt[]         /* in  */, 
             double        A[]                       ,
             double        local_A[]        /* out */, 
             int           m                /* in  */, 
             int           local_n          /* in  */, 
             int           n                /* in  */,
             MPI_Datatype  block_col_mpi_t  /* in  */,
             int           my_rank          /* in  */,
             MPI_Comm      comm             /* in  */, 
             double* distribution_time) {
   int local_ok = 1;
   int i, j;

   if (my_rank == 0) {
      if (A == NULL) local_ok = 0;
      Check_for_error(local_ok, "Read_matrix",
                  "Can't allocate temporary matrix", comm);

      // printf("Enter the matrix %s\n", prompt);
      for (i = 0; i < m; i++)
         for (j = 0; j < n; j++)
            A[i * n + j] = -1 + rand() * 1.0 / RAND_MAX * 2;
   }
   else {
      Check_for_error(local_ok, "Read_matrix", "Can't allocate temporary matrix", comm);
   }
   double local_start = MPI_Wtime();
   MPI_Scatter(A, 1, block_col_mpi_t, local_A, m*local_n, MPI_DOUBLE, 0, comm);
   *distribution_time += MPI_Wtime() - local_start;

}  /* Read_matrix */


/*-------------------------------------------------------------------*/
void Print_matrix(char title[], double local_A[], int m, int local_n, 
      int n, MPI_Datatype block_col_mpi_t, int my_rank, MPI_Comm comm) {
   double* A = NULL;
   int local_ok = 1;
   int i, j;

   if (my_rank == 0) {
      A = malloc(m*n*sizeof(double));
      if (A == NULL) local_ok = 0;
      Check_for_error(local_ok, "Print_matrix",
                  "Can't allocate temporary matrix", comm);

      MPI_Gather(local_A, m*local_n, MPI_DOUBLE, A, 1, block_col_mpi_t,
            0, comm);

      printf("The matrix %s\n", title);
      for (i = 0; i < m; i++) {
         for (j = 0; j < n; j++)
            printf("%.2f ", A[i*n+j]);
         printf("\n");
      }
      
      free(A);
   } else {
      Check_for_error(local_ok, "Print_matrix",
                  "Can't allocate temporary matrix", comm);
      MPI_Gather(local_A, m*local_n, MPI_DOUBLE, A, 1, block_col_mpi_t,
            0, comm);
   }
}  /* Print_matrix */


/*-------------------------------------------------------------------*/
void Read_vector(
             char      prompt[]     /* in  */, 
             double    vec[]                   ,
             double    local_vec[]  /* out */, 
             int       n            /* in  */,
             int       local_n      /* in  */,
             int       my_rank      /* in  */,
             MPI_Comm  comm         /* in  */, 
             double* distribution_time) {
   int i, local_ok = 1;
   
   if (my_rank == 0) {
      if (vec == NULL) local_ok = 0;
      Check_for_error(local_ok, "Read_vector",
                  "Can't allocate temporary vector", comm);
      // printf("Enter the vector %s\n", prompt);
      for (i = 0; i < n; i++)
         vec[i] = -1 + rand() * 1.0 / RAND_MAX * 2;
   } else {
      Check_for_error(local_ok, "Read_vector",
                  "Can't allocate temporary vector", comm);
   }
   double local_start = MPI_Wtime();
   MPI_Scatter(vec, local_n, MPI_DOUBLE, local_vec, local_n, MPI_DOUBLE, 0, comm);
   *distribution_time += MPI_Wtime() - local_start;
}  /* Read_vector */


/*-------------------------------------------------------------------*/
void Mat_vect_mult(
               double    parallel_result[]   ,
               double    local_A[]  /* in  */, 
               double    local_x[]  /* in  */, 
               double    local_y[]  /* out */,
               int       local_m    /* in  */, 
               int       m          /* in  */,
               int       n         /* in  */,
               int       local_n    /* in  */, 
               int       comm_sz,
               MPI_Comm  comm       /* in  */, 
               double* parallel_calc_time, 
               double* distribution_time) {
   
   double* my_y;
   int i, loc_j;
   int local_ok = 1;
   
   my_y = malloc(n*sizeof(double));
   if (my_y == NULL) local_ok = 0;
   Check_for_error(local_ok, "Mat_vect_mult",
               "Can't allocate temporary arrays", comm);
   MPI_Barrier(comm);
   double local_calc_start = MPI_Wtime();
   for (i = 0; i < m ; i++) {
      my_y[i] = 0.0;
      for (loc_j = 0; loc_j < local_n ; loc_j++)
         my_y[i] += local_A[i*local_n + loc_j]*local_x[loc_j];
   }
   double local_calc_time = MPI_Wtime() - local_calc_start;
   MPI_Reduce(&local_calc_time, parallel_calc_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
   
   double local_reduce_start = MPI_Wtime();
   MPI_Reduce(my_y, parallel_result, n, MPI_DOUBLE, MPI_SUM, 0, comm);
   *distribution_time += MPI_Wtime() - local_reduce_start;

   free(my_y);
}  /* Mat_vect_mult */

void Serial_mat_vect_mult(double* result, double* A, double* x, int n) {
   int i = 0, j = 0;
   for (i = 0; i < n; i++) {
      result[i] = 0;
   }
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         result[i] += A[i * n + j] * x[j];
      }
   }
   // printf("The serial vector is:\n");
   // for (i = 0; i < n; i++) {
   //    printf("%f ", result[i]);
   // }
   // printf("\n");
}

/*-------------------------------------------------------------------*/
void Print_vector(
              char      title[]     /* in */, 
              double*   vec                 ,
              double    local_vec[] /* in */, 
              int       n           /* in */,
              int       local_n     /* in */,
              int       my_rank     /* in */,
              MPI_Comm  comm        /* in */) {
   int local_ok = 1;
   
   if (my_rank == 0) {
      if (vec == NULL) local_ok = 0;
      Check_for_error(local_ok, "Print_vector",
                  "Can't allocate temporary vector", comm);
      MPI_Gather(local_vec, local_n, MPI_DOUBLE,
               vec, local_n, MPI_DOUBLE, 0, comm);
      printf("\nThe vector %s\n", title);
      int i = 0;
      for (i = 0; i < n; i++)
         printf("%f ", vec[i]);
      printf("\n");
   }  else {
      Check_for_error(local_ok, "Print_vector",
                  "Can't allocate temporary vector", comm);
      
      MPI_Gather(local_vec, local_n, MPI_DOUBLE,
               vec, local_n, MPI_DOUBLE, 0, comm);
   }
}  /* Print_vector */

void Gather_result(double* vec, double local_vec[], int n,
      int local_n, int my_rank, MPI_Comm comm) {
   int local_ok = 1;
   
   if (my_rank == 0) {
      if (vec == NULL) local_ok = 0;
      Check_for_error(local_ok, "Gather_vector",
                  "NULL ptr here for gather", comm);
      MPI_Gather(local_vec, local_n, MPI_DOUBLE,
               vec, local_n, MPI_DOUBLE, 0, comm);

   }  else {
      Check_for_error(local_ok, "Gather_vector",
                  "Can't allocate temporary vector", comm);
      
      MPI_Gather(local_vec, local_n, MPI_DOUBLE,
               vec, local_n, MPI_DOUBLE, 0, comm);
   }
      
}