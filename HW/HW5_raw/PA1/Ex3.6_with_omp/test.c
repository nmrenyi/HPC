#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

// const int ORDER = 20000;
int main(void) {


   int my_rank, comm_sz;
   MPI_Comm comm;
   double local[2] = {0, 0};
    double global[6] = {1,2,3,4,5,6};
    int displacements[6] = {0, 2, 4, 0, 2, 4};
    int data_count[6] = {2,2,2,2,2,2};


    MPI_Init(NULL, NULL);
    comm = MPI_COMM_WORLD;

   MPI_Comm_size(comm, &comm_sz);
   MPI_Comm_rank(comm, &my_rank);
    MPI_Scatterv(global, data_count, displacements, MPI_DOUBLE, local, data_count[my_rank],MPI_DOUBLE,0, comm);
   int i = 0;
   for (i = 0; i < 2; i++) {
       printf("rank = %d, data[%d] = %f\n", my_rank, i, local[i]);
   }

   MPI_Finalize();
   return 0;
}  /* main */