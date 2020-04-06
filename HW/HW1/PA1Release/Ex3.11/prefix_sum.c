#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

const int SIZE = 10;


void ErrorMessage() {
    printf("mpiexec -n <p> prefix_sum <g|i>\n- p: the number of processes\n- g: generate random, distributed list\n- i: user will input list on process 0\n ");
}

void Get_args(int argc, char** argv, int* local_arr, int my_rank, int comm_sz, MPI_Comm comm) {
    int randomGenerate = 0;
    int error = 0;
    if (my_rank == 0) {
        if (argc != 2) {
            ErrorMessage();
            error = 1;
        } else {
            if (argv[1][1] == 'g') {
                randomGenerate = 1;
            } else if (argv[1][1] == 'i') {
                randomGenerate = 0;
            } else {
                ErrorMessage();
                error = 1;
            }
        }
    }

    MPI_Bcast(&error, 1, MPI_INT, 0, comm);
    if (error) {
        MPI_Finalize();
        exit(1);
    }

    MPI_Bcast(&randomGenerate, 1, MPI_INT, 0, comm);
    if (randomGenerate) {
        srandom((unsigned int)time(0));
        int i = 0;
        for (i = 0; i < SIZE; i++) {
            local_arr[i] = random() % 10;
        }
    } else {
        int* tmp_buff;
        if (my_rank == 0) {
            int i = 0;
            tmp_buff = malloc(comm_sz * SIZE * sizeof(int));
            for (i = 0; i < comm_sz * SIZE; i++) {
                scanf("%d", &tmp_buff[i]);
            }
        }
        MPI_Scatter(tmp_buff, SIZE, MPI_INT, local_arr, SIZE, MPI_INT, 0, comm);
        if (my_rank == 0) {
            free(tmp_buff);
        }
    }
}


int main(int argc, char* argv[]) {
    int comm_sz = 0, my_rank = 0;
    
    MPI_Comm comm;
    MPI_Init(NULL, NULL);
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &comm_sz);
    MPI_Comm_rank(comm, &my_rank);
    int local_arr[SIZE];  // local arr, with 10 elements
    int i = 0;

    // get input
    Get_args(argc, argv, local_arr, my_rank, comm_sz, comm);

    int local_sum = 0;  //  local sum
    for (i = 0; i < SIZE; i++) {
        local_sum += local_arr[i];
    }

    int process_prefix_sum = 0;  // the prefix sum of the local_sum
    MPI_Scan(&local_sum, &process_prefix_sum, 1, MPI_INT, MPI_SUM, comm);

    int prefix_sum[10];  //  the global prefix sum in each process, waiting to be combined to the main process
    for (i = SIZE - 1; i >= 0; i--) {
        prefix_sum[i] = process_prefix_sum;
        process_prefix_sum -= local_arr[i];
    }


    if (my_rank == 0) {
        int* global_sum = malloc(comm_sz * SIZE * sizeof(int));
        MPI_Gather(prefix_sum, SIZE, MPI_INT, global_sum, SIZE, MPI_INT, 0, comm);
        for (i = 0; i < SIZE * comm_sz; i++) {
            printf("%d ", global_sum[i]);
        }
        printf("\n");
        free(global_sum);
    } else {
        MPI_Gather(prefix_sum, SIZE, MPI_INT, NULL, SIZE, MPI_INT, 0, comm);
    }

    MPI_Finalize();
    
    return 0;
}
