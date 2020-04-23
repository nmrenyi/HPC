#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#define LEN 5

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_size);

    MPI_Group world_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);

    int n = 3;
    const int ranks[3] = {1,3,5};
    const int ori1[1] = {1};
    const int ori2[1] = {0};
    int root1, root2;

    // 从world_group进程组中构造出来两个进程组
    MPI_Group group1, group2;
    MPI_Group_incl(world_group, n, ranks, &group1);
    MPI_Group_excl(world_group, n, ranks, &group2);
    // 根据group1 group2分别构造两个通信域
    MPI_Comm comm1, comm2;
    MPI_Comm_create(MPI_COMM_WORLD, group1, &comm1);
    MPI_Comm_create(MPI_COMM_WORLD, group2, &comm2);

    // 维护发送缓冲区和接受缓冲区
    int i;
    double *sbuf, *rbuf1, *rbuf2, *rbuf3;
    sbuf = malloc(LEN*sizeof(double));
    rbuf1 = malloc(LEN*sizeof(double));
    rbuf2 = malloc(LEN*sizeof(double));
    rbuf3 = malloc(LEN*sizeof(double));
    srand(world_rank*100);
    for(i=0; i<LEN; i++) sbuf[i] = (1.0*rand()) / RAND_MAX;
    fprintf(stderr,"rank %d:\t", world_rank);
    for(i=0; i<LEN; i++) fprintf(stderr,"%f\t",sbuf[i]);
    fprintf(stderr,"\n");
    MPI_Group_translate_ranks(world_group, 1, ori1, group1, &root1);
    MPI_Group_translate_ranks(world_group, 1, ori2, group2, &root2);
    // MPI_COMM_WORLD comm1 comm2分别执行不同的归约操作
    if (MPI_COMM_NULL!=comm1) { // comm1
        MPI_Reduce(sbuf, rbuf1, LEN, MPI_DOUBLE, MPI_MAX, root1, comm1);
        int rank_1;
        MPI_Comm_rank(comm1, &rank_1);
        if (root1==rank_1) {
            fprintf(stderr,"MAX:\t");
            for(i=0; i<LEN; i++) fprintf(stderr,"%f\t",rbuf1[i]);
            fprintf(stderr,"\n");
        }
    } 
    else if (MPI_COMM_NULL!=comm2) { // comm2
        MPI_Reduce(sbuf, rbuf2, LEN, MPI_DOUBLE, MPI_MIN, root2, comm2);
        int rank_2;
        MPI_Comm_rank(comm2, &rank_2);
        if (root2==rank_2) {
            fprintf(stderr,"MIN:\t");
            for(i=0; i<LEN; i++) fprintf(stderr,"%f\t",rbuf2[i]);
            fprintf(stderr,"\n");
        }
    }
    MPI_Reduce(sbuf, rbuf3, LEN, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // MPI_COMM_WORLD 
    if (0==world_rank) {
        fprintf(stderr,"SUM:\t");
        for(i=0; i<LEN; i++) fprintf(stderr,"%f\t",rbuf3[i]);
        fprintf(stderr,"\n");
    }
    // 清理进程组和通信域
    if(MPI_GROUP_NULL!=group1) MPI_Group_free(&group1);
    if(MPI_GROUP_NULL!=group2) MPI_Group_free(&group2);
    if(MPI_COMM_NULL!=comm1) MPI_Comm_free(&comm1);
    if(MPI_COMM_NULL!=comm2) MPI_Comm_free(&comm2);
    MPI_Finalize();
}