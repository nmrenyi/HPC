#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>

void CountSort(int a[], int n) {
    int i, j, count;
    int* temp = malloc(n * sizeof(int));
    for (i = 0; i < n; i++) {
        count = 0;
        for (j = 0; j < n; j++) {
            if (a[j] < a[i])
                count++;
            else if (a[j] == a[i] && j < i)
                count++;
        }
        temp[count] = a[i];
    }
    memcpy(a, temp, n * sizeof(int));
    free(temp);
}
void ParallelCountSort(int a[], int n, int thread_count, int chunk_size) {
    int i, j, count;
    int* temp = malloc(n * sizeof(int));
    # pragma omp parallel for num_threads(thread_count) schedule(static, chunk_size)\
        default(none) private(i, j, count) shared(a, n, temp, chunk_size)
    for (i = 0; i < n; i++) {
        count = 0;
        for (j = 0; j < n; j++) {
            if (a[j] < a[i])
                count++;
            else if (a[j] == a[i] && j < i)
                count++;
        }
        temp[count] = a[i];
    }

    memcpy(a, temp, n * sizeof(int));
    free(temp);
}
int compare (const void * a, const void * b) {
  return ( *(int*)a - *(int*)b );
}

int main() {
    int range = 100;
    int length = 1000000;
    int thread_count = 10;
    int chunkSize = 100;
    int* arr = malloc(length * sizeof(int));
    int* countSortArr = malloc(length * sizeof(int));
    int* countSortParallelArr = malloc(length * sizeof(int));
    int* qsortArr = malloc(length * sizeof(int));
    int i = 0;
    printf("Input length of the array\n");
    scanf("%d", &length);
    printf("Input thread count\n");
    scanf("%d", &thread_count);
    printf("Input Chunk Size\n");
    scanf("%d", &chunkSize);
    srand((unsigned)time(NULL));
    for (i = 0; i < length; i++) {
        arr[i] = rand() % range;
    }
    memcpy(countSortArr, arr, length * sizeof(int));
    memcpy(countSortParallelArr, arr, length * sizeof(int));
    memcpy(qsortArr, arr, length * sizeof(int));
    printf("Doing Serial Counting Sort...\n");
    double start = omp_get_wtime();
    CountSort(countSortArr, length);
    double CountSortTime = omp_get_wtime() - start;
    printf("Doing Parallel Counting Sort...\n");
    start = omp_get_wtime();
    ParallelCountSort(countSortParallelArr, length, thread_count, chunkSize);
    double CountSortParallelTime = omp_get_wtime() - start;
    printf("Doing QuickSort...\n");
    start = omp_get_wtime();
    qsort(qsortArr, length, sizeof(int), compare);
    double qsortTime = omp_get_wtime() - start;

    int countSortOK = 1, countSortParallelOK = 1;
    for (i = 0; i < length; i++) {
        if (countSortArr[i] != qsortArr[i]) {
            countSortOK = 0;
            printf("Serial Counting Sort Wrong at %d, countSort[%d] == %d, qsortArr[%d] == %d\n", i, i, countSortArr[i], i, qsortArr[i]);
            break;
        }
    }
    for (i = 0; i < length; i++) {
        if (countSortParallelArr[i] != qsortArr[i]) {
            countSortParallelOK = 0;
            printf("Parallel Counting Sort Wrong at %d, ParallelCountSort[%d] == %d, qsortArr[%d] == %d\n", i, i, countSortParallelArr[i], i, qsortArr[i]);
            break;
        }
    }

    if (countSortOK) {
        printf("Serial Counting Sort Correct!\n");
    }
    if (countSortParallelOK) {
        printf("Parallel Counting Sort Correct!\n");
    }
    printf("QuickSort Time : %f\n", qsortTime);
    printf("Serial Counting Sort Time : %f\n", CountSortTime);
    printf("Parallel Counting Sort Time : %f\n", CountSortParallelTime);
    free(arr);
    free(countSortArr);
    free(countSortParallelArr);
    free(qsortArr);
    if ((countSortOK + countSortParallelOK) != 2) {
        return 1;
    }
    return 0;
}
