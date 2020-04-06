#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include "timer.h"

const int MAX_THREADS = 1024;
int done = 0;  //  finished threads count
long thread_count;  // total threads count
long long n;  //  fib(n)

pthread_mutex_t mutex;  // mutex used in thread writing
pthread_mutex_t done_mutex;  // mutex used in finished threads counting
pthread_mutex_t main_mutex;  // mutex used for main_cond
pthread_cond_t main_cond;  // conditional variable used to block the main thread

pthread_barrier_t barrier_p;  // barrier used in timing
double parallel_time = 0.0;  // parallel calculation time
int* calc_count;  // how many matrix multiplication a thread need to do
long long matrix[4] = {1, 0, 0, 1};  // global matrix, initialized as an identity matrix

void* Thread_result(void* rank);  // thread function
long long Serial_fib(long long n);  /* Only executed by main thread */

int main(int argc, char* argv[]) {
   long       thread;  /* Use long in case of a 64-bit system */
   pthread_t* thread_handles;
   double start, finish, elapsed;

   printf("Please put in the number of threads\n");
   scanf("%ld", &thread_count);
   if (thread_count > MAX_THREADS) {
       thread_count = MAX_THREADS;
       printf("thread count has exceeded the max thread limit %d, thread count has been set as %d\n", MAX_THREADS, MAX_THREADS);
   }
   // printf("threads = %ld\n", thread_count);
   printf("Please put in which fib do you want to calculate, for example put in k for fib(k)\n");
   scanf("%lld", &n);
   // printf("n = %lld\n", n);
   calc_count = malloc(thread_count * sizeof(int));
   thread_handles = (pthread_t*) malloc (thread_count*sizeof(pthread_t)); 
   pthread_mutex_init(&mutex, NULL);
   pthread_mutex_init(&main_mutex, NULL);
   pthread_mutex_init(&done_mutex, NULL);
   pthread_cond_init(&main_cond, NULL);
   pthread_barrier_init(&barrier_p, NULL, thread_count);

   // calculating the matrix multiplication times for each thread
   int quotient = (n - 1) / thread_count;
   int remainder = (n - 1) % thread_count;
   int i = 0;
   for (i = 0; i < thread_count; i++) {
      calc_count[i] = quotient;
   }
   for (i = 0; i < remainder; i++) {
      calc_count[i] += 1;
   }
   long long paralell_result;


   if (n == 0) {
      paralell_result = 0;
   } else {
      for (thread = 0; thread < thread_count; thread++)  
         pthread_create(&thread_handles[thread], NULL,
            Thread_result, (void*)thread);  

      // wait for the threads to finish
      pthread_mutex_lock( &main_mutex );
      while (done < thread_count) {
         if (done < thread_count)
            while(pthread_cond_wait( &main_cond, &main_mutex ) != 0); 
         printf("done = %d / %ld\n", done, thread_count);
      }
      paralell_result = matrix[0];
      pthread_mutex_unlock( &main_mutex );
      for (thread = 0; thread < thread_count; thread++) 
         pthread_join(thread_handles[thread], NULL); 
   }

   printf("Parallel fib(%lld) = %lld\n", n, paralell_result);
   printf("The elapsed time is %e seconds\n", parallel_time);

   GET_TIME(start);
   long long serial_result = Serial_fib(n);
   GET_TIME(finish);
   elapsed = finish - start;
   printf("Serial  fib(%lld) = %lld\n", n, serial_result);
   printf("The elapsed time is %e seconds\n", elapsed);
   int res = 0;

   printf("Threads = %ld, n = %lld\n", thread_count, n);
   if (serial_result == paralell_result) {
      printf("Parallel == Serial, OK!\n");
      res = 0;
   } else {
      printf("Parallel != Serial, Wrong!\n");
      res = 1;
   }
   pthread_mutex_destroy(&mutex);
   pthread_mutex_destroy(&main_mutex);
   pthread_mutex_destroy(&done_mutex);
   pthread_cond_destroy(&main_cond);
   pthread_barrier_destroy(&barrier_p);

   free(thread_handles);
   free(calc_count);

   // right or wrong can be judged from the return value
   // if right, return 0, if wrong return 1
   return res;
}  /* main */

// thread function
void* Thread_result(void* rank) {
   pthread_barrier_wait(&barrier_p);
   double start = 0.0;
   GET_TIME(start);
   long my_rank = (long) rank;
   if (calc_count[my_rank] != 0) {
      long long local_matrix[4] = {1, 0, 0, 1};
      long long i = 0;
      // calculate local result
      for (i = 0; i < calc_count[my_rank]; i++) {
         long long a11 = local_matrix[0];
         long long a12 = local_matrix[1];
         long long a21 = local_matrix[2];
         long long a22 = local_matrix[3];
         local_matrix[0] = a11 + a12;
         local_matrix[1] = a11;
         local_matrix[2] = a21 + a22;
         local_matrix[3] = a21;
      }

      // calculating global result
      pthread_mutex_lock(&mutex);
      long long a11 = matrix[0];
      long long a12 = matrix[1];
      long long a21 = matrix[2];
      long long a22 = matrix[3];
      matrix[0] = a11 * local_matrix[0] + a12 * local_matrix[2];
      matrix[1] = a11 * local_matrix[1] + a12 * local_matrix[3];
      matrix[2] = a21 * local_matrix[0] + a22 * local_matrix[2];
      matrix[3] = a21 * local_matrix[1] + a22 * local_matrix[3];
      pthread_mutex_unlock(&mutex);
   }
   double finish = 0.0;
   GET_TIME(finish);
   double local_elapsed = finish - start;
   // calculating done and time
   pthread_mutex_lock(&done_mutex);
   done++;
   if (parallel_time < local_elapsed) {
      parallel_time = local_elapsed;
   }
   pthread_cond_signal(&main_cond); 
   pthread_mutex_unlock(&done_mutex);
   return NULL;
}  /* Thread_result */


// Serial calculation for fib(n)
long long Serial_fib(long long n) {
   if (n == 0) {
      return 0;
   }
   long long f = 1;
   long long g = 0;
   long long i = 0;
   for (i = 0; i < n - 1; i++) {
      f = f + g;
      g = f - g;
   }
   return f;
}  /* Serial_fib */
