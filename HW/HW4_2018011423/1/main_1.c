
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>  /* Semaphores are not part of Pthreads */

const int MAX_THREADS = 1024;
const int MSG_MAX = 100;

/* Global variables:  accessible to all threads */

char* message;
int message_available;
pthread_mutex_t mutex;

int thread_count;


void Usage(char* prog_name);
void *Send_msg(void* rank);  /* Thread function */

/*--------------------------------------------------------------------*/
int main(int argc, char* argv[]) {
   long       thread;
   pthread_t* thread_handles; 
   message_available = 0;
   thread_count = 2;
   message = malloc(MSG_MAX * sizeof(char));
   pthread_mutex_init(&mutex, NULL);

   thread_handles = malloc (thread_count*sizeof(pthread_t));

   for (thread = 0; thread < thread_count; thread++)
      pthread_create(&thread_handles[thread], (pthread_attr_t*) NULL,
          Send_msg, (void*) thread);

   for (thread = 0; thread < thread_count; thread++) {
      pthread_join(thread_handles[thread], NULL);
   }

   pthread_mutex_destroy(&mutex);

   free(thread_handles);
   free(message);
   return 0;
}  /* main */


/*--------------------------------------------------------------------
 * Function:    Usage
 * Purpose:     Print command line for function and terminate
 * In arg:      prog_name
 */
void Usage(char* prog_name) {

   fprintf(stderr, "usage: %s <number of threads>\n", prog_name);
   exit(0);
}  /* Usage */


/*-------------------------------------------------------------------
 * Function:       Send_msg
 * Purpose:        Create a message and ``send'' it by copying it
 *                 into the global messages array.  Receive a message
 *                 and print it.
 * In arg:         rank
 * Global in:      thread_count
 * Global in/out:  messages, semaphores
 * Return val:     Ignored
 * Note:           The my_msg buffer is freed in main
 */

void *Send_msg(void* rank) {
   long my_rank = (long) rank;
   while (1) {
      pthread_mutex_lock(&mutex);
      if (my_rank == 0) {
         if (message_available) {
            printf("%s\n", message);
            pthread_mutex_unlock(&mutex);
            break;
         }
      } else {
         sprintf(message, "Hello from producer");
         message_available = 1;
         pthread_mutex_unlock(&mutex);
         break;
      }
      pthread_mutex_unlock(&mutex);
   }
   return NULL;
}  /* Send_msg */
