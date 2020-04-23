
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>  /* Semaphores are not part of Pthreads */

const int MAX_THREADS = 1024;
const int MSG_MAX = 100;

/* Global variables:  accessible to all threads */

char* message;
int* message_available;
int* sent;
int* recved;

pthread_mutex_t mutex;

int thread_count;


void Usage(char* prog_name);
void *Send_msg(void* rank);  /* Thread function */

/*--------------------------------------------------------------------*/
int main(int argc, char* argv[]) {
   long       thread;
   pthread_t* thread_handles; 
   printf("Enter thread number\n");
   scanf("%d", &thread_count);
   message = malloc(MSG_MAX * sizeof(char));
   message_available = malloc(thread_count * (sizeof(int)));
   sent = malloc(thread_count * sizeof(int));
   recved = malloc(thread_count * sizeof(int));
   for (thread = 0; thread < thread_count; thread++) {
       message_available[thread] = 0;
       sent[thread] = 0;
       recved[thread] = 0;
   }

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
   free(sent);
   free(recved);
   free(message_available);
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

// indicates whether there're unreceived message in the buffer
int messageBusy = 0; 
void *Send_msg(void* rank) {
   long my_rank = (long) rank;
   while (1) {
      pthread_mutex_lock(&mutex);
      if (!messageBusy) {
         if (!sent[my_rank]) {
            sent[my_rank] = 1;
            sprintf(message, "Hello from %ld", my_rank);
            message_available[(my_rank + 1) % thread_count] = 1;
            messageBusy = 1;
            if (sent[my_rank] && recved[my_rank]) {
               pthread_mutex_unlock(&mutex);
               break;
            }
         }
      } else {
         if (message_available[my_rank]) {
            printf("%s, in %ld\n", message, my_rank);
            recved[my_rank] = 1;
            messageBusy = 0;
            if (sent[my_rank] && recved[my_rank]) {
               pthread_mutex_unlock(&mutex);
               break;
            }
         }
      }
      pthread_mutex_unlock(&mutex);
   }

   return NULL;
}  /* Send_msg */
