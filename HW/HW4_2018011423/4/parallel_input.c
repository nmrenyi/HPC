#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <time.h>

struct list_node_s {
   int    data;
   struct list_node_s* next;
};

void Usage(char* prog_name);
int  Is_empty(struct list_node_s* head_p);
int  Insert(int value, struct list_node_s** head_p);
void Print(struct list_node_s* head_p);
int  Member(int value, struct list_node_s* head_p);
int  Delete(int value, struct list_node_s** head_p);
void Free_list(struct list_node_s** head_p);
int  Get_value(void);

int noWork = 0;
int thread_count = 0;
int* queue;
int head = 0, tail = 0;
int firstTime = 1;

struct list_node_s* head_p = NULL;  /* start with empty list */

pthread_mutex_t main_mutex;
pthread_cond_t main_cond;
pthread_mutex_t mutex;
pthread_cond_t cond;
int threadWorkOk = 0;

void* task(void* rank) {
    long my_rank = (long) rank;
    pthread_mutex_lock( &mutex );
    while (1) {
         //  in case when broadcast, some thread is still running
         if (noWork)
            break;
         if (firstTime)
            firstTime = 0;
         while(pthread_cond_wait(&cond, &mutex) != 0); 
      
        // when broadcast, the threads waiting above can break out
        if (noWork)
            break;
        // do some list work
        if (queue[head] == 0) {
            Insert(Get_value(), &head_p);
        } else if (queue[head] == 1) {
            Delete(Get_value(), &head_p);
        } else if (queue[head] == 2) {
            Member(Get_value(), head_p);
        } else if (queue[head] == 3) {
            Print(head_p);
        } else {
            printf("Undefined task\n");
        }

        // signal to the main thread so that the main thread can move on
        threadWorkOk = 1;
        pthread_cond_signal( &main_cond );
    }
    pthread_mutex_unlock( & mutex );
    return NULL;
}


int main(int argc, char* argv[]) {
    pthread_t* thread_handles;
    long   thread;  /* Use long in case of a 64-bit system */
    int n = 0;
    if(argc != 3) 
        Usage(argv[0]);

    thread_count = strtol(argv[1], NULL, 10);
    n = strtol(argv[2], NULL, 10);
    thread_handles = malloc(thread_count*sizeof(pthread_t));

    queue = malloc(n * sizeof(int));
    head = 0;
    tail = n;
    pthread_mutex_init(&mutex, NULL);
    pthread_cond_init(&cond, NULL);
    pthread_mutex_init(&main_mutex, NULL);
    pthread_cond_init(&main_cond, NULL);
    
    int i = 0;
    for (i = 0; i < n; i++) {
       scanf("%d", &queue[i]);
    }


    for (thread = 0; thread < thread_count; thread++) {
       pthread_create(&thread_handles[thread], NULL, task, (void*)thread);  
    }
    pthread_mutex_lock( &main_mutex );
    for (head = 0; head < tail; head++) {
       // for the first time, let the main thread wait for a while, waiting for the threads to arrive
       // in case the signal is sent but the threads haven't arrived at the wait so that the program runs into deadlock
        while(firstTime);
        pthread_cond_signal( &cond );  //  send signal to the threads, then the threads can move on to work
        // if the thread's signal has sent before the main thread wait, this "if" can avoid deadlock
        if (!threadWorkOk) {
            while (pthread_cond_wait(&main_cond, &main_mutex) != 0);  // wait for the thread to finish
        }
        threadWorkOk = 0;
        usleep(1000);  // sleep for 1ms
    }
    pthread_mutex_unlock( & main_mutex );
    
    // release all the threads
    noWork = 1;
    pthread_cond_broadcast( &cond );
    for (thread = 0; thread < thread_count; thread++) 
        pthread_join(thread_handles[thread], NULL); 
    Free_list(&head_p);

    // /* Destroy mutex and conditional variables*/
    pthread_mutex_destroy(&mutex);
    pthread_cond_destroy(&cond);
    pthread_mutex_destroy(&main_mutex);
    pthread_cond_destroy(&main_cond);
    free(queue);
    free(thread_handles);
    return 0;
} /* main */


int Insert(int value, struct list_node_s** head_pp) {
   struct list_node_s* curr_p = *head_pp;
   struct list_node_s* pred_p = NULL;
   struct list_node_s* temp_p;
   
   while (curr_p != NULL && curr_p->data < value) {
      pred_p = curr_p;
      curr_p = curr_p->next;
   }

   if (curr_p == NULL || curr_p->data > value) {
      temp_p = malloc(sizeof(struct list_node_s));
      temp_p->data = value;
      temp_p->next = curr_p;
      if (pred_p == NULL)
         *head_pp = temp_p;
      else
         pred_p->next = temp_p;
      return 1;
   } else { /* value in list */
      printf("%d is already in the list\n", value);
      return 0;
   }
}  /* Insert */

void Print(struct list_node_s* head_p) {
   struct list_node_s* curr_p;

   printf("list = ");

   curr_p = head_p;
   while (curr_p != (struct list_node_s*) NULL) {
      printf("%d ", curr_p->data);
      curr_p = curr_p->next;
   }
   printf("\n");
}  /* Print */

int  Member(int value, struct list_node_s* head_p) {
   struct list_node_s* curr_p;

   curr_p = head_p;
   while (curr_p != NULL && curr_p->data < value)
      curr_p = curr_p->next;

   if (curr_p == NULL || curr_p->data > value) {
      printf("%d is not in the list\n", value);
      return 0;
   } else {
      printf("%d is in the list\n", value);
      return 1;
   }
}  /* Member */

int Delete(int value, struct list_node_s** head_pp) {
   struct list_node_s* curr_p = *head_pp;
   struct list_node_s* pred_p = NULL;

   /* Find value */
   while (curr_p != NULL && curr_p->data < value) {
      pred_p = curr_p;
      curr_p = curr_p->next;
   }
   
   if (curr_p != NULL && curr_p->data == value) {
      if (pred_p == NULL) { /* first element in list */
         *head_pp = curr_p->next;
#        ifdef DEBUG
         printf("Freeing %d\n", value);
#        endif
         free(curr_p);
      } else { 
         pred_p->next = curr_p->next;
#        ifdef DEBUG
         printf("Freeing %d\n", value);
#        endif
         free(curr_p);
      }
      return 1;
   } else {
   printf("%d is not in the list\n", value);
      return 0;
   }
}  /* Delete */


void Free_list(struct list_node_s** head_pp) {
   struct list_node_s* curr_p;
   struct list_node_s* succ_p;

   if (Is_empty(*head_pp)) return;
   curr_p = *head_pp; 
   succ_p = curr_p->next;
   while (succ_p != NULL) {
#     ifdef DEBUG
      printf("Freeing %d\n", curr_p->data);
#     endif
      free(curr_p);
      curr_p = succ_p;
      succ_p = curr_p->next;
   }
#  ifdef DEBUG
   printf("Freeing %d\n", curr_p->data);
#  endif
   free(curr_p);
   *head_pp = NULL;
}  /* Free_list */
int  Is_empty(struct list_node_s* head_p) {
   if (head_p == NULL)
      return 1;
   else
      return 0;
}  /* Is_empty */

void Usage(char* prog_name) {

   fprintf(stderr, "usage: %s <number of threads> <number of tasks>\n", prog_name);
   exit(0);
}  /* Usage */
int  Get_value(void) {
   int val;

   printf("Please enter a value:  ");
   scanf("%d", &val);
   return val;
}  /* Get_value */
